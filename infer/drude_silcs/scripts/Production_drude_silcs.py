#!/usr/bin/python
###########################################
# Runs a MD simulation with Drude FF.
#
# Outputs DCD files, last frame (CRD) and RST files with positions and
# velocities. Requires a PSF, a CRD and RST from equilibration. It allows
# restarting a simulation at any given point (using both .chk and .rst files).
# It also makes a backup of restarted .log file in order to save previous saved
# energies.
#
# USAGE: python Production_drude.py -h
#
# mdpoleto@vt.edu -> version 2023
###########################################
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
from openmmplumed import *
import MDAnalysis as mda
import warnings
import parmed as pmd
import time

from sys import stdout, exit, stderr
import os, math, fnmatch
import argparse

warnings.filterwarnings("ignore", message="Found no information for attr:")
warnings.filterwarnings("ignore", message="Found missing chainIDs")
warnings.filterwarnings("ignore", message="Supplied AtomGroup was missing the following attributes")
warnings.filterwarnings("ignore", message="DCDReader currently makes independent timesteps")


# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-pdb', type=str, default=None, required=True,
                help='Input coordinate file (.pdb)')
ap.add_argument('-psf', type=str, default=None, required=True,
                help='Topology file in XPLOR format (.psf)')
ap.add_argument('-toppar', type=str, default='toppar.str', required=True,
                help='Force field stream file (ex. "toppar.str").')
ap.add_argument('-state', type=str, required=True,
                help='XML file to read positions/velocities from (.rst).')
ap.add_argument('-fc_silcs', default=50.208, type=float, required=True,
                help='Force constant for SILCS positional restraints (kJ/mol/nm^2). Default: 50.208')
ap.add_argument('-restraint_file', default="restraint.dat", type=str, required=True,
                help='File containing heavy atoms to restraint in SILCS')


# Options
ap.add_argument('-outname', type=str, default="output",
                help='Default name for output files. Default is "output".')

ap.add_argument('-runtime', default=100, type=float,
                help='Simulation length of each stride (in ps). Default 100.')
ap.add_argument('-nstride', default=5, type=int,
                help='Number of strides/chunks in which the simulation will be splitted. Default 5.')
ap.add_argument('-dt', default=1, type=int,
                help='Integration step (in fs). Default 1.')

ap.add_argument('-savefreq', type=float, default=5,
                help='Frequency (in ps) to save coordinates, checkpoints and trajectory.')
ap.add_argument('-printfreq', type=float, default=5,
                help='Frequency (in ps) to print and write in .log file.')

ap.add_argument('-firststride', type=int, default=1,
                help='First stride number. Default 1.')


ap.add_argument('-temp', default=298, type=float,
                help='Target temperature, in Kelvin. Default is 298.')
ap.add_argument('-pressure', default=1.0, type=float,
                help='Target pressure, in bar. Default is 1.0.')

ap.add_argument('-plumed', default=None, type=str,
                help='Plumed script to use Plumed with OpenMM')

#print(sys.argv)
cmd = ap.parse_args()

with open('Production_drude.dat', 'w') as out:
	out.write(' '.join(sys.argv[1:]))

#############################################
# Although we use checkpoints to restart simulations, an unexpected crash may
# harm the dcd integrity beyond repair. It is rare, but that may happenself.
# Therefore, we use 10 simulation strides over a loop, creating 10 dcd files that
# are concatenated afterwards.

jobname = cmd.outname

stride_time = cmd.runtime*picosecond		# in ps
dt = cmd.dt*femtosecond						#fs
nstride = cmd.nstride

print_freq  = cmd.printfreq*picosecond
savcrd_freq = cmd.savefreq*picosecond

temperature = cmd.temp*kelvin
pressure	= cmd.pressure*bar

fc_silcs = float(cmd.fc_silcs)
restraint_file = cmd.restraint_file

nsteps = int(stride_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nprint  = int(print_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nsavcrd = int(savcrd_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))

plumedscript = cmd.plumed

#############################################
# Defining functions to use below:
def backup_old_log(pattern, string):
	result = []
	for root, dirs, files in os.walk("./"):
		for name in files:
			if fnmatch.fnmatch(name, pattern):

				try:
					number = int(name[-2])
					avail = isinstance(number, int)
					#print(name,avail)
					if avail == True:
						result.append(number)
				except:
					pass

	if len(result) > 0:
		maxnumber = max(result)
	else:
		maxnumber = 0

	backup_file = "\#" + string + "." + str(maxnumber + 1) + "#"
	os.system("mv " + string + " " + backup_file)
	return backup_file

def get_cubic_box(psf, rstfile):

	f = open(rstfile, 'r')

	box = {}

	while True:
		line = f.readline()
		if not line: break

		if line.split()[0] == "<A":
			size = line.split()[1].strip('x="')
			box['A'] = float(size)
		elif line.split()[0] == "<B":
			size = line.split()[2].strip('y="')
			box['B'] = float(size)
		elif line.split()[0] == "<C":
			size = line.split()[3].strip('z="').strip('"/>')
			box['C'] = float(size)
		else:
			pass

	boxX = box['A']*nanometer
	boxY = box['B']*nanometer
	boxZ = box['C']*nanometer

	psf.setBox(boxX, boxY, boxZ)

	return psf

def read_toppar(filename):
	extlist = ['rtf', 'prm', 'str']

	parFiles = ()
	for line in open(filename, 'r'):
		if '!' in line: line = line.split('!')[0]
		parfile = line.strip()
		if len(parfile) != 0:
			ext = parfile.lower().split('.')[-1]
			if not ext in extlist: continue
			parFiles += ( parfile, )

	params = CharmmParameterSet( *parFiles )
	return params, parFiles

def restraints(system, crd, fc_silcs, restraint_file):

	boxlx = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
	boxly = system.getDefaultPeriodicBoxVectors()[1][1].value_in_unit(nanometers)
	boxlz = system.getDefaultPeriodicBoxVectors()[2][2].value_in_unit(nanometers)

	if fc_silcs > 0:
		# positional restraints for protein
		posresPROT = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2;')
		posresPROT.addPerParticleParameter('k')
		posresPROT.addPerParticleParameter('x0')
		posresPROT.addPerParticleParameter('y0')
		posresPROT.addPerParticleParameter('z0')
		for line in open(restraint_file, 'r'):
			segments = line.strip().split()
			atom1 = int(segments[0])
			state = segments[1]
			xpos  = crd.positions[atom1].value_in_unit(nanometers)[0]
			ypos  = crd.positions[atom1].value_in_unit(nanometers)[1]
			zpos  = crd.positions[atom1].value_in_unit(nanometers)[2]
			if state == 'SILCS' and fc_silcs > 0:
				fc_ppos = fc_silcs
				posresPROT.addParticle(atom1, [fc_ppos, xpos, ypos, zpos])
		system.addForce(posresPROT)

	return system
##############################################

#############################################
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tPDB file = " + str(cmd.pdb))
print("\tPSF file = " + str(cmd.psf))
print("\tToppar stream file = " + str(read_toppar(cmd.toppar)[1]))

print("\n\tSimulation_time = " + str(stride_time*nstride))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps*nstride))
print("\tNumber of strides = " + str(cmd.nstride) + " (" + str(stride_time) + " in each stride)")

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tSave checkpoint each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure))
#############################################

print("\n> Setting the system:\n")
print("\t- Reading force field directory...")
charmm_params = read_toppar(cmd.toppar)[0]

print("\t- Reading topology and structure file...")
psf = CharmmPsfFile(cmd.psf)
crd = PDBFile(cmd.pdb)

print("\t- Setting box (using information on -state file)...")
psf = get_cubic_box(psf, cmd.state)

print("\t- Creating system and setting parameters...")
system = psf.createSystem(charmm_params, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=HBonds)

nbforce = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), NonbondedForce)][0]
nbforce.setNonbondedMethod(NonbondedForce.PME)
nbforce.setEwaldErrorTolerance(0.0001)
nbforce.setCutoffDistance(1.2*nanometer)
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(1.0*nanometer)

# not every system has NBFIX terms, so check
cstnb = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), CustomNonbondedForce)]
if cstnb:
	nbfix = cstnb[0]
	nbfix.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
	nbfix.setCutoffDistance(1.2*nanometer)
	nbfix.setUseSwitchingFunction(True)
	nbfix.setSwitchingDistance(1.0*nanometer)

print("\t- Applying SILCS restraints... (using " + str(restraint_file) + ")")
system = restraints(system, crd, fc_silcs, restraint_file)

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

if plumedscript != None:
	print("\t- Loading Plumed script...")
	with open(plumedscript, 'r') as plumedfile:
		plumed_param = plumedfile.read()
	system.addForce(PlumedForce(plumed_param))
else:
	pass

print("\t- Setting integrator and thermostat...")
integrator = DrudeLangevinIntegrator(temperature, 5/picosecond, 1*kelvin, 20/picosecond, 0.001*picoseconds)
integrator.setMaxDrudeDistance(0.02) # Drude Hardwall


print('\t- Setting simulation context...')
simulation = Simulation(psf.topology, system, integrator)#, platform, properties)
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()

print('\t- Using platform:', simulation.context.getPlatform().getName())

# Opening a loop of extension NSTRIDE to simulate the entire STRIDE_TIME*NSTRIDE
for n in range(cmd.firststride, nstride + 1):

	print("\n\n>>> Simulating Stride #" + str(n) + " <<<")

	dcd_file = jobname + "_" + str(n) + ".dcd"
	chk_file = jobname + "_" + str(n) + ".chk"
	log_file = jobname + "_" + str(n) + ".log"
	rst_file = jobname + "_" + str(n) + ".rst"
	prv_rst_file = jobname + "_" + str(n-1) + ".rst"
	pdb_file = jobname + "_" + str(n) + ".pdb"

	if os.path.exists(rst_file):
		print("> Stride #" + str(n) + " finished (" + rst_file + " present). Moving to next stride... <")
		continue


	if os.path.exists(chk_file):
		simulation.loadCheckpoint(chk_file)

		print("> Restarting from checkpoint > " + chk_file + " <")

		# Calculating elapsed time
		chk_time = simulation.context.getState().getTime()
		chk_time_val = chk_time.value_in_unit(picosecond)
		chk_step = math.ceil(chk_time.in_units_of(picosecond)/dt.in_units_of(picosecond))
		print("\t- Elapsed simulation time = " + str(round(chk_time_val,2)*picosecond) + " (step = " + str(chk_step) + ")")

		# Calculating remaining running time
		remaining_time = (stride_time*n - chk_time).in_units_of(picosecond)
		remaining_nsteps = int(math.ceil(remaining_time.in_units_of(picosecond)/dt.in_units_of(picosecond)))

		if remaining_nsteps < 0:
			sys.exit("\n>>>> WARNING: checkpoint last step saved is beyond the last nstep requested. Maybe increase simulation time?<<<<<<\n")

		# Adjust remaining running time
		simulation.currentStep = chk_step
		simulation.context.setTime(chk_time.value_in_unit(picosecond))
		print("\t- Restarting from step " + str(chk_step) + " ...")

		dcd=DCDReporter(dcd_file, nsavcrd, append=True)
		print("\t- Appending to file " + jobname + ".dcd ...")

		backup_file = backup_old_log("*" + log_file + "*", log_file)
		print("\t- Backuping old log to " + backup_file + " ...")

		simulation.reporters.append(dcd)
		simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=(nsteps*nstride), remainingTime=True, separator='\t\t'))
		simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
		simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))

		print("\n> Simulating " + str(remaining_nsteps) + " steps... (Stride #" + str(n) + ")")
		simulation.step(remaining_nsteps)

	else:
		if n == cmd.firststride:
			equil_rst_file = cmd.state
			print("\n> Loading previous state from equilibration > " + equil_rst_file + " <")
			with open(equil_rst_file, 'r') as f:
				simulation.context.setState(XmlSerializer.deserialize(f.read()))
				currstep = int((n-1)*nsteps)
				currtime = currstep*dt.in_units_of(picosecond)
				simulation.currentStep = currstep
				simulation.context.setTime(currtime)
				print("> Current time: " + str(currtime) + " (Step = " + str(currstep) + ")")

		else:
			print("> Loading previous state from > " + prv_rst_file + " <")
			with open(prv_rst_file, 'r') as f:
				simulation.context.setState(XmlSerializer.deserialize(f.read()))
				currstep = int((n-1)*nsteps)
				currtime = currstep*dt.in_units_of(picosecond)
				simulation.currentStep = currstep
				simulation.context.setTime(currtime)
				print("> Current time: " + str(currtime) + " (Step = " + str(currstep) + ")")


		dcd = DCDReporter(dcd_file, nsavcrd)
		firstdcdstep = (currstep) + nsavcrd
		dcd._dcd = DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0

		simulation.reporters.append(dcd)
		simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=(nsteps*nstride), remainingTime=True, separator='\t\t'))
		simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
		simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))

		print("\n> Simulating " + str(nsteps) + " steps... (Stride #" + str(n) + ")")
		simulation.step(nsteps)

	simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


	##################################
	# Writing last frame information of stride
	print("\n> Writing stride state file (" + str(rst_file) + ")...")
	state = simulation.context.getState( getPositions=True, getVelocities=True )
	with open(rst_file, 'w') as f:
		f.write(XmlSerializer.serialize(state))

	last_frame = int(nsteps/nsavcrd)
	print("> Writing last coordinate (" + str(pdb_file) + ", frame = " + str(last_frame) + ")...")
	u = mda.Universe(cmd.psf, dcd_file)
	system = u.select_atoms('all')
	for ts in u.trajectory[(len(u.trajectory)-1):len(u.trajectory)]:
		with mda.Writer(str(pdb_file), system.n_atoms) as W:
			W.write(system)

try:
	quote = Quotes.getquote()
	print(quote)
except:
	pass

print("\n> Finished!\n")
