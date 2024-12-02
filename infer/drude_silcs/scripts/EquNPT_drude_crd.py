#!/usr/bin/python
###########################################
# Runs a NPT equilibration with Drude FF after a CHARMM36 NPT.
#
# Outputs DCD files, last frame (CRD) and RST files with positions and
# velocities. Requires a PSF, a CRD and a (CHARMM) RST. It allows restarting a
# simulation at any given point (using .chk). It also  makes a backup of
# restarted .log file in order to save previous saved energies.
#
# USAGE: python EquNPT_drude.py -h
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


# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-crd', type=str, default=None, required=True,
                help='Input coordinate file (.crd)')
ap.add_argument('-psf', type=str, default=None, required=True,
                help='Topology file in XPLOR format (.psf)')
ap.add_argument('-toppar', type=str, default='toppar.str', required=True,
                help='Force field stream file (ex. "toppar.str").')
ap.add_argument('-fc_bb', default=500, type=int,
                help='Force constant for backbone heavy atoms (kJ/mol/nm^2). Default: 500')
ap.add_argument('-fc_sc', default=500, type=int,
                help='Force constant for sidechain heavy atoms (kJ/mol/nm^2). Default: 500')
ap.add_argument('-restraint_file', default="restraint.dat", type=str, required=True,
                help='File containing heavy atoms to restraint separated by BB and SC.')
ap.add_argument('-state_pequil', type=str, required=False, default=None,
                help='XML file to read positions/velocities from previous equilibration (.rst).')
ap.add_argument('-state_c36', type=str, required=False,
                help='XML file to read box size from (.rst).')


# Options
ap.add_argument('-outname', type=str, default="output.drude.npt",
                help='Default name for output files. Default is "output".')

ap.add_argument('-runtime', default=100, type=float,
                help='Simulation length (in ps). Default = 100.')

ap.add_argument('-dt', default=1, type=int,
                help='Integration step (in fs). Default = 1.')
ap.add_argument('-savefreq', type=float, default=5,
                help='Frequency (in ps) to save coordinates, checkpoints and trajectory. Default = 5')
ap.add_argument('-printfreq', type=float, default=5,
                help='Frequency (in ps) to print and write in .log file. Default = 5')

ap.add_argument('-temp', default=298, type=float,
                help='Target temperature, in Kelvin. Default = 298.')
ap.add_argument('-pressure', default=1.0, type=float,
                help='Target pressure, in bar. Default is 1.0.')

ap.add_argument('-fc_dih', default=4, type=int,
                help='Force constant for dihedral restraint (kJ/mol/nm^2). Default: 4')
ap.add_argument('-dih_restraint_file', default=None, type=str,
                help='File containing dihedrals to restraint.')


cmd = ap.parse_args()

with open('Equilibration_NPT_drude.dat', 'w') as out:
	out.write(' '.join(sys.argv[1:]))

#############################################

jobname = cmd.outname

simulation_time = cmd.runtime*picosecond		# in ps
dt = cmd.dt*femtosecond						#fs

print_freq  = cmd.printfreq*picosecond
savcrd_freq = cmd.savefreq*picosecond

temperature = cmd.temp*kelvin
pressure	= cmd.pressure*bar

fc_bb = float(cmd.fc_bb)
fc_sc = float(cmd.fc_sc)
fc_dih = float(cmd.fc_dih)
restraint_file = cmd.restraint_file
dih_restraint_file = cmd.dih_restraint_file

nsteps = int(simulation_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nprint  = int(print_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nsavcrd = int(savcrd_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))

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

def restraints(system, crd, fc_bb, fc_sc, restraint_file):

	boxlx = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
	boxly = system.getDefaultPeriodicBoxVectors()[1][1].value_in_unit(nanometers)
	boxlz = system.getDefaultPeriodicBoxVectors()[2][2].value_in_unit(nanometers)

	if fc_bb > 0 or fc_sc > 0:
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
			if state == 'BB' and fc_bb > 0:
				fc_ppos = fc_bb
				posresPROT.addParticle(atom1, [fc_ppos, xpos, ypos, zpos])
			if state == 'SC' and fc_sc > 0:
				fc_ppos = fc_sc
				posresPROT.addParticle(atom1, [fc_ppos, xpos, ypos, zpos])
		system.addForce(posresPROT)

	return system

def dih_restraints(system, crd, fc_dih, restraint_file):

	# dihedral restraints
	dihres = CustomTorsionForce('fc_dih*max(0, abs(diff+wrap) - rwidth)^2; \
	                                 wrap = 2*pi*(step(-diff-pi)-step(diff-pi)); \
	                                 diff = theta - rtheta0; \
	                                 rtheta0 = theta0*pi/180; \
	                                 rwidth = width*pi/180;')
	dihres.addGlobalParameter('fc_dih', fc_dih)
	dihres.addGlobalParameter('pi', 3.141592653589793)
	dihres.addPerTorsionParameter('width')
	dihres.addPerTorsionParameter('theta0')
	if os.path.isfile(restraint_file):
		for line in open(restraint_file, 'r'):
			segments = line.strip().split()
			atom1  = int(segments[0])
			atom2  = int(segments[1])
			atom3  = int(segments[2])
			atom4  = int(segments[3])
			theta0 = float(segments[4])
			width  = float(segments[5])
			dihres.addTorsion(atom1,atom2,atom3,atom4,[width,theta0])
	if os.path.isfile(restraint_file):
		for line in open(restraint_file, 'r'):
			segments = line.strip().split()
			atom1  = int(segments[0])
			atom2  = int(segments[1])
			atom3  = int(segments[2])
			atom4  = int(segments[3])
			theta0 = float(segments[4])
			width  = float(segments[5])
			dihres.addTorsion(atom1,atom2,atom3,atom4,[width,theta0])
	system.addForce(dihres)

	return system

##############################################


#############################################
print("\n> NPT Drude Equilibration! Let's do this:")
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCRD file = " + str(cmd.crd))
print("\tPSF file = " + str(cmd.psf))
print("\tToppar stream file = " + str(read_toppar(cmd.toppar)[1]))

print("\n\tSimulation_time = " + str(simulation_time))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps))

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tSave checkpoint each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure) + " (NPT ensemble)")
print("\tBB Restraint constant = " + str(fc_bb) + " (in kJ/mol/nm^2)")
print("\tSC Restraint constant = " + str(fc_sc) + " (in kJ/mol/nm^2)")
if dih_restraint_file != None:
	print("\tDihedral Restraint constant = " + str(fc_dih) + " (kJ/mol/rad^2)")

#############################################


print("\n> Setting the system:\n")
print("\t- Reading force field directory...")
charmm_params = read_toppar(cmd.toppar)[0]

print("\t- Reading topology and structure file...")
psf = CharmmPsfFile(cmd.psf)
crd = CharmmCrdFile(cmd.crd)

if cmd.state_pequil is None:
	print("\t- Setting box (using information on -state_c36 file)...")
	psf = get_cubic_box(psf, cmd.state_c36)
else:
	print("\t- Setting box (using information on -state_pequil file)...")
	psf = get_cubic_box(psf, cmd.state_pequil)

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

print("\t- Applying restraints... (using " + str(restraint_file) + ")")
if dih_restraint_file != None:
	print("\t- Applying dihedral restraints... (using " + str(dih_restraint_file) + ")")
	system = dih_restraints(system, crd, fc_dih, dih_restraint_file)
else:
	system = restraints(system, crd, fc_bb, fc_sc,restraint_file)

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

print("\t- Setting integrator and thermostat...")
integrator = DrudeLangevinIntegrator(temperature, 5/picosecond, 1*kelvin, 20/picosecond, 0.001*picoseconds)
integrator.setMaxDrudeDistance(0.02) # Drude Hardwall

simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()

equil_rst_file = cmd.state_pequil
if equil_rst_file is not None:
	print("\t- Reading previous equilibration state...")
	with open(equil_rst_file, 'r') as f:
		simulation.context.setState(XmlSerializer.deserialize(f.read()))
		simulation.context.setTime(0.0)
		simulation.currentStep = 0

print('\t- Using platform:', simulation.context.getPlatform().getName())

##########################################

dcd_file = jobname + ".dcd"
chk_file = jobname + ".chk"
log_file = jobname + ".log"
rst_file = jobname + ".rst"
prv_rst_file = jobname + ".rst"
crd_file = jobname + ".crd"


if os.path.exists(chk_file):
	simulation.loadCheckpoint(chk_file)

	print("> Restarting from checkpoint > " + chk_file + " <")

	# Calculating elapsed time
	chk_time = simulation.context.getState().getTime()
	chk_time_val = round(chk_time.value_in_unit(picosecond),4)
	chk_step = math.ceil(chk_time_val/dt.value_in_unit(picosecond))
	print("\t- Elapsed simulation time = " + str(round(chk_time_val,2)*picosecond) + " (step = " + str(chk_step) + ")")

	# Calculating remaining running time
	remaining_time = (simulation_time - chk_time).in_units_of(picosecond)
	remaining_nsteps = int(math.ceil(remaining_time.in_units_of(picosecond)/dt.in_units_of(picosecond)))

	# Adjust remaining running time
	simulation.currentStep = chk_step
	print("\t- Restarting from step " + str(chk_step) + " ...")

	dcd=DCDReporter(dcd_file, nsavcrd, append=True)
	print("\t- Appending to file " + jobname + ".dcd ...")

	backup_file = backup_old_log("*" + log_file + "*", log_file)
	print("\t- Backuping old log to " + backup_file + " ...")

	simulation.reporters.append(dcd)
	simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=nsteps, remainingTime=True, separator='\t\t'))
	simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
	simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))

	print("\n> Simulating " + str(remaining_nsteps) + " steps...")
	simulation.step(remaining_nsteps)

else:
	dcd = DCDReporter(dcd_file, nsavcrd)
	firstdcdstep = (nsteps) + nsavcrd
	dcd._dcd = DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0

	simulation.reporters.append(dcd)
	simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=nsteps, remainingTime=True, separator='\t\t'))
	simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
	simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))

	print("\n> Simulating " + str(nsteps) + " steps...")
	simulation.step(nsteps)

simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


##################################
# Writing last frame information
print("\n> Writing state file (" + str(rst_file) + ")...")
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open(rst_file, 'w') as f:
	f.write(XmlSerializer.serialize(state))

last_frame = int(nsteps/nsavcrd)
print("> Writing last coordinate (" + str(crd_file) + ", frame = " + str(last_frame) + ")...")
u = mda.Universe(cmd.psf, dcd_file)
system = u.select_atoms('all')
for ts in u.trajectory[(len(u.trajectory)-1):len(u.trajectory)]:
	with mda.Writer(str(crd_file), system.n_atoms, extended = True) as W:
		W.write(system)

try:
	quote = Quotes.getquote()
	print(quote)
except:
	pass

print("\n> Finished!\n")
