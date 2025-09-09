#!/usr/bin/env python3

import sys
import re
import os

# Function to format the output in PDB format
def outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype):
	# Use PDB format
	#return f"{atom:<6}{anum:>5}{name:>5} {resn:<4}{chain}{resnr:>4}	{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}		{segid:<4}\n"
	# Use CHARMM CRD format
	return f"{atoms:10}{resnr:10}  {resn:<10}{name:<10}{x:18.10f}{y:20.10f}{z:20.10f}  {segid:<10}{resnr:<10}{0:18.10f}\n"

if len(sys.argv) == 5:
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	ions = sys.argv[3]
	nstrands = sys.argv[4]
else:
	print("Usage: python convert_gmx2drude_na.py input.pdb output.pdb ion_resid num_strands")
	sys.exit(1)

with open(input_file, "r") as infile:
	lines = infile.readlines()

# allow for multiple strands
if int(nstrands) > 26:
	print("Number of strands cannot exceed number of letters in the alphabet")
	sys.exit(1)
else:
	alpha = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	alphamap = dict(enumerate(alpha))

# starting counter for strand index
istrand = 0

crd_files = {}

# Open the output files
with open(output_file, "w") as fullout, open("system_properties.str", "w") as strout:
	strout.write("* system properties for SILCS\n*\n\n")

	# Counters for atoms and residue renumbering
	atoms = 0
	offset = 0

	# Save the Nucleic Acid sequence
	seq = ''

	# Counters for keeping track of the number of each molecule/chain length
	nres = nheta = nsolv = nmolatoms = nacey = nbenx = ndmee = nform = nimid = nmamy = nmeoh = nprpx = 0

	# Box information
	boxx = boxy = boxz = alp = bet = gam = 0
	shape = ''

	for i, line in enumerate(lines):
		# Get box size and shape information
		if line.startswith("CRYST1"):
			tmp = line.split()
			boxx = float(tmp[1])
			boxy = float(tmp[2])
			boxz = float(tmp[3])
			alp = float(tmp[4])
			bet = float(tmp[5])
			gam = float(tmp[6])

			if boxx == boxy:
				if boxy == boxz:
					shape = "cubic"
				else:
					shape = "tetragonal"
			else:
				if alp == 90.0 and gam == 90.0:
					if bet == 90.0:
						shape = "orthorhombic"
					else:
						shape = "monoclinic"
				else:
					shape = "triclinic"

		# Parse atoms
		if line.startswith("ATOM"):
			atoms += 1
			atom, anum, name, altlocation, resn, chain, resnr, codeinsert, x, y, z, occ, bfac, segid, atomtype= re.match(
				"(.{6})(.{5})(.{5})(.{1})(.{4})(.{1})(.{4})(.{4})(.{8})(.{8})(.{8})(.{6})(.{6})(.{4})(.{3})", line
			).groups()
			#print(atom, anum, name, altlocation, resn, chain, resnr, codeinsert, x, y, z, occ, bfac, segid, atomtype)
			anum = int(anum)
			resnr = int(resnr)
			x = float(x)
			y = float(y)
			z = float(z)
			bfac = float(bfac)
			occ = float(occ)

			# Remove whitespace from residue name
			resn = resn.strip()
			#print(str(name)+"yay")
			name = name.strip()
			segid = ''
		
			# If there are multiple strands, map back to an alpha character to designate
			# the last character, using base of DNA or RNA (below)
			# Since we start with istrand = 0 (A), then we increment the counter after finding
			# the 3'-terminus, as indicated by the H3T atom, which is unique to the 3TER patch
			char = alphamap[istrand]	# default to 'A'
			
			# for DNA
			if resn in ['DA', 'DC', 'DG', 'DT']:
				segid = "DNA" + char
				moltype = "DNA" + char
				resn = resn.replace('DA', 'ADE').replace('DC', 'CYT').replace('DG', 'GUA').replace('DT', 'THY')
				if "C5'" in name:
					nres += 1
					seq += " " + resn
				if "H3T" in name:
					istrand += 1
			# for RNA
			elif resn in ['ADE', 'CYT', 'GUA', 'URA']:
				segid = "RNA" + char
				moltype = "RNA" + char
				if "C5'" in name:
					nres += 1
					seq += " " + resn
				if "H3T" in name:
					istrand += 1

			elif resn == str(ions): 
				segid = "HETA"
				nheta += 1
				resnr = nheta

			elif resn == 'SOL':
				segid = "SOLV"
				resn = "SWM4"
				name = name.replace('OW', 'OH2').replace('HW1', 'H1').replace('HW2', 'H2')
				if "OH2" in name:
					nsolv += 1
					resnr = nsolv
				else: 
					resnr = nsolv

			elif resn == "BENX":
				segid = resn
				if "CG" in name:
					nbenx += 1
					resnr = nbenx
				else:
					resnr = nbenx

			elif resn == "PRPX":
				segid = resn
				if "H11" in name:
					nprpx += 1
					resnr = nprpx
				else:
					resnr = nprpx

			elif resn == "DMEE":
				segid = resn
				if "C1" in name:
					ndmee += 1
					resnr = ndmee
				else:
					resnr = ndmee

			elif resn == "FORM":
				segid = resn
				if "HA" in name:
					nform += 1
					resnr = nform
				else:
					resnr = nform

			elif resn == 'IMIA':
				resn = 'IMID'
				segid = resn
				if "CG" in name:
					nimid += 1
					resnr = nimid
				else:
					resnr = nimid

			elif resn == 'MEOH':
				segid = resn
				name = name.replace('CB', 'C1').replace('OG', 'O1').replace('HG1', 'HO1').replace('HB1', 'H1A').replace('HB2', 'H1B').replace('HB3', 'H1C')
				if "C1" in name:
					nmeoh += 1
					resnr = nmeoh
				else:
					resnr = nmeoh

			elif resn == "ACEY":
				segid = resn
				if "C1" in name:
					nacey += 1
					resnr = nacey
				else:
					resnr = nacey

			elif resn == "MAMY":
				segid = resn
				if "CE" in name:
					nmamy += 1
					resnr = nmamy
				else:
					resnr = nmamy

			else:
				segid = resn

			# Check for change in segid
			if atoms > 1:
				prevline = lines[i-1]
				if prevline.startswith("ATOM"):
					prevatom, prevanum, prevname, prevalt, prevresn, prevchain, prevresnr, previnsert, prevx, prevy, prevz, prevocc, prevb, prevsegid = re.match(
						"(.{6})(.{5})(.{5})(.{1})(.{4})(.{1})(.{4})(.{4})(.{8})(.{8})(.{8})(.{6})(.{6})(.{4})", prevline
					).groups()
					prevresnr = int(prevresnr)
					prevresn = prevresn.strip()
					prevsegid = ''
					# jal - please check

					if prevresn in ['DA', 'DC', 'DG', 'DT']:
						prevsegid = "DNA" + char
					elif prevresn in ['ADE', 'CYT', 'GUA', 'URA']:
						prevsegid = "RNA" + char										
					elif prevresn == 'POT':
						prevsegid = "HETA"
					elif prevresn in ['SOL', 'SWM4']:
						prevsegid = "SOLV"			
					else:
						if prevresn == 'IMIA':
							prevresn = "IMID"
							prevsegid = "IMID"
						prevsegid = prevresn

			# Element name for printing last field
			atomtype = name[0]
			fullout.write(outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype))

			# Count number of atoms and print associated segid lines to segid.crd	
			if segid.startswith("DNA") or segid.startswith("RNA"):
				nmolatoms += 1
				crd_filename = f"{segid[:3].lower()}{char.lower()}.crd"
				#print(crd_filename)
				
				if char not in crd_files:
					crd_files[char] = crd_filename
					with open(crd_filename, "w") as crd_file:
						crd_file.write(outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype))	
				else:
					# Append to the existing file
					with open(crd_files[char], "a") as crd_file:
						crd_file.write(outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype))	   
	
	# Print to stream file
	strout.write("set boxx	" + str(boxx) + "\n")
	strout.write("set boxy	" + str(boxy) + "\n")
	strout.write("set boxz	" + str(boxz) + "\n")
	strout.write("set alp	" + str(alp) + "\n")
	strout.write("set bet	" + str(bet) + "\n")
	strout.write("set gam	" + str(gam) + "\n")
	strout.write("set shape " + str(shape) + "\n\n")
	# account for multiple moltypes
	# compare first three characters of moltype variable
	strout.write("set nstrands = " + nstrands + "\n")

	for x in range(0,int(nstrands)):
		base = moltype[:3]	# account for DNA or RNA
		newmoltype =  base + alphamap[x]
		counter = x + 1		# fix zero-based numbering for readability
		strout.write("set moltype" + str(counter) + " " + str(newmoltype) + "\n")
	strout.write("set mol atoms " + str(nmolatoms) + "\n")
	strout.write("set nres " + str(nres) + "\n")
	strout.write("set ions " + str(ions) + "\n")
	strout.write("set nheta " + str(nheta) + "\n")
	strout.write("set nsolv " + str(nsolv) + "\n")
	strout.write("set nacey " + str(nacey) + "\n")
	strout.write("set nbenx " + str(nbenx) + "\n")
	strout.write("set ndmee " + str(ndmee) + "\n")
	strout.write("set nform " + str(nform) + "\n")
	strout.write("set nimid " + str(nimid) + "\n")
	strout.write("set nmamy " + str(nmamy) + "\n")
	strout.write("set nmeoh " + str(nmeoh) + "\n")
	strout.write("set nprpx " + str(nprpx) + "\n")

# reopen full.crd a bunch of times to add in header
with open(output_file,'r') as full_contents:
	  save = full_contents.read()

with open(output_file, 'w')as full_contents:
	header = "* CONVERTED FROM " + str(input_file) + " FOR DRUDE SILCS"

	full_contents.write
	full_contents.write(header + "\n")
	full_contents.write("*\n")
	full_contents.write(str(atoms).rjust(10) + " EXT\n")

with open(output_file,'a') as full_contents:
	full_contents.write(save)
	full_contents.write("TER\n")
	full_contents.write("ENDMDL\n")



for filename in crd_files.values():
	strandname =  filename[:-4]
	#print(strandname)
	# reopen to read content
	with open(filename,'r') as mol_contents:
		save = mol_contents.read()
		atomcount = save.count(strandname.upper())
		#print(atomcount)

	with open(filename, 'w') as mol_contents:
		
		header = "* CONVERTED FROM " + str(input_file) + " FOR DRUDE SILCS"

		mol_contents.write
		mol_contents.write(header + "\n")
		mol_contents.write("*\n")
		mol_contents.write(str(atomcount).rjust(10) + " EXT\n")

	with open(filename,'a') as mol_contents:
		mol_contents.write(save)
		mol_contents.write("TER\n")
		mol_contents.write("ENDMDL\n")

## TODO: update for multiple strands 
	os.rename(filename, filename.lower())	 
