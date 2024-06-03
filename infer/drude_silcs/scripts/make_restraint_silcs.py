#!/usr/bin/python
###########################################
# Creates a restraint file to be used in positional restraints
# Atom naming is CHARMM compatible!!! Bear in mind your atom naming!
#
# USAGE: python make_restraint.py -pdb filename.pdb (-use_custom)
#
# By default, all Drude particles and lone pairs are ignored.
#
# Optionally, you can add add a customized restraining rule by adding
# patterns in the custom_dict object in line 30. MIND THE FORMAT!
#
# mdpoleto@vt.edu -> version 2023
###########################################
import os
import sys
import argparse
import MDAnalysis as mda
import warnings

warnings.filterwarnings("ignore", message="Element information is missing")
############################################################
bb_ref_list = ["CA", "C1'"]

sc_ref_dict = {
    "ADE": ["N1"],
    "GUA": ["N1"],
    "CYT": ["N3"],
    "THY": ["N3"],
    "URA": ["N3"]
}

res_ref_list = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER",
                "THR", "ASN", "GLN", "CYS", "GLY", "PRO", "ARG", "HIS", "LYS",
                "ASP", "GLU", "HSD", "HIE", "HSE", "HSP", "CYX", "HID", "GUA",
                "ADE", "URA", "THY", "CYT"]

ions_ref_list = ["POT", "CLA", "MG", "CAL", "SOD"]
exclude_list  = ["SWM4", "TIP3", "TIP", "WAT"]
############################################################

# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-pdb', type=str, default=None, required=True,
                help='Input coordinate file (.pdb)')

ap.add_argument('-use_custom', type=str, default=None, required=False,
                help='Whether or not add a custom restraining rule')

cmd = ap.parse_args()

file = cmd.pdb
rest = cmd.use_custom
#print(cmd.pdb)

if file != None:
	try:
		u = mda.Universe(file)
	except:
		raise ValueError(""">>> ERROR: Topology could not be loaded. Check your inputs!\n""")
else:
	raise ValueError(""">>> ERROR: Topology not found. Aren't you forgetting something?!\n""")

if rest != None:
	custom_selection = u.select_atoms(cmd.use_custom.replace('"',''))

################################################
def use_pdb(file):
	f = open(file,"r")

	bb_list = []
	sc_list = []
	custom_list = []

#	print("\n>>> Reading " + file + " file...")
	while True:
		line = f.readline()
		if not line: break

		if line.split()[0] == "ATOM" or line.split()[0] == "HETATM":
			atom_nr = int(line[6:11].strip(" "))
			atom_name = line[12:16].strip(" ")
			resname = line[17:21].strip(" ")
			resnumber = line[22:26].strip(" ")
			segid = line[72:76]

			atom_nr = str(atom_nr - 1) # 0-based numbering

			# ignoring Drude and lone pairs
			if atom_name[0] == "D" or atom_name[0:2] == "LP":
				continue
			else:

				# Ignoring water molecules based on those residue names
				if resname in exclude_list:
					continue

				else:
					#print(atom_nr, atom_name,resnumber, segid)
					if resname in res_ref_list:
						if atom_name in bb_ref_list:
							bb_list.append(atom_nr)
						if atom_name in sc_ref_dict[resname]:
							sc_list.append(atom_nr)
						if atom_name in ions_ref_list and atom_name[0] != "D" and atom_name[0:2] != "LP":
							bb_list.append(atom_nr)
					else:
						if resname not in exclude_list:
							#print(atom_nr, atom_name,resnumber, segid)
							#print(">>>>>> Residue " + str(resname+"-"+resnumber) + " not found in the RESIDUES/IONS list! Not restraining...")
							pass

	return bb_list, sc_list



################################################
bb_list, sc_list= use_pdb(file)

o = open("restraint_silcs.dat", "w")

for index in bb_list:
	string = " " + index + " SILCS\n"
	o.write(string)

for index in sc_list:
	string = " " + index + " SILCS\n"
	o.write(string)

if cmd.use_custom:
	for index in custom_selection.indices:
		string = " " + str(index) + " SILCS\n"
		o.write(string)

o.close()
#print(">>> Done!")
