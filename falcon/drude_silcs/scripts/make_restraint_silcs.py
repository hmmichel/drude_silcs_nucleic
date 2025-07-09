#!/usr/bin/python
###########################################
# Creates a restraint file to be used in positional restraints
# Atom naming is CHARMM compatible!!! Bear in mind your atom naming!
#
# USAGE: python make_restraint.py -crd filename.crd (-use_custom)
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

ap.add_argument('-crd', type=str, default=None, required=False,
                help='Input coordinate file (.crd)')

ap.add_argument('-pdb', type=str, default=None, required=False,
                help='Input coordinate file (.pdb)')

ap.add_argument('-use_custom', type=str, default=None, required=False,
                help='Whether or not add a custom restraining rule')

cmd = ap.parse_args()

file = cmd.pdb
filecrd = cmd.crd
rest = cmd.use_custom
#print(cmd.pdb)

if filecrd != None:
	try:
		u = mda.Universe(filecrd)
	except:
		raise ValueError(""">>> ERROR: Topology could not be loaded. Check your inputs!\n""")
else:
	raise ValueError(""">>> ERROR: Topology not found. Aren't you forgetting something?!\n""")

if rest != None:
	custom_selection = u.select_atoms(cmd.use_custom.replace('"',''))

################################################
def use_crd(file):
        f = open(file,"r")

        line = ""

        bb_list = []
        sc_list = []
        
        flag=False
        print("\n>>> Reading " + file + " file...")
        while True:
                line = f.readline()
                if not line: break

                #print(len(line.split()), line)

                if len(line.split()) > 1:
                        if "*" in line:
                                continue
                        elif "EXT" in line.split():
                                flag = True
                                continue

                        if flag == True:
                                atom_nr = int(line.split()[0])
                                resname = line.split()[2]
                                atom_name = line.split()[3]
                                segid = line.split()[7]
                                resnumber = line.split()[8]

                                #print(atom_nr, resnumber, resname, atom_name, segid)

                                atom_nr = str(atom_nr - 1)

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
                                                        if resname in sc_ref_dict: 
                                                                if atom_name in sc_ref_dict[resname]:
                                                                        sc_list.append(atom_nr)
                                                        
                                                       #if atom_name in ions_ref_list and atom_name[0] != "D" and atom_name[0:2] != "LP":
                                                       #         bb_list.append(atom_nr)
                                                else:
                                                        if resname not in exclude_list:
                                                                #print(atom_nr, atom_name,resnumber, segid)
                                                                #print(">>>>>> Residue " + str(resname+"-"+resnumber) + " not found in the RESIDUES/IONS list! Not restraining...")
                                                                pass

        return bb_list, sc_list

def use_pdb(file):
	f = open(file,"r")

	bb_list = []
	sc_list = []

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
						if resname in sc_ref_dict:
				  			if atom_name in sc_ref_dict[resname]:
							        sc_list.append(atom_nr)
						#if atom_name in ions_ref_list and atom_name[0] != "D" and atom_name[0:2] != "LP":
						#	bb_list.append(atom_nr)
					else:
						if resname not in exclude_list:
							#print(atom_nr, atom_name,resnumber, segid)
							#print(">>>>>> Residue " + str(resname+"-"+resnumber) + " not found in the RESIDUES/IONS list! Not restraining...")
							pass

	return bb_list, sc_list



################################################
if filecrd != None:
        bb_list, sc_list = use_crd(filecrd)
else:
        bb_list, sc_list = use_pdb(file)

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
