#!/usr/bin/env python3

import sys
import re
import os

# set number of strands as a variable based on user input
strands =  sys.argv[1:]

# Function to format the output in PDB format
def outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype):
    # Use PDB format
    #return f"{atom:<6}{anum:>5}{name:>5} {resn:<4}{chain}{resnr:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}      {segid:<4}\n"
    # Use CHARMM CRD format
    return f"{atoms:10}{resnr:10}  {resn:<10}{name:<10}{x:18.10f}{y:20.10f}{z:20.10f}  {segid:<10}{resnr:<10}{0:18.10f}\n"

if len(sys.argv) == 4:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ions = sys.argv[3] 
elif len(sys.argv) == 3: 
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ions = ""
else:
    print("Usage: python convert_gmx2drude_na.py input.pdb output.pdb ion_resid")
    sys.exit(1)

with open(input_file, "r") as infile:
    lines = infile.readlines()

# Open the output files
with open(output_file, "w") as fullout, open("system_properties.str", "w") as strout, open("seg.crd", "w") as segcrd:
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

            # for DNA
            if resn in ['DA', 'DC', 'DG', 'DT']:
                segid = "DNAA"
                moltype = "DNAA"
                resn = resn.replace('DA', 'ADE').replace('DC', 'CYT').replace('DG', 'GUA').replace('DT', 'THY')
                if "C5'" in name:
                    nres += 1
                    seq += " " + resn
            # for RNA
            elif resn in ['ADE', 'CYT', 'GUA', 'URA']:
                segid = "RNAA"
                moltype = "RNAA"
                if "C5'" in name:
                    nres += 1
                    seq += " " + resn
            
            elif resn == str(ions): 
                segid = "HETA"
                nheta += 1

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

                    if prevresn in ['DA', 'DC', 'DG', 'DT']:
                        prevsegid = "DNAA"
                    elif prevresn in ['ADE', 'CYT', 'GUA', 'URA']:
                        prevsegid = "RNAA"
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

            # Count number of protein atoms and print "PROA" lines to proa.crd
            if segid == "DNAA":
                nmolatoms += 1
                segcrd.write(outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype))       
            elif segid == "RNAA":
                nmolatoms += 1
                segcrd.write(outfmt(atom, anum, name, resn, resnr, x, y, z, occ, bfac, segid, atomtype))
    
    # Write terminal to out files
    fullout.write("TER\n")
    fullout.write("ENDMDL\n")
    
    segcrd.write("TER\n")
    segcrd.write("ENDMDL\n")

    # Print to stream file
    strout.write("set boxx  " + str(boxx) + "\n")
    strout.write("set boxy  " + str(boxy) + "\n")
    strout.write("set boxz  " + str(boxz) + "\n")
    strout.write("set alp   " + str(alp) + "\n")
    strout.write("set bet   " + str(bet) + "\n")
    strout.write("set gam   " + str(gam) + "\n")
    strout.write("set shape " + str(shape) + "\n\n")
    strout.write("set moltype " + str(moltype) + "\n")
    strout.write("set mol atoms " + str(nmolatoms) + "\n")
    strout.write("set seq   " + str(seq) + "\n")
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

moltype = str(moltype)

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


# reopen proa.crd a bunch of times to add in header
with open("seg.crd",'r') as mol_contents:
      save = mol_contents.read()

with open("seg.crd", 'w')as mol_contents:
    header = "* CONVERTED FROM " + str(input_file) + " FOR DRUDE SILCS"

    mol_contents.write
    mol_contents.write(header + "\n")
    mol_contents.write("*\n")
    mol_contents.write(str(nmolatoms).rjust(10) + " EXT\n")

with open("seg.crd",'a') as mol_contents:
    mol_contents.write(save)

os.rename("seg.crd", moltype.lower() + ".crd")   
