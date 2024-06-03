#!/usr/bin/env python3

import sys
import re

# Function to format the output in PDB format
def pdbfmt(atom, anum, name, resn, chain, resnr, x, y, z, occ, bfac, segid, atomtype):
    return f"{atom:<6}{anum:>5}{name:>5} {resn:<4}{chain}{resnr:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}      {segid:<4}\n"

def get_next_resname(lines, current_index):
    """
    Get the next residue name that is not 'ACE'.
    """
    for d in range(current_index + 1, len(lines)):
        line = lines[d]
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            resnr = line[22:26].strip()
            if resname != "ACE":
                return resname, resnr
    return None


if len(sys.argv) != 3:
    print("Usage: python {sys.argv[0]} input.pdb output.pdb")
    sys.exit(1)

prevresname = None
prevresnum = None
nextresn = None

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as infile:
    lines = infile.readlines()

# Open the output files
with open(output_file, "w") as outfile, open("system_properties.str", "w") as strout:
    strout.write("* system properties for SILCS\n*\n\n")

    # Counters for atoms and residue renumbering
    atoms = 0
    offset = 0

    # Save the PROTEIN sequence
    seq = ''

    # Counters for keeping track of the number of each molecule/chain length
    nprot = nheta = nsolv = nacey = nbenx = ndmee = nform = nimid = nmamy = nmeoh = nprpx = 0

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
            
            if resn in ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER",
                "THR", "ASN", "GLN", "CYS", "GLY", "PRO", "ARG", "HIS", "LYS",
                "ASP", "GLU", "HSD", "HIE", "HSE", "HSP", "CYX", "HID", "GUA"]:
                segid = "PROA"

                if "CA" in name:
                    nprot += 1
                    seq += " " + resn

            elif resn == 'ACE':
                nextresn = None

                if nextresn is None:
                    #print(i)
                    nextresn = get_next_resname(lines, i)[0]
                    #print(get_next_resname(lines, i))
                    resn = nextresn
                    nextresnr = get_next_resname(lines, i)[1]
                    resnr = int(nextresnr)
                    print(nextresnr)
                segid = "PROA"
            
            elif resn == "NME":
                if prevresname is not None:
                    #print(resn)
                    resn = prevresname
                    resnr = prevresnum
                    print(prevresnum)
                    segid = "PROA"
                    #print(prevresname, resn)

                #resn = resn.replace('DA', 'ADE').replace('DC', 'CYT').replace('DG', 'GUA').replace('DT', 'THY')

            elif resn == 'IMIA':
                resn = 'IMID'
                segid = resn

            elif resn == 'MEOH':
                segid = resn
                name = name.replace('CB', 'C1').replace('OG', 'O1').replace('HG1', 'HO1').replace('HB1', 'H1A').replace('HB2', 'H1B').replace('HB3', 'H1C')

            elif resn == 'POT':
                segid = "HETA"
                nheta += 1

            elif resn == 'SOL':
                segid = "SOLV"
                resn = "SWM4"
                name = name.replace('OW', 'OH2').replace('HW1', 'H1').replace('HW2', 'H2')
                if "OH2" in name:
                    nsolv += 1

            else:
                segid = resn

            prevresname = resn
            prevresnum = resnr

            if atoms > 1:
                prevline = lines[i-1]
                if prevline.startswith("ATOM"):
                    prevatom, prevanum, prevname, prevalt, prevresn, prevchain, prevresnr, previnsert, prevx, prevy, prevz, prevocc, prevb, prevsegid = re.match(
                        "(.{6})(.{5})(.{5})(.{1})(.{4})(.{1})(.{4})(.{4})(.{8})(.{8})(.{8})(.{6})(.{6})(.{4})", prevline
                    ).groups()
                    prevresnr = int(prevresnr)
                    prevresn = prevresn.strip()
                    prevsegid = ''

                    if prevresn in ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER",
                "THR", "ASN", "GLN", "CYS", "GLY", "PRO", "ARG", "HIS", "LYS",
                "ASP", "GLU", "HSD", "HIE", "HSE", "HSP", "CYX", "HID", "GUA", "ACE","NME"]:
                        prevsegid = "PROA"
                    elif prevresn == 'POT':
                        prevsegid = "HETA"
                    elif prevresn in ['SOL', 'SWM4']:
                        prevsegid = "SOLV"
                    else:
                        if prevresn == 'IMIA':
                            prevresn = "IMID"
                            prevsegid = "IMID"
                        prevsegid = prevresn

                    if prevsegid != segid:
                        offset = resnr - 1

                    if resn == "BENX":
                        if prevresnr != resnr:
                            nbenx += 1

                    elif resn == "PRPX":
                        if prevresnr != resnr:
                            nprpx += 1

                    elif resn == "DMEE":
                        if prevresnr != resnr:
                            ndmee += 1

                    elif resn == "MEOH":
                        if prevresnr != resnr:
                            nmeoh += 1

                    elif resn == "FORM":
                        if prevresnr != resnr:
                            nform += 1

                    elif resn == "IMID":
                        if prevresnr != resnr:
                            nimid += 1

                    elif resn == "ACEY":
                        if prevresnr != resnr:
                            nacey += 1

                    elif resn == "MAMY":
                        if prevresnr != resnr:
                            nmamy += 1
                        

            resnr -= offset
            

            # Element name for printing last field
            atomtype = name[0]

            outfile.write(pdbfmt(atom, anum, name, resn, chain, resnr, x, y, z, occ, bfac, segid, atomtype))

        else:
            outfile.write(line)

    outfile.write("TER\n")
    outfile.write("ENDMDL\n")

    # Print to stream file
    strout.write("set boxx  " + str(boxx) + "\n")
    strout.write("set boxy  " + str(boxy) + "\n")
    strout.write("set boxz  " + str(boxz) + "\n")
    strout.write("set alp   " + str(alp) + "\n")
    strout.write("set bet   " + str(bet) + "\n")
    strout.write("set gam   " + str(gam) + "\n")
    strout.write("set shape " + str(shape) + "\n\n")
    strout.write("set seq   " + str(seq) + "\n")
    strout.write("set nprot " + str(nprot) + "\n")
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
