# script to offset the Drude SILCS GFE values to zero based on user input
#
#  usage: python offset_gfe.py system_name offser_value
# example: python offset_gfe.py 2lwk 1.3
#
# Happy simulating! HMM 10/19/2024

import sys
import os
import subprocess

# assign system name and offset value from command line
sysname = sys.argv[1]
offset = sys.argv[2]

# generate list of SILCS solutes
frags = ["swm4o", "acec", "mamn", "aceo", "meoo", "hbdon", "imin", "forn", "prpc", "apolar", "gehc", "benc", "hbacc", "forc", "foro", "iminh", "dmeo"]

# rename original map files with no offsert to "old" so we don't overwrite them
subprocess.run("rename s/map/map_old/ " +  sysname + "*gfe.map", shell=True)

# iterate through each solute and generate a new map file containing the offset gfe values
for frag in frags:
    
    dx_file=sysname + "." + frag + ".gfe.dx"
    #print(dx_file)
    
    # if dx files are present in maps/ pymol will use them for visualization
    # must delete any old ones in order to view the new gfe maps
    if os.path.isfile(dx_file):
        os.remove(dx_file)

    # specify names of old map and new map files
    old_map = open(sysname + "." + frag + ".gfe.map_old", "r")
    
    new_map = open(sysname + "." + frag +".gfe.map", "w")
    
    # for the header lines copy them over to the new file
    # ensures correct box size and location of occupancies
    for line in old_map.readlines():
        if line[0] != " ":
            new_map.write(line)
        else:
            old_gfe = line[4:10]
            new_gfe = float(old_gfe) + float(offset)
            new_gfe = ("%.3f" % new_gfe)
            #print(new_gfe)


 
            newline = "    " + str(new_gfe) + "\n"
            #print (newline)

            new_map.write(newline)
            continue
    
    new_map.close()

# the fragmap pymol plugin is hard coded to read "tipo" for water
# create symlink to read "swm4o" file as the "tipo" file
if os.path.isfile(sysname + ".tipo.gfe.map"):
    pass
else:
    os.symlink(sysname + ".swm4o.gfe.map", sysname + ".tipo.gfe.map")




