import sys

in_file = sys.argv[1]
sys_name = sys.argv[2]
o=open(sys_name + "_boxsize.rst","w")

with open(in_file) as f:
    lines = f.readlines()
    for line in lines:
        if "boxx" in line:
            boxx=float(line.split()[2]) * 0.1
#            print(boxx)
        if "boxy" in line:
            boxy=float(line.split()[2]) * 0.1
#            print(boxy)
        if "boxz" in line:
            boxz=float(line.split()[2]) * 0.1
#            print(boxz)

string='<?xml version="1.0" ?>\n\
<State>\n\
\t<PeriodicBoxVectors>\n\
\t\t<A x="{}" y="0" z="0"/>\n\
\t\t<B x="0" y="{}" z="0"/>\n\
\t\t<C x="0" y="0" z="{}"/>\n\
\t<PeriodicBoxVectors>\
'.format(boxx,boxy,boxz)

o.write(string)