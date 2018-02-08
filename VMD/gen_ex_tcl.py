#!/usr/bin/python
# This is used to generate a tcl script for VMD which is used to generate a dcd format trajectory file for each clusters.
# Usage: python gen_ex_tcl.py

openname = 'cluster_c2.5.dat'
f = open(openname,'r')
lines = f.readlines()
outname = 'extract_cluster.tcl'
fout = open(outname, 'w')
s = ''
s +="""set Syspdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]


#Read in the pdb and dcd file
mol new $Syspdb
animate read dcd $TrjDCD waitfor all 0

set numframe [molinfo 0 get numframes]

"""

count = 0
neighborlist = []
center = 0
cluster_size = 0
cluster_id = 1
for line in lines:
	line = line.split()
	outfile = 'clusters/cluster%d.dcd' % cluster_id
	center = int(line[1])
	cluster_size = int(line[3])
	neighborlist = [int(x) for x in line[4:]]
	frame = 'animate dup frame %d 0\n' % center
	s += frame
	for item in neighborlist:
		frame = 'animate dup frame %d 0\n' % item
        	s += frame
	s += 'set numadd %d\n' % (cluster_size+1)
	s += 'set last [expr $numframe + $numadd -1]\n'
	s += 'set sel [atomselect top all]\n'
	writeframe = 'animate write dcd %s beg $numframe end $last sel $sel 0\n' % outfile
	s += writeframe
	deloldframe = 'animate delete beg $numframe end $last\n'
	s += deloldframe + '\n'
	cluster_id += 1
s += 'exit\n'
fout.write(s)
