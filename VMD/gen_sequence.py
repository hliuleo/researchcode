# Use to generate a part of tcl script 'salt_dist.tcl'

f = open('/home/hliu/Project/BBL_Salt/1w4h_li/doc/representative_atom_potential.dat','r')
s = ''
blank = '          '
lines = f.readlines()
for line in lines:
	line = line.strip()
	select_atom = blank + 'set sel1 [atomselect 0 "%s"]\n' % line
	s += select_atom
	gorf = blank + """set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl\n"""
	s += gorf
print s	
