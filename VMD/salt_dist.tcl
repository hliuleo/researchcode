# Calculate the salt distribution along sequence and save the value at 6A. The atoms used to represent residues are consistent with the ones used in assignment of electrostatic potentials.

set SysPdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]

mol new $SysPdb
animate read dcd $TrjDCD waitfor all 0
set nf [molinfo 0 get numframes]
set clout [open "Cl_dist.dat" w]
set liout [open "Li_dist.dat" w]
set nr 45

## Calculate salt distribution along sequence

set Cldist_list [list]
set Lidist_list [list]
proc saltdis_cal {nf} {
     global Lidist_list
     global Cldist_list
     for  {set i 1} {$i < $nf} {incr i} {
	  set Lidist {}
          set Cldist {}
          set sel1 [atomselect 0 "protein and resid 1 and name NE2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 2 and name ND2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 3 and name ND2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 4 and name OD2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 5 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 6 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 7 and name OG"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 8 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 9 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 10 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 11 and name NH1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 12 and name NH1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 13 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 14 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 15 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 16 and name OE2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 17 and name ND1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 18 and name ND2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 19 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 20 and name OD2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 21 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 22 and name OG"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 23 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 24 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 25 and name NZ"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 26 and name CA"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 27 and name OG1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 28 and name CA"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 29 and name CG1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 30 and name CA"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 31 and name CA"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 32 and name NH1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 33 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 34 and name OG1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 35 and name NH1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 36 and name OE2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 37 and name OD2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 38 and name CG1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 39 and name OE2"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 40 and name NZ"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 41 and name ND1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 42 and name CD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 43 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 44 and name NZ"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 45 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
	  lappend Lidist_list $Lidist
          lappend Cldist_list $Cldist
          unset Lidist
          unset Cldist
     }
}

saltdis_cal $nf

set count 1
foreach i $Cldist_list j $Lidist_list {
	foreach l $i h $j {
		puts $clout $l
		puts $liout $h
	}
}

exit
          

