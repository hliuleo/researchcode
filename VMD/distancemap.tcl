set SysPdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]
set NumClu [lindex $argv 2]

mol new $SysPdb
for {set i 0} {$i < $NumClu} {incr i} {
     set file_name "cluster"
     append file_name [expr $i+1]
     append file_name ".pdb"
     puts $file_name
     mol addfile $file_name
}

set nf [molinfo 0 get numframes]
set outfile [open "center_info.dat" w]
set nr 45

## calculate surface area for each residues and store in a list
set sarealist [list]
proc sarea_cal {nf nr} {
     global sarealist
     for {set i 1} {$i < $nf} {incr i} {
         set sum 0
         set sarea {}
         for {set j 1} {$j <= $nr} {incr j} {
              set allsel [atomselect top protein frame $i]
              set sel [atomselect top "protein and residue $j" frame $i]
              set rsasa [measure sasa 1.4 $allsel -restrict $sel]
	      lappend sarea $rsasa
         }
         lappend sarealist $sarea
         unset sarea
     }
}

set rmsdlist [list]
proc rmsd_cal {nf} {
     global rmsdlist
     ## Calculate Ca (8-45) RMSD. Reference structure is the native state.
     for {set i 1} {$i < $nf} {incr i} {
         set ref [atomselect top "protein and resid 8 to 45 and name CA" frame 0]
         set comp [atomselect top "protein and resid 8 to 45 and name CA" frame $i]
         set rmsd [measure rmsd $ref $comp]
         lappend rmsdlist $rmsd
     }
}

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
          set sel1 [atomselect 0 "protein and resid 2 and name ND"]
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
          set sel1 [atomselect 0 "protein and resid 4 and name OD1"]
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
          set sel1 [atomselect 0 "protein and resid 6 and name CB"]
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
          set sel1 [atomselect 0 "protein and resid 10 and name CB"]
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
          set sel1 [atomselect 0 "protein and resid 13 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 14 and name CB"]
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
          set sel1 [atomselect 0 "protein and resid 16 and name NE2"]
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
          set sel1 [atomselect 0 "protein and resid 19 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 20 and name OD1"]
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
          set sel1 [atomselect 0 "protein and resid 24 and name CB"]
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
          set sel1 [atomselect 0 "protein and resid 29 and name CB"]
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
          set sel1 [atomselect 0 "protein and resid 33 and name CB"]
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
          set sel1 [atomselect 0 "protein and resid 36 and name OE1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 37 and name OD1"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 38 and name CB"]
          set sel2 [atomselect 0 "resname LI"]
          set sel3 [atomselect 0 "resname CL"]
          set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $i last $i step 1]
          set igrli [lindex [lindex $gr1 2] 61]
          set igrcl [lindex [lindex $gr2 2] 61]
          lappend Lidist $igrli
          lappend Cldist $igrcl
          set sel1 [atomselect 0 "protein and resid 39 and name OE1"]
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
          set sel1 [atomselect 0 "protein and resid 42 and name CB"]
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

rmsd_cal $nf
sarea_cal $nf $nr
saltdis_cal $nf

set count 1
foreach i $rmsdlist j $sarealist k $Cldist_list o $Lidist_list {
     puts $outfile "Cluster$count"
     puts $outfile "RMSD"
     puts $outfile $i
     puts $outfile "Surface area along sequence"
     puts $outfile $j
     puts $outfile "Cl distribution along sequence"
     puts $outfile $k
     puts $outfile "Li distribution along sequence"
     puts $outfile $o
     puts $outfile "\n"
     incr count
}

exit
