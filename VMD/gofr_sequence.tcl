# Use to calculate gofr and integral of gofr between Li, Cl and selected protein atoms. protein atoms showed here are consistent with representative_atom.dat
## check if there are correct number of arguments
if { [llength $argv] != 2} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysPdb [lindex $argv 0];
set TrjDCD [lindex $argv 1];

## first creat a new molecule at $SysMae and then append trajecotries to the molecule.
mol new $SysPdb
animate read dcd $TrjDCD waitfor all 0

## print number of frames in the current trajectory
puts [molinfo 0 get numframes]

## calculate gofr for each residues. Only one or two atoms are selected to represent the residues (from Qiang Cui,JACS,2006)
set sel1 [atomselect 0 "protein and resid 15 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res15_Na_all.dat w]
set outfile2 [open res15_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 21 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res21_Na_all.dat w]
set outfile2 [open res21_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 23 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res23_Na_all.dat w]
set outfile2 [open res23_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 43 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res43_Na_all.dat w]
set outfile2 [open res43_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 45 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res45_Na_all.dat w]
set outfile2 [open res45_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 5 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res5_Na_all.dat w]
set outfile2 [open res5_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 9 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res9_Na_all.dat w]
set outfile2 [open res9_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 11 and (name NH1 or name NH2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res11_Na_all.dat w]
set outfile2 [open res11_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 12 and (name NH1 or name NH2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res12_Na_all.dat w]
set outfile2 [open res12_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 32 and (name NH1 or name NH2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res32_Na_all.dat w]
set outfile2 [open res32_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 35 and (name NH1 or name NH2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res35_Na_all.dat w]
set outfile2 [open res35_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 18 and (name ND2 or name OD1)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res18_Na_all.dat w]
set outfile2 [open res18_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 2 and (name ND2 or name OD1)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res2_Na_all.dat w]
set outfile2 [open res2_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 3 and (name ND2 or name OD1)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res3_Na_all.dat w]
set outfile2 [open res3_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 20 and (name OD1 or name OD2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res20_Na_all.dat w]
set outfile2 [open res20_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 37 and (name OD1 or name OD2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res37_Na_all.dat w]
set outfile2 [open res37_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 4 and (name OD1 or name OD2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res4_Na_all.dat w]
set outfile2 [open res4_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 1 and (name NE2 or name OE1)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res1_Na_all.dat w]
set outfile2 [open res1_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 16 and (name OE2 or name OE1)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res16_Na_all.dat w]
set outfile2 [open res16_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 36 and (name OE1 or name OE2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res36_Na_all.dat w]
set outfile2 [open res36_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 39 and (name OE1 or name OE2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res39_Na_all.dat w]
set outfile2 [open res39_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 26 and name CA"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res26_Na_all.dat w]
set outfile2 [open res26_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 28 and name CA"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res28_Na_all.dat w]
set outfile2 [open res28_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 30 and name CA"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res30_Na_all.dat w]
set outfile2 [open res30_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 31 and name CA"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res31_Na_all.dat w]
set outfile2 [open res31_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 17 and (name ND1 or name NE2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res17_Na_all.dat w]
set outfile2 [open res17_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 41 and (name ND1 or name NE2)"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res41_Na_all.dat w]
set outfile2 [open res41_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 10 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res10_Na_all.dat w]
set outfile2 [open res10_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 24 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res24_Na_all.dat w]
set outfile2 [open res24_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 13 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res13_Na_all.dat w]
set outfile2 [open res13_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 14 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res14_Na_all.dat w]
set outfile2 [open res14_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 19 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res19_Na_all.dat w]
set outfile2 [open res19_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 33 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res33_Na_all.dat w]
set outfile2 [open res33_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 42 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res42_Na_all.dat w]
set outfile2 [open res42_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 6 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res6_Na_all.dat w]
set outfile2 [open res6_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 25 and name NZ"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res25_Na_all.dat w]
set outfile2 [open res25_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
save the dcd trajectory file into the crd format. This is used for calculating native contacts by Haijun's code.t i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 40 and name NZ"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res40_Na_all.dat w]
set outfile2 [open res40_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 44 and name NZ"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res44_Na_all.dat w]
set outfile2 [open res44_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 8 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res8_Na_all.dat w]
set outfile2 [open res8_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 22 and name OG"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res22_Na_all.dat w]
set outfile2 [open res22_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 7 and name OG"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res7_Na_all.dat w]
set outfile2 [open res7_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 27 and name OG1"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res27_Na_all.dat w]
set outfile2 [open res27_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 34 and name OG1"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res34_Na_all.dat w]
set outfile2 [open res34_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 29 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res29_Na_all.dat w]
set outfile2 [open res29_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2
set sel1 [atomselect 0 "protein and resid 38 and name CB"]
set sel2 [atomselect 0 "resname LI"]
set sel3 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open res38_Na_all.dat w]
set outfile2 [open res38_Cl_all.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
set r [lindex $gr2 0]
set grd [lindex $gr2 1]
set igr [lindex $gr2 2]
set i 0
foreach j $r k $grd l $igr {
   puts $outfile2 "$j $k $l"
}
close $outfile2

exit
