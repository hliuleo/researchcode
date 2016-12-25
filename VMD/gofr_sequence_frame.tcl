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
set outfile1 [open Li_all_byFrame.dat w]
set outfile2 [open Cl_all_byFrame.dat w]

## print number of frames in the current trajectory

set num_frames [expr [molinfo 0 get numframes] - 1]

proc calIonDistByFrame {idx outfile1 outfile2} {
    ## calculate gofr for each residues. Only one or two atoms are selected to represent the residues (from Qiang Cui,JACS,2006)
    set sel1 [atomselect 0 "protein and resid 1 and (name NE2 or name OE1)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 2 and (name ND2 or name OD1)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 3 and (name ND2 or name OD1)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 4 and (name OD1 or name OD2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 5 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 6 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 7 and name OG"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 8 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 9 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 10 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 11 and (name NH1 or name NH2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 12 and (name NH1 or name NH2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 13 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 14 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 15 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 16 and (name OE2 or name OE1)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 17 and (name ND1 or name NE2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 18 and (name ND2 or name OD1)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 19 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 20 and (name OD1 or name OD2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 21 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 22 and name OG"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 23 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 24 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 25 and name NZ"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 26 and name CA"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 27 and name OG1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 28 and name CA"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 29 and name CG1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 30 and name CA"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 31 and name CA"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 32 and (name NH1 or name NH2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 33 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 34 and name OG1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 35 and (name NH1 or name NH2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 36 and (name OE1 or name OE2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 37 and (name OD1 or name OD2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 38 and name CG1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 39 and (name OE1 or name OE2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 40 and name NZ"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 41 and (name ND1 or name NE2)"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 42 and name CD1"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 43 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 44 and name NZ"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }
    set sel1 [atomselect 0 "protein and resid 45 and name CB"]
    set sel2 [atomselect 0 "resname LI"]
    set sel3 [atomselect 0 "resname CL"]
    set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set gr2 [measure gofr $sel1 $sel3 delta .1 rmax 10 usepbc 1 selupdate 0 first $idx last $idx step 1]
    set r [lindex $gr1 0]
    set grd [lindex $gr1 1]
    set igr [lindex $gr1 2]
    foreach j $r k $grd l $igr {
        puts $outfile1 "$j $k $l"
    }
    set r [lindex $gr2 0]
    set grd [lindex $gr2 1]
    set igr [lindex $gr2 2]
    foreach j $r k $grd l $igr {
        puts $outfile2 "$j $k $l"
    }

}

for {set idx 1} {$idx <= $num_frames} {incr idx} {
    calIonDistByFrame $idx $outfile1 $outfile2
}
close $outfile1
close $outfile2
exit
