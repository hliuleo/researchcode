## Calculate the radius distribution function and its integral between Li+ and Cl-

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


set sel1 [atomselect 0 "resname LI"]
set sel2 [atomselect 0 "resname CL"]
set gr1 [measure gofr $sel1 $sel2 delta .1 rmax 10 usepbc 1 selupdate 0 first 5000 last -1 step 1]
set outfile1 [open gofr.dat w]
set r [lindex $gr1 0]
set grd [lindex $gr1 1]
set igr [lindex $gr1 2]
foreach j $r k $grd l $igr {
   puts $outfile1 "$j $k $l"
}
close $outfile1
exit
