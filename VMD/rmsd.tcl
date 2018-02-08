## Calculate RMSD for each frames in a trajectory

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


set nf [molinfo 0 get numframes]
set outstk [open "rmsd.dat" w]
set native [atomselect top "resid 8 to 45 and name CA" frame 0]
for {set i 1} {$i < $nf} {incr i} {   
    set sel [atomselect top "resid 8 to 45 and name CA" frame $i]
    set rmsd  [measure rmsd $native $sel]
    puts $outstk "$rmsd"
}

exit
