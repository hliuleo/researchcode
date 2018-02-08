## Calculate surface area for each residues of each frames.

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

set first [atomselect top protein frame 0]
set nf [molinfo top get numframes]
set outfile [open sarea.dat w]
set nr 45

## Note in tcl script, the index of residue starts from 0.
for {set i 1} {$i < $nf} {incr i} {
    for {set j 0} {$j <= $nr-1} {incr j} {
         set allsel [atomselect top protein frame $i]
         set sel [atomselect top "protein and residue $j" frame $i]
         set rsasa [measure sasa 1.4 $allsel -restrict $sel]
	 puts $outfile [format "%.2f" $rsasa]
    }
}

exit     
