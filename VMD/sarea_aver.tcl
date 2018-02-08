## Calculate surface area for each residues and take an average number over the trajecotry.

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

## calculate surface area for each residues and take an average number over the trajecotry.
for {set j 1} {$j <= $nr} {incr j} {
    set sum 0
    for {set i 1} {$i < $nf} {incr i} {
         set allsel [atomselect top protein frame $i]
         set sel [atomselect top "protein and residue $j" frame $i]
         set rsasa [measure sasa 1.4 $allsel -restrict $sel]
         set sum [expr {$sum + $rsasa}]
    }
    set sum [expr {$sum/$nf}]
    puts $outfile "$sum"
}

exit     
