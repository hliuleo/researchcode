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
set outfile [open 2nd_struct.dat w]
set nr 45


## calculate the secondary structure for each frame and save them into a file
for {set i 1} {$i < $nf} {incr i} {
    animate goto $i
    vmd_calculate_structure 0 
    set sel [atomselect 0 "name CA"]
    set 2nd [$sel get structure]
    puts $outfile "$2nd"

}

exit     
