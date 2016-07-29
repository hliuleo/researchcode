# Use to extract cluster snapshot and save to pdb

## check if there are correct number of arguments
if { [llength $argv] != 2} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysPdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]

mol new $SysPdb
animate read dcd $TrjDCD waitfor all 0


# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]
for {set i 1} {$i < $nf} {incr i} {
    set out [atomselect 0 "protein or resname NME" frame $i]
    set outname $i
    append outname ".pdb"
    $out writepdb $outname
}
exit
