# save the dcd trajectory file into the crd format. This is used for calculating native contacts by Haijun's code.

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

set nf [molinfo 0 get numframes]
set sel [atomselect 0 protein]
set last [expr $nf -1]
animate write crd "all.crd" beg 1 end $last sel $sel waitfor all 0
exit
