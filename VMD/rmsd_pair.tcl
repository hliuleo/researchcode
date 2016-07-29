# function to calculate pair-wise RMSD, output N-1 + N-2 + N-3 + .... + 1 pair-wise RMSD.
proc rmsd {nf} {
	set k 1
	set last [expr {$nf -1}]
	set out [open rmsd_pair.dat w]
        for {set j 1} {$j < $last} {incr j} {
	    set ref [atomselect 0 all frame $j]
            for {set i [expr {$j+$k}]} {$i < $nf} {incr i} {
                 set comp [atomselect 0 all frame $i]
		 set all [atomselect 0 all frame $i]
		 set mfit [measure fit $comp $ref]
                 $all move $mfit
                 set rmsd [measure rmsd $ref $comp]
                 puts $out [format "%.2f" $rmsd]
                 $all delete
                 $comp delete
#                 $mfit delete
                 unset comp
		 unset all
                 unset rmsd
		 unset mfit
            }
	    $ref delete
	    unset ref
        }
}


#Check number of argument
if { [llength $argv] != 2} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysPdb [lindex $argv 0];
set TrjDCD [lindex $argv 1];

#Read in the pdb and dcd file
mol new $SysPdb
animate read dcd $TrjDCD waitfor all 0

set nf [molinfo 0 get numframes]
rmsd $nf

exit
