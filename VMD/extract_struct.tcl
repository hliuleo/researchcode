#Usage vmd -dispdev text -e extract_struct.tcl -args {structure_file} {traj_file} {nameofnewtrajecory}


## check if there are correct number of arguments
if { [llength $argv] != 3} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysPdb [lindex $argv 0];
set TrjDCD [lindex $argv 1];
set Prefix [lindex $argv 2];

## first creat a new molecule at $SysMae and then append trajecotries to the molecule.
mol new $SysPdb
animate read dcd $TrjDCD waitfor all 0


if { [file exists ${Prefix}.dcd] } {
		puts "ERROR: ${Prefix}.dcd exists, delete it or specify a diffrent prefix"
		exit
	}


set nf [molinfo 0 get numframes ]
pbc wrap -compound residue -center com -centersel "protein" -all

## Align all the structures to frame 0 (native state) and calculate RMSD. If the current snapshot reaches certain criteria, then save it.
set numadd 0

for {set i 0} {$i < $nf} {incr i} {
     set sel1 [atomselect top "(resid 8 to 45) and name CA" frame 0]
     set sel2 [atomselect top "(resid 8 to 45) and name CA" frame $i]
     set mfit [measure fit $sel2 $sel1]
     set sel2all [atomselect top all frame $i]
     $sel2all move $mfit
     set rmsd [measure rmsd $sel2 $sel1]
     if {$rmsd < 5 && $i != 0} {
        animate dup frame $i 0
# save the current snapshot into the end of the trajectory.

	incr numadd
     }
     $sel2 delete
     $sel1 delete
     $sel2all delete
#               mfit delete
     unset sel2
     unset sel1
     unset rmsd
     unset mfit
     unset sel2all
}

set last [expr $nf +$numadd -1]
set sel [atomselect top all]

#write out all the saved snapshots into a file.
animate write dcd ${Prefix}.dcd beg $nf end $last sel $sel 0

exit
