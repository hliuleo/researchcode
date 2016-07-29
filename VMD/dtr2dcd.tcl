##  apply pbc and align trajectories to the native structure then write trj in the dcd format

## import pbctool plugin
package require pbctools

set SysPdb [lindex $argv 0]
set TrjSTK [lindex $argv 1]
set Prefix [lindex $argv 2]

mol new $SysPdb
set instk [open $TrjSTK r]
while { [gets $instk dtr] >= 0} {
        animate read dtr $dtr waitfor all 0
}

set nf [molinfo 0 get numframes]
set last [expr $nf -1]

## wrap water and salts around protein
pbc wrap -compound residue -center com -centersel "protein" -all



for {set i 1} {$i < $nf} {incr i} {
     set sel1 [atomselect top "(resid 8 to 45) and name CA" frame 0]
     set sel2 [atomselect top "(resid 8 to 45) and name CA" frame $i]
     set mfit [measure fit $sel2 $sel1]
     set sel2all [atomselect top all frame $i]
     $sel2all move $mfit
     $sel2 delete
     $sel1 delete
     $sel2all delete
# delete the data stored in an array and release memory
     unset sel2
     unset sel1
     unset mfit
     unset sel2all
}

set outdcd [atomselect top all]
animate write dcd ${Prefix}.dcd beg 1 end $last sel $outdcd 0

quit
