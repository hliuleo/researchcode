#calculate RMSD for trajectory, align current snapshot with the one previous then calculate it. 
set nf [expr [molinfo top get numframes] -1]
set out [open rmsd_next.dat w]
set k 1
for {set j 2} {$j <= $nf} {incr j} {
    set i [expr {$j-$k}]
    set sel1 [atomselect top protein frame $i]
    set sel2 [atomselect top protein frame $j]
    set mfit [measure fit $sel2 $sel1]
    $sel2 move $mfit
    set rmsd [measure rmsd $sel2 $sel1]
    puts $out "$rmsd"
}
