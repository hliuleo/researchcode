#For given snapshots, calculate pair-wise RMSD between each other. 
#This script is only compatible with the clustering file generated from MSMBuilder and June's code "cluster_info.py"

set Syspdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]
set Prefix [lindex $argv 2]

#Read in the pdb and dcd file (the whole trajectory, not the one of cluster)
mol new $Syspdb
animate read dcd $TrjDCD waitfor all 0

#Open a file which saves the index of conformations (like file microstate_*)
set fp [open $Prefix r]
set file_data [read $fp]
set data [split $file_data "\n"]

# Define an empty list
set index [list]                                   
foreach i $data {

#Get the index of conformations and store in a list called index
	set item [split $i " "]
	set ll [expr [lindex $item 1] + 1]
#	puts $ll
	lappend index $ll
}

set file_name "rmsd_pair"
append file_name "_" $Prefix
append file_name "" ".dat"

#Delete last element in the list, because when read in the file the end of file will be treated as " " and stored in the list.
set last [expr [llength $index] -1]
set index [lreplace $index $last $last]
puts [llength $index]
set out [open $file_name w]

#Calcualte pair-wise RMSD
foreach i $index {
	set ref [atomselect top all frame $i]
	foreach j $index {
		set comp [atomselect top all frame $j]
		set mfit [measure fit $comp $ref]
		$comp move $mfit
		set rmsd [measure rmsd $ref $comp]
		puts $out [format "%.2f" $rmsd]
		$comp delete
		unset comp
		unset mfit
		unset rmsd
	}
	$ref delete
	unset ref
}

exit
