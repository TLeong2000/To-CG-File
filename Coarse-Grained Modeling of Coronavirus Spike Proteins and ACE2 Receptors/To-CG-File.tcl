# Timothy Leong
# The pdb file must be in your current directory

package require cggui
package require cgtools
if {[lsearch -exact $auto_path [pwd]] == -1} {
	lappend $auto_path [pwd]
	package require cgtools_patch 1.0
	lreplace $auto_path [lsearch -exact $auto_path [pwd]] [lsearch -exact $auto_path [pwd]]
} else {
	package require cgtools_patch 1.0
}

# All-atom PDB I/O

set move_on 0

while {$move_on == 0} {
	puts "Enter the path of the pdb file to be converted (or exit to exit): \n"

	gets stdin filename

	if {$filename == "exit"} {
		puts "Terminating operation."
		return
	}
	
	set error_value [catch "mol new $filename" error_desc]
	if {$error_value == 1} {
		puts "The following error occurred while loading in the PDB file: \n"
		puts $error_desc
		puts "\nPlease try again."
	} else {
		set move_on 1
	}
	
}

set move_on 0

while {$move_on == 0} {

	puts "Do you have a Protein Structure File ready for this pdb file?"
	puts "Enter the path of the psf file if you do, or no if you do not."
	puts "If you do not have a psf file, autoPSF will generate one for you. \n"

	gets stdin prefab_psf

	set trimmed_filename [string trim $filename ".pdb"]
	# set trimmed_autopsf "$trimmed_filename\_autopsf"
		
	if {$prefab_psf == "no"} {

		# Catch for pre-existing autoPSF-generated PSF
		
		set trimmed_autopsf "$trimmed_filename\_autopsf"

		if {[file exists "$trimmed_autopsf\.psf"] == 1} {
			puts "A psf file of the pdb file to be converted already exists."
			puts "Continuing will result in an irreversible overwrite of such file."
			puts "Do you wish to proceed (yes/no)?"
			
			set move_on 0
			
			while {$move_on == 0} {
				gets stdin choice
				
				if {$choice == "yes"} {
					puts "Continuing with operation."
					# mol new $filename
					# set orig_ID [molinfo top]
					set move_on 1	
				} elseif {$choice == "no"} {
					puts "Terminating operation."
					return
				} else {
					puts "Your choice was neither yes nor no. Please re-enter:"
				}
			}
		} else {
			# mol new $filename
			# set orig_ID [molinfo top]
			set move_on 1
		}
		
		set orig_ID [molinfo top]

		# PSF Creation
		# To-do: create a user input prompt to allow for custom file prefix

		autopsf -mol $orig_ID -protein -regen -patch

		mol delete $orig_ID
		
	} else {
		
		# To-do: create a catch for invalid file name
		# mol new $filename
		set error_value [catch "mol addfile $prefab_psf" error_desc]
		if {$error_value == 1} {
			puts "The following error occurred while loading in the PSF file: \n"
			puts $error_desc
			puts "\nPlease try again."
		} else {
			set move_on 1
		}

	}
}

set aa_autopsf_ID [molinfo top]

# CG Model Creation

# To do: create a user input prompt to allow for parameter customization
# Implement a catch for pre-existing CG files

set numBeads 15
set inMolID [molinfo top]
set cgResName $::cggui::toCGshapeResName
set cgPrefix $::cggui::toCGshapeNamePrefix
set outPDB "cg_$trimmed_filename\.pdb"
set outAllAtom "aa_ref_$trimmed_filename\.pdb"
set outTop "cg_$trimmed_filename\.top"
set outParm "cg_$trimmed_filename\.par"
set numSteps [expr $numBeads * $::cggui::kStepsPerBead]
set epsInit $::cggui::toCGshapeEpsInit
set epsFinal $::cggui::toCGshapeEpsFinal
set lambdaInit [expr $numBeads * $::cggui::kLambdaMult]
set lambdaFinal $::cggui::toCGshapeLambdaFinal
set bondCutMethod 1
# set bondCutDist $::cggui::toCGshapeBondCutoff
set bondCutDist 85
set massValue -1

::cgnetworking::networkCGMolecule puts $inMolID $cgResName $cgPrefix $outPDB \
$outAllAtom $outTop $outParm $numBeads $numSteps $epsInit $epsFinal\
$lambdaInit $lambdaFinal $bondCutMethod $bondCutDist $massValue

mol delete $aa_autopsf_ID

mol new $outAllAtom
set aa_ref_ID [molinfo top]

# Extracting Dihedral and LJ Parameters
set out_LJ "[string trim $outParm ".par"]\_updated_LJ.par"
::cgtools_patch::sasa_LJ_networking puts $outParm $aa_ref_ID $out_LJ 150.0 1.0

mol delete $aa_ref_ID

set par_in [open $out_LJ r]
# set dihed_crude [split [read $par_in] \n]
set vdw_data [split [read $par_in] \n]
close $par_in

# set dihed_reached [lsearch $dihed_crude "DIHEDRALS"]
# set end_dihed [lsearch $dihed_crude "NONBONDED"]

# for {set i 0} {$i < [llength $dihed_crude]} {incr i} {
# lset dihed_crude $i [string trim [lindex $dihed_crude $i]]

# }

# set dihed_data {}

# for {set i $dihed_reached} {$i < $end_dihed} {incr i} {
	# lappend dihed_data [lindex $dihed_crude $i]
# }

# unset dihed_crude

for {set i 0} {$i < [llength $vdw_data]} {incr i} {
lset vdw_data $i [string trim [lindex $vdw_data $i]]

}

# vdw data will contain both dihedral and LJ parameters
set vdw_reached 0
while {$vdw_reached == 0} {
	if {[lindex $vdw_data 0] == "DIHEDRALS"} {
		set vdw_reached 1
	} else {
		set vdw_data [lreplace $vdw_data 0 0]
	}
}

mol new $outPDB
set cg_PDB_ID [molinfo top]

# To do: Implement a catch for pre-existing CG PSFs

autopsf -mol $cg_PDB_ID -top $outTop -regen -patch
set cg_PSF_ID [molinfo top]

mol delete $cg_PDB_ID

# Extracting Bond/Angle Params from All-Atom Simulation
# To do: create a user input prompt to allow for parameter customization
puts \n
puts "Enter the path to the equilibriation trajectory of the all-atom model: "
gets stdin dcd_ref

file mkdir Params

set aa_par_basename "from-aa-[string trim $outPDB ".pdb"]"

set psf_CG "[string trim $outPDB ".pdb"]\_autopsf.psf"
set pdb_CG "[string trim $outPDB ".pdb"]\_autopsf.pdb"
set pdb_ref $outAllAtom
# set dcd_ref "[string trim $filename ".pdb"]\.dcd"
set temperature 310.0
set f_out "Params/$aa_par_basename\.par"
set dat_bond "Params/from-aa-[string trim $outPDB ".pdb"]\-bonds.dat"
set dat_angle "Params/from-aa-[string trim $outPDB ".pdb"]\-angles.dat"

::cgnetworking::all_ba puts $psf_CG $pdb_CG $pdb_ref $dcd_ref $temperature $f_out $dat_bond $dat_angle

# Appending Dihedral and LJ parameters to extracted bond/angle param file

set from_aa_par [open $f_out a]
# puts $from_aa_par "\n"
# foreach line $dihed_data {
	# puts $from_aa_par $line
# }
puts $from_aa_par "\n"
foreach line $vdw_data {
	puts $from_aa_par $line
}
close $from_aa_par

# Map CG beta-field to CG bead numbers
mol new $pdb_CG

set A_all [atomselect top all]
set N [$A_all num]

for {set i 1} {$i <= $N} {incr i} {
  set A [atomselect top "serial $i"]
  $A set beta [expr $i * 1.0]
  $A delete
}

set pdb_bfield_CG "[string trim $pdb_CG ".pdb"]\-beta.pdb"

$A_all writepdb $pdb_bfield_CG
$A_all delete

mol delete all

# puts "Enter the amount of Boltzmann inversions you want to apply to the CG molecule: "
# gets stdin numIversions

set numInversions 4

file mkdir NAMD

set iter_template [open iteration1.conf r]
set iter_lines [split [read $iter_template] \n]
close $iter_template

for {set i 1} {$i <= $numInversions} {incr i} {
	puts "\nCreating NAMD iteration file directories"
	file mkdir "NAMD/iteration$i"
	file mkdir "NAMD/iteration$i\/input"
	file mkdir "NAMD/iteration$i\/output"
	file mkdir "Params/iteration$i"
	# file copy iteration1.conf "NAMD/iteration$i"
	
	puts {Copying pdb and psf CG files to NAMD input directory}
	file copy -force $psf_CG "NAMD/iteration$i\/input" 
	file copy -force $pdb_CG "NAMD/iteration$i\/input" 
	
	set trimmed_min "cg_$trimmed_filename\_minimization_$i"
	
	# Set up NAMD config file I/O parameters
	puts {Setting up NAMD config file I/O parameters}
	if {$i == 1} {
		set idxStruct [lsearch $iter_lines "structure 	input/cg_monomer-psfgen.psf"]
		lset iter_lines $idxStruct "structure 	input/$psf_CG"
		set idxCoords [lsearch $iter_lines "coordinates	input/cg_monomer-psfgen.pdb"]
		lset iter_lines $idxCoords "coordinates	input/$pdb_CG"
		set idxOutput [lsearch $iter_lines "set outputname     output/iteration1"]
		
		file copy -force $f_out NAMD/iteration1/input
		set idxParams [lsearch $iter_lines "parameters	    input/from-aa.par"]
		lset iter_lines $idxParams "parameters	    input/$aa_par_basename\.par"
		
		file copy -force $f_out "NAMD/iteration1/input" 
	} else {
		# Copy the scaled parameter file to the new iteration input folder
		lset iter_lines $idxParams "parameters	    input/$scale_name\.par"
		file copy -force $outParFilename "NAMD/iteration$i\/input" 
	}
	
	lset iter_lines $idxOutput "set outputname     output/$trimmed_min"
	puts {opening NAMD config file for newest iteration}
	set new_iter_conf [open "NAMD/iteration$i/iteration$i\.conf" w+]
	puts {writing to NAMD config file for newest iteration}
	foreach line $iter_lines {
		puts $new_iter_conf $line
	}
	close $new_iter_conf
	
	puts {executing NAMD simulation}
	exec namd2 "NAMD/iteration$i\/iteration$i\.conf"
	
	# set psf_CG "[string trim $outPDB ".pdb"]\_autopsf.psf"
	# set pdb_CG "[string trim $outPDB ".pdb"]\_autopsf.pdb"	
	set pdb_ref $pdb_bfield_CG
	set dcd_ref "NAMD/iteration$i\/output/$trimmed_min\.dcd"
	# set temperature 310.0
	set new_f_out "Params/iteration$i\/from-iter$i\-[string trim $outPDB ".pdb"]\.par"
	set dat_bond "Params/iteration$i\/from-iter$i\-[string trim $outPDB ".pdb"]\-bonds.dat"
	set dat_angle "Params/iteration$i\/from-iter$i\-[string trim $outPDB ".pdb"]\-angles.dat"

	::cgnetworking::all_ba puts $psf_CG $pdb_CG $pdb_ref $dcd_ref $temperature $new_f_out $dat_bond $dat_angle
	
	set scale_name "scaled-from-iter$i\-[string trim $outPDB ".pdb"]"
	
	# Set up Bond/Angle Spring Constant Scaling input parameters
	set bondScaleFactor  0.3
	set bondCutoff 3.5
	set angleScaleFactor 0.3 
	set angleCutoff 170.0
	set outParFilename "Params/iteration$i\/$scale_name\.par"

	if {$i == 1} {
		set inParFilename $f_out
	} else {
		set inParFilename "Params/iteration[expr $i - 1]\/from-iter[expr $i - 1]\-[string trim $outPDB ".pdb"]\.par"
	}
	
	::cgnetworking::scaleParameterConstants puts $inParFilename \
	$bondScaleFactor $bondCutoff $angleScaleFactor $angleCutoff $outParFilename
	
	# Append Dihedral and LJ Parameters from all-atom model to both CG scaled and extracted model
	# set scale_CG_par [open $outParFilename a]
	# # puts $scale_CG_par "\n"
	# # foreach line $dihed_data {
		# # puts $scale_CG_par $line
	# # }
	# puts $scale_CG_par "\n"
	# foreach line $vdw_data {
		# puts $scale_CG_par $line
	# }
	# close $scale_CG_par
	
	set from_CG_par [open $new_f_out a]
	# puts $from_CG_par "\n"
	# foreach line $dihed_data {
		# puts $from_CG_par $line
	# }
	puts $from_CG_par "\n"
	foreach line $vdw_data {
		puts $from_CG_par $line
	}
	close $from_CG_par
	
}