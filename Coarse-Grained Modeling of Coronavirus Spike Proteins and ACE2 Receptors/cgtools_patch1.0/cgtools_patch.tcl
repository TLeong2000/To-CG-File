# sasa_LJ_patch is a modification of cgtools.tcl intended to fix memory leak
# issues in the sasa_LJ_networking procedure.
#
# This patch is developed in accordance with the VMD Plugin Library License.
#
# cgtools is the work of Anton Arkhipov, Peter L. Freddolino, and Klaus Schulten
#
# sasa_LJ_patch is the work of Timothy Leong
#
package provide cgtools_patch 1.0

namespace eval ::cgtools_patch:: {
	namespace export sasa_LJ_networking
}

# -----------------------------------------------------
# This proc assigns Lennard-Jones (LJ) parameters for a
# coarse-grained (CG) structure based on the all-atom one.
# It extracts solvent accessible surface area (SASA)
# for each atomic domain representing a CG bead from an
# all-atom structure corresponding to the CG structure.
# The values extracted are SASA for the whole domain
# (in the context of the rest of the structure) and SASA for
# hydrophobic residues only. These values are used to assign the
# LJ well depths E_i to individual CG beads, as
# E_i = ELJ * (SASA_i[hydrophobic]/SASA_i[total])^2, where
# ELJ is an adjustable constant.
# Radius of gyration is used to compute LJ radii; LJ radius R_i is
# obtained as R_i = r_gyr_i + RLJ, where r_gyr_i is
# the radius of gyration for all atoms represented by the i-th CG bead
# and RLJ is an adjustable constant.
##################################################
# INPUT DATA:
#
# CG parameter file, "par_CG";
#
# ID ("pdbrefID") of the all-atom reference structure that should be
# loaded in VMD already; this structure should have beta values
# filled with numbers corresponding to the CG beads' numbers;
# since radius of gyration is computed, it is better to load
# the PSF file into pdbrefID too, so that correct masses are used;
#
# the file to where the output is written, "f_out";
#
# maximum energy value for the LJ well depth, "ELJ" (kcal/mol);
#
# an addition to the LJ radius RLJ (A).
##################################################
proc ::cgtools_patch::sasa_LJ_networking {statusProc par_CG pdbrefID f_out ELJ RLJ} {

# Check if we can write the output file.
set outfile [open $f_out w]
close $outfile

# # Find out number of CG beads.
# # puts "Checking the number of CG beads..."
# set NAA [[atomselect $pdbrefID all] num]
# set NCG 0
# for {set i 0} {$i < $NAA} {incr i} {
  # set Ntmp [expr int([[atomselect $pdbrefID "index $i"] get beta])]
  # if {$Ntmp > $NCG} {
    # set NCG $Ntmp
  # }
# }
# # puts "done."
# # puts "According to the reference molecule (ID $pdbrefID) NCG = $NCG."
# # puts ""

##################################################
# this section is done by Tim Leong b/c of memory
# issues concering the above CG bead check
##################################################

# Using only one atom selection prevents the memory from being clogged up
# by numerous undeleted object pointers (atomselections)
set beta_list [[atomselect $pdbrefID all] get beta]
set NCG 0
# jdex is for debugging purposes
# set jdex -1
foreach {beta_val} $beta_list {
	# incr jdex
	if {$beta_val > $NCG} {
		set NCG $beta_val
	}
	# puts $jdex
}
unset beta_list
# There should not be any adverse effects from leaving the integer conversion
# until the very end. The benefit of this is decreased computational load.
set NCG [expr int($NCG)]

##################################################

# Read the CG parameter file and extract LJ parameters.
##################################################
set fdata ""
set par_CG_file [open $par_CG r];

set outfile [open $f_out w]
# Find where non-bonded entries start.
gets $par_CG_file fdata
set tmp [lindex $fdata 0]
while {$tmp != "NONBONDED"} {
  puts $outfile "$fdata"
  gets $par_CG_file fdata
  set tmp [lindex $fdata 0]
}

# Read current LJ parameters
# (basically, this is done to initialize arrays.).
set k 0
set i 0
while {($k < $NCG) && ($i < 100000000)} {
  incr i
  gets $par_CG_file fdata
  if {[string range [lindex $fdata 0] 0 0] != "!"} {
    incr k
    set CG($k.name) [lindex $fdata 0]
    set CG($k.E) [lindex $fdata 2]
    set CG($k.r) [lindex $fdata 3]
  }
}

##################################################


##################################################
# Loop over all CG beads; find SASA and r_gyr.   #
##################################################
set A_all [atomselect $pdbrefID all]
for {set k 1} {$k <= $NCG} {incr k} {
  # puts "Bead $k of $NCG"
  set A [atomselect $pdbrefID "beta $k"]
  set A_hphob [atomselect $pdbrefID "beta $k and hydrophobic"]
  set tmp [measure sasa 0.0 $A_all -restrict $A]
  set CG($k.sasa) [measure sasa 0.0 $A_all -restrict $A]
  set CG($k.sasa_hphob) [measure sasa 0.0 $A_all -restrict $A_hphob]
  
  # r_gyr
  if {[$A num] > 0} {
    set CG($k.r) [expr [measure rgyr $A] + $RLJ]
  }
  
  $A delete
  $A_hphob delete
  if {$CG($k.sasa) == 0} {
    set CG($k.sasa_ratio) 0.0
  } else {
    set CG($k.sasa_ratio) [expr $CG($k.sasa_hphob)/$CG($k.sasa)]
  }
}
$A_all delete
##################################################
##################################################


##################################################
# Append the output to a file.                  #
##################################################

puts $outfile "NONBONDED"
puts $outfile "!"
puts $outfile "!V(Lennard-Jones) = Eps,i,j\[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6\] "
puts $outfile "!"
puts $outfile "!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j) "
puts $outfile "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j "
puts $outfile "!"
puts $outfile "!atom  ignored    epsilon      Rmin/2 "
puts $outfile "!"

##################################################
for {set k 1} {$k <= $NCG} {incr k} {
  puts $outfile [format "%-5s%10f%11.6f%12.6f" \
                            $CG($k.name) 0.0 \
			    [expr (-1.0)*$ELJ*$CG($k.sasa_ratio)*$CG($k.sasa_ratio)] \
			    $CG($k.r)]
}
##################################################
puts $outfile ""
puts $outfile "END"
puts $outfile ""

close $outfile
##################################################
##################################################

}