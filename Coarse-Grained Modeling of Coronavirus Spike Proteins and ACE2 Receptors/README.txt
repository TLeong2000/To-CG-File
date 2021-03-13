Place your protein pdb file and its equilibriation trajectory
in the same directory as the one containing To-CG-File.tcl,
iteration1.conf. Make sure that if you are using a protein pdb
file that has non-protein molecules in it (e.g. solvent), you
have a protein-only version of the equilibriation trajectory
available. The current version of To-CG-File does not have the
functionality of stripping non-protein molecules from an
equilibriation trajectory currently implemented.

The cgtools_patch1.0 folder should be placed in your default
VMD library directory OR should be placed in the same directory
as the one containing To-CG-File.tcl and iteration1.conf

Always make backups of your files before you this script, as
this script will result in the overwrite of some files/folders
located within the same directory of this script.

To-CG-File is the work of Timothy Leong under the instruction
of Zhangli Peng, Ph.D.