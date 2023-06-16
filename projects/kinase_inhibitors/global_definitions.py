import os
import reeds
from pygromos.euler_submissions import FileManager as fM
from pygromos.utils import bash

name = "NAME_OF_SYSTEM" # example "PIM_complex_openff_subset2"
root_dir = os.getcwd()

# 
# change the following paths to correspond to the system with which you will be working
#

input_folder =    root_dir+"/0_input"
in_cnf_file =     input_folder+"/PIM_subset2_hybrid_complex.cnf"  
in_top_file =     input_folder+"/PIM_subset2_openff_hybrid_complex.top"
in_pert_file =    input_folder+"/PIM_subset2_openff_hybrid.ptp"


in_disres_file =  None 

gromosXX_bin = "/path/to/gromosXX/bin/"
gromosPP_bin = "/path/to/gromos++/bin/"
ene_ana_lib =  input_folder +"/ene_ana.md++.lib"

in_template_eds_imd = input_folder+"/template_md.imd"

undersampling_frac_thresh = 0.90
