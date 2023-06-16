#!/usr/bin/env python3
import os, sys, glob

import reeds
import pygromos

sys.path.append(os.getcwd())
from global_definitions import *

#spefici parts
from reeds.modules import do_RE_EDS_findLowerBound as findLowerBound

out_lowerBound_dir = root_dir+"/b_lowerBound"
in_name = name+"_find_lower_bound"

if(not os.path.exists(out_lowerBound_dir)):
    os.mkdir(out_lowerBound_dir)

topology_lower_bound = fM.Topology(top_path=in_top_file, disres_path=in_disres_file, perturbation_path=in_pert_file)
system_lower_bound = fM.System(coordinates=in_cnf_file, name=in_name,top=topology_lower_bound)
print(system_lower_bound)

# Here I will just run the default time

findLowerBound.do(in_simSystem=system_lower_bound,template_imd=in_template_eds_imd,
out_root_dir=out_lowerBound_dir,gromosXX_bin=gromosXX_bin, 
              ene_ana_lib = ene_ana_lib, 
              job_duration= "24:00", 
              submit = True
              )


