#!/usr/bin/env python3
import os, sys, glob
sys.path.append(os.getcwd())

from global_definitions import fM, bash
from global_definitions import name, root_dir
from global_definitions import gromosXX_bin, gromosPP_bin, ene_ana_lib
from global_definitions import in_top_file, in_pert_file, in_disres_file
from global_definitions import undersampling_frac_thresh

from reeds.modules import do_RE_EDS_production as production

for i in range(1, 6):
            
    #Paths
    in_name = name+"_prod_seed" + str(i)
    next_production_dir = root_dir+"/path/to/d_mixed_opt/mixedOptX/analysis/next" #CHANGE HERE
    out_production_dir = root_dir+"/e_prod_seed" + str(i)
    optimized_states_dir = root_dir + "/a_optimizedStates/analysis/next"
    lower_bound_dir = root_dir + "/b_lowerBound/analysis/next"
     
    ##make folder
    out_production_dir = bash.make_folder(out_production_dir)

    #In-Files
    topology = fM.Topology(top_path=in_top_file, disres_path=in_disres_file, perturbation_path=in_pert_file)

    coords = glob.glob(next_production_dir+"/coord/*cnf") 

    in_template_reeds_imd = glob.glob(next_production_dir+"/*imd")[0]

    system = fM.System(coordinates=coords, name=in_name, top=topology)

    #Additional Options
    ## Simulation Params

    # We want 5ns total with 1ns equil. (set according to length desired)

    num_equilibration_runs = 2
    num_simulation_runs= 8

    job_duration="24:00:00"
    nmpi_per_replica = 8
    memory = 1000

    randomize = True

    #Do:
    last_jobID = production.do(out_root_dir=out_production_dir,
                                in_simSystem=system,
                                in_template_imd=in_template_reeds_imd,
                                gromosXX_bin_dir = gromosXX_bin,
                                gromosPP_bin_dir = gromosPP_bin,
                                in_ene_ana_lib_path = ene_ana_lib,
                                undersampling_fraction_threshold = undersampling_frac_thresh,
                                num_equilibration_runs = num_equilibration_runs,
                                num_simulation_runs = num_simulation_runs,
                                nmpi_per_replica = nmpi_per_replica,
                                memory = memory,
                                optimized_states_dir = optimized_states_dir,
                                lower_bound_dir = lower_bound_dir,
                                duration_per_job = job_duration, 
                                randomize = randomize, 
                                submit = True)


