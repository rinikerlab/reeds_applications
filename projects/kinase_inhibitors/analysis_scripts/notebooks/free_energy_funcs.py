import reeds
from reeds.function_libs.file_management import file_management as fM

import os

import numpy as np
import pandas as pd

from scipy import constants as const
from scipy import special
from scipy.stats import pearsonr
from scipy.stats import kendalltau
from scipy.stats import spearmanr

import matplotlib.pyplot as plt
from matplotlib import cm

#active_qualitative_map_mligs = lambda num_ligs: plt.cm.viridis(np.linspace(0,1,num_ligs))
#active_qualitative_map_mligs = lambda num_ligs: plt.cm.cividis(np.linspace(0,1,num_ligs))
# color_blind_frienfly = active_qualitative_map_mligs(5)

def plot_free_energies(exp, sim, labels_sim, title=None, colors=None, markers=None, 
                       two_greyed_regions=True, relative=False, manual_lims=None):
    """
    This function will plot the results in a square scatter plot
    
    """
    
    active_qualitative_map_mligs = lambda num_ligs: plt.cm.viridis(np.linspace(0,1,num_ligs))
    active_qualitative_map_mligs = lambda num_ligs: plt.cm.cividis(np.linspace(0,1,num_ligs))
    color_blind_frienfly = active_qualitative_map_mligs(5)
    
    def plot_errorbar(ax, x, y, y_err, color, label, zorder, alpha, marker = 's'):
        markers, caps, bars = ax.errorbar(x, y, y_err, fmt=marker, markersize = 6, 
                                          color = color, mec = 'black', 
                                          ecolor = color, label = label, 
                                          capthick=1.5, capsize=4, linewidth = 3, 
                                          zorder=zorder, alpha = 1
                                         )
        [bar.set_alpha(alpha) for bar in bars]
        [cap.set_alpha(alpha) for cap in caps]
    
    xmin = np.min(exp)
    xmax = np.max(exp)
    
    fig, axes = plt.subplots(ncols=1, nrows=1, figsize =[9,9])
    
    if colors is None:
        colors = ['royalblue', 'firebrick', 'skyblue', 'darkorange']
    
    if markers is None:
        markers = ['D'] * len(labels_sim)
    
    if len(sim) == 2: # We we have gaff and openFF
        plot_errorbar(axes, exp, sim[0][0], sim[0][1], color = color_blind_frienfly[1],
                      label = labels_sim[0], zorder = 10, alpha = 1, marker = 's')
        
        plot_errorbar(axes, exp, sim[1][0], sim[1][1], color = color_blind_frienfly[-1],
              label = labels_sim[1], zorder = 10, alpha = 1, marker = 's')
    
    else:
        for i, (s, label) in enumerate(zip(sim, labels_sim)):
            plot_errorbar(axes, exp, s[0], s[1], color = colors[i],
                          label = label, zorder = 10, alpha = 1, marker = markers[i])
        
        if np.min(s[0]) < xmin:
            xmin = np.min(s[0])
        if np.max(s[0]) > xmax:
            xmax = np.max(s[0])
        
    leg = axes.legend(loc = 'best', fontsize=16, fancybox=True)
    leg.get_frame().set_edgecolor('black')
    
    x = np.linspace(-100, 100, 50)

    plt.plot(x, x, color = 'black')
    plt.plot(x, x-4.185, color = 'lightgrey')
    plt.plot(x, x+4.185, color = 'lightgrey')
    #
    plt.fill_between(x, x-4.185, x+4.185, color='lightgrey')
        
    if two_greyed_regions:
        plt.fill_between(x, x-2*4.185, x+2*4.185, color='lightgrey', alpha=0.5)

    if title is not None:
        plt.title(title, fontsize = 18)
    
    if relative:
        plt.xlabel(r'$\Delta\Delta G_{bind}$ - Experimental Values [kJ/mol]', fontsize = 18)
        plt.ylabel(r'$\Delta\Delta G_{bind}$ - RE-EDS Simulation [kJ/mol]', fontsize = 18)
    else: 
        plt.xlabel(r'$\Delta G_{bind}$ - Experimental [kJ/mol]', fontsize = 18)
        plt.ylabel(r'$\Delta G_{bind}$ - RE-EDS Simulation [kJ/mol]', fontsize = 18)
    
    # Automatically find the value in x closest to largest and smallest values
    plot_max = np.min(x[x>xmax])
    plot_min = np.max(x[x<xmin])
    
    if manual_lims is None:
        axes.set_ylim([plot_min, plot_max])
        axes.set_xlim([plot_min, plot_max])
    else:
        axes.set_ylim(manual_lims)
        axes.set_xlim(manual_lims)
    # Formatting
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(1)

    axes.tick_params(axis="y",direction="in", length = 6, width = 1, labelsize=18)
    axes.tick_params(axis="x",direction="in", length = 6, width = 1, labelsize=18)
    
    return fig, axes
            
#
# Data analysis
#

def direct_comparison(ic50s, alch_dg_complex, alch_dg_water):
    """
    This function performs a direct comparison where the simulated delta_delta_G values for 
    transformations 1->2, 1->3, ..., 1->N are compared to the experimental values (converted from IC_50s).
    
    Here this hence returns (N-1) data points corresponding to ddGs
    
    Args
    -----
        ic50s:
            list of ic50s
        alch_dg_complex:
            list of alchemical deltaGs (1->1, 1->2, ..., 1->N) (average and standard deviations) in complex
        alch_dg_water
            list of alchemical deltaGs (1->1, 1->2, ..., 1->N) (average and standard deviations) in water
    
    Returns
    -------
        ddgs_exp, ddgs_sim
    """
    
    # All numbers are going to be given w.r.t. the first state. 
    ddgs_exp = - (8.314 / 1000 * 298.15) * np.log(ic50s[0]/ic50s)
    ddgs_exp = ddgs_exp[1:]
    
    # Simulated vales
    ddg_binding = alch_dg_complex[0] - alch_dg_water[0]
    ddg_err = np.sqrt(np.power(alch_dg_complex[1], 2) + np.power(alch_dg_water[1], 2))
   
    ddgs_sim = np.array([ddg_binding[1:], ddg_err[1:]])

    # Here we don't return the obvious value for 1->1. 
    
    return (ddgs_exp, ddgs_sim)

def offset_comparison(ic50s, alch_dg_complex, alch_dg_water):
    """
    This function performs a direct comparison where the simulated delta_G values for 
    all states are calculated and offset by following the procedure described in:
    
    in eq. 2 of https://pubs.acs.org/doi/10.1021/acs.jctc.6b01141
     
    Here this hence returns N data points corresponding to dGs
    
    The error bars for state 1 are calculated as the average of all other errors. 
    
    
    Args
    -----
        ic50s:
            list of ic50s
        alch_dg_complex:
            list of alchemical deltaGs (1->1, 1->2, ..., 1->N) (average and standard deviations) in complex
        alch_dg_water
            list of alchemical deltaGs (1->1, 1->2, ..., 1->N) (average and standard deviations) in water
    
    Returns
    -------
        ddgs_exp, ddgs_sim
    """
    
    dgs_exp = (8.314 / 1000 * 298.15) * np.log(ic50s)
    
    # Simulated vales
    ddgs_binding = alch_dg_complex[0] - alch_dg_water[0]
    
    # Offsetting all data points 
    n = len(ddgs_binding)
    dgs_binding = ddgs_binding - (np.sum(ddgs_binding) / n - np.sum(dgs_exp) / n )
    
    # and propagating errors (water + complex) for other states.
    dgs_err = np.sqrt(np.power(alch_dg_complex[1], 2) + np.power(alch_dg_water[1], 2))
    
    # This sets an approximate error (when we don't use errors w.r.t. ref state )
    if dgs_err[0] == 0:
        dgs_err[0] = np.average(dgs_err[1:])

    dgs_sim = np.array([dgs_binding, dgs_err])
    return (dgs_exp, dgs_sim)
    
    


def additional_analysis(exp, sim, label):
    """
    Calculating MUEs and Pearson coefficients
    
    """
    
    print ('working on dataset for: ' + label + f' n={len(exp)}')
    
    mue = np.average(np.abs(exp-sim[0]))
    print (f'MUE calculated: {mue:.2f} kJ/mol' )
    print (f'MUE calculated: {mue/4.1815:.2f} kcal/mol' )
    
    corr, _ = pearsonr(exp, sim[0])
    print(f'Pearsons correlation: {corr:.2f}')
    
    tau, p_value = kendalltau(exp, sim[0])
    print(f'Kendall Tau: {tau:.2f}')
    
    rho, _ = spearmanr(exp, sim[0])
    print (f'Spearman rho: {rho:.2f}')
    
    print('\n')

    return mue, corr, tau   
    
    

#
# Calculating with all legs (older functions I am not really using anymore)
#

def make_complete_ddg_dict(ddg_binding):
    """
    Give as input a list of delta_delta_gs [1->1, 1->2, ..., 1->N]
    
    """
    # Put all of these vlaues in a dict 
    ddg_all_legs = {}
    for i, v in enumerate(ddg_binding):
        ddg_all_legs[f'{1}_{i+1}'] = v
    
    num_states = len(ddg_binding)
    for j in range(1, num_states):
        for k in range(j, num_states):
            # This is always the 1-> j and 1-> k free energy
            # to get the j -> k subtract them
            ddg_all_legs[f'{j+1}_{k+1}'] = ddg_binding[k] - ddg_binding[j]
        
    return ddg_all_legs
    

def calc_abs_ddG_all_paths(ic50s, alch_dg_complex, alch_dg_water):
    """
    This function will calculate the absolute binding energies for all N states
    by considering all paths (A->B, A->C, etc...)
    
    Write the equation here
    
    """
    
    num_states = len(ic50s)
    
    # Convert the ic50s to deltaGs 
    dgs_exp = (8.314 / 1000 * 298.15) * np.log(ic50s)
    
    # Simulated vales
    ddgs_binding = alch_dg_complex[0] - alch_dg_water[0]
    
    # Create the dict with all pair-wise combinations
    
    ddg_all_paths = make_complete_ddg_dict(ddgs_binding)
    
    sim_dGs = np.zeros(num_states)
    
    # Make a list of all pairs    
    for i in range(1, num_states+1):
        
        # Make a list of the partners of state i
        states = np.arange(1, num_states+1).tolist()
        
        # We will calculate it as  dG_x AVERAGE ( ddG_xy(path) + dG_exp_y)
        ddG_tmp = np.zeros(num_states)
        ddG_errs = np.zeros(num_states)
        
            
        for k, j in enumerate(states):
            #print ('working with pair: ' +str(i) + ' and ' + str(j))
            #print (ddG_tmp)
            
            # Find which value to use from the list
            if str(i) + '_' + str(j) in ddg_all_paths.keys():
                ddG_tmp[k] = ddg_all_paths[str(i) + '_' + str(j)] 
                #ddG_errs[k] = err[idx]
            else : # Here we will have to reverse the sign !!!!
                ddG_tmp[k] = - ddg_all_paths[str(j) + '_' + str(i)]
                #ddG_errs[k] = err[idx]
        
        #print(np.round(dgs_exp, 2))
        #print (np.round(ddG_tmp, 2))
        
        sim_dGs[i-1] = np.round(np.average(dgs_exp - ddG_tmp),2)  
    
    # For now just put 0 error for all of this
    
    sim_dGs_witherr = np.array([sim_dGs, np.zeros(num_states)])
    
    return (sim_dGs_witherr, dgs_exp)

#
# Extracting free energy values and errors from the gromos output files
#

def calc_dgs(energy_trajectory: pd.DataFrame, temp:float, trim_beg:float, num_states:int, sub_sample:int) -> np.array:

    """
    This function applies eqn. 6  of Sidler et al., J. Chem. Phys. 2016, 145, 154114 
    to estimate the energy offsets for a specific replica.
    Note that the original offset do not go into this equation!
    Note: We use special.logsumexp, (logsum and not log average) as 
          all states have the name number of data points, and log(1/N) cancels out.

    Parameters
    ----------
        energy_trajectory: pandas DataFrame 
            contains the potential energies of the end state and the ref. state
        temp: float
            temperature in Kelvin
        trim_beg: float
            fraction of the values to remove at the begining for "equilibration"
        num_states: int
            number of end states in our RE-EDS simulation
                
    Returns
    -------
        new_eoffs: np.array
            Energy offsets estimated for this replica scaled such that
            ligand 1 has an offset of 0. 
    """

    beta =  1000 * 1 / (temp * const.k * const.Avogadro)
    new_eoffs = np.zeros(num_states)

    initial_offsets = np.zeros(num_states)

    # note exp_term is a vector (each element is a timestep) 
    # containing the term to be exponentiated
    trim_beg = int(trim_beg*len(energy_trajectory['e1']))

    for i in range(num_states):
            v_i = np.array(energy_trajectory['e' + str(i+1)])[trim_beg::sub_sample]
            v_r = np.array(energy_trajectory['eR'])[trim_beg::sub_sample]
            exp_term = - beta * (v_i - v_r)
            new_eoffs[i] =  -(1/beta) * special.logsumexp(exp_term)
      
    return (new_eoffs - new_eoffs[0]), new_eoffs


def extract_alchemical_dgs(ene_trajs, num_states, out_path, verbose = True, sub_sample = 1, trim_beg = 0):
    """
    This function will calculate the free energies from the RE-EDS data of 
    n independent runs and save the results in a numpy array .npy file
    
    The array contains the average and standard dev. accross the 5 runs.
    
    """
    num_seeds = len(ene_trajs)
    
    ddgs = []
    dgs = [] 

    for i in range(num_seeds):
        dd_gs, d_gs = calc_dgs(ene_trajs[i], temp = 298.15, trim_beg=trim_beg, num_states = num_states, sub_sample = sub_sample)
    
        ddgs.append(dd_gs)
        dgs.append(d_gs)
    
        tmp = ""
        for v in dd_gs:
            tmp += " " + str(f'{v:.2f}')
        
        if verbose:
            print (tmp)
            print ('\n')

    # Saving results     
    # errors calcualted with the std of the dG (x -> R)
    arr = [np.average(ddgs, axis = 0), np.std(dgs, axis = 0)]

    if os.path.exists(out_path):
        print ('\n' + out_path + '  allready exits. We are not overwriting.\n')
    else:    
        np.save(out_path, arr)
    
    return arr


def bootstrap_analysis_replacement(sim_data, exp_data, with_noise=True, nboot=1000, ci_lvl = 0.95):
    """
    Do the bootstrap analysis with replacement
    """
    
    def print_ci_res(data, npoints):
        ci_low = np.mean(data) - ci_lvl * np.std(data) / np.sqrt(npoints)
        ci_up = np.mean(data) + ci_lvl * np.std(data) / np.sqrt(npoints)
        print (f'{np.mean(data):.2f} [{ci_low:.2f}, {ci_up:.2f}]')
        
        return None
        
    npoints = len(sim_data)
    
    mues = np.zeros(nboot)
    rhos = np.zeros(nboot)
    taus = np.zeros(nboot)
    
    for i in range(nboot):
        idx = np.random.randint( 0, npoints, npoints)
        # Resample 
        tmp_sim = sim_data[idx]
        tmp_exp = exp_data[idx]
        
        if with_noise:
            tmp_sim += np.random.normal(loc=0.0, scale=1, size=1)

        # Calc props        
        mues[i] = np.average(np.abs(tmp_exp-tmp_sim))
        rhos[i], _ = spearmanr(tmp_exp, tmp_sim)
        taus[i], _ = kendalltau(tmp_exp, tmp_sim)
        
    
    print_ci_res(mues, npoints)
    print_ci_res(rhos, npoints)
    print_ci_res(taus, npoints)
    
    return mues

# Make a list of all relative energies:

def get_rel_ddGs(dgs):
    """
    Convert a list of dg values into a list (linearized matrix)
    of all possible n * (n-1) relative values (excluding self terms)
    
    """
    n = len(dgs)
    num_rel = int(n*(n-1)/2)
    ddgs = np.zeros(num_rel)
    
    beg = 0
    for i in range(n-1):
        delta = (n - (i+1))
        end = beg + delta
        ddgs[beg:end] = (dgs - dgs[i])[i+1:]
        beg = end
        
    return np.array(ddgs)

def bootstrap_analysis_fromseed(sim_data, sim_err, exp_data, nboot=1000, ci_lvl = 0.95):
    """
    Do the bootstrap analysis by adding noise to the data from the standard deviations 
    propagated from the simulations
    """

    def print_ci_res(data, npoints, varname):

        ci_low = np.mean(data) - ci_lvl * np.std(data) / np.sqrt(npoints)
        ci_up = np.mean(data) + ci_lvl * np.std(data) / np.sqrt(npoints)
        print (f'{varname} = {np.mean(data):.1f} [{ci_low:.2f} -- {ci_up:.2f}]')

        return None

    npoints = len(sim_data)

    mues = np.zeros(nboot)
    rhos = np.zeros(nboot)
    taus = np.zeros(nboot)

    rmses = np.zeros(nboot)
    
    ddgs_exp = get_rel_ddGs(exp_data)
    
    for i in range(nboot):
        # Add noise to data 
        tmp_sim = sim_data + np.random.normal(loc=0, scale=sim_err, size = len(sim_err))
        
        ddgs_tmp = get_rel_ddGs(tmp_sim)
        
        # Calc props        
        mues[i] = np.average(np.abs(exp_data-tmp_sim))
        rhos[i], _ = spearmanr(exp_data, tmp_sim)
        taus[i], _ = kendalltau(exp_data, tmp_sim)

        rmses[i] = np.sqrt(np.mean(np.power(ddgs_tmp-ddgs_exp, 2)))
        
        
    print_ci_res(mues, npoints=npoints, varname= 'MUE')
    print_ci_res(taus, npoints=npoints, varname= r'$\tau$')
    print_ci_res(rhos, npoints=npoints, varname= r'$\rho$')
    print_ci_res(rmses, npoints=npoints, varname= r'RMSE')
    
    return None
