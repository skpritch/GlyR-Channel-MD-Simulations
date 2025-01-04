import numpy as np
import MDAnalysis
import MDAnalysis.analysis.rms
from MDAnalysis.lib.distances import capped_distance
import nglview as nv
import matplotlib.pyplot as plt
import pandas
import seaborn as sns
import math
import pickle     # save and load python variables
import time

# Load libraries:
import json                             # read in JSON files
import argparse                         # parse command line arguments

import sys
import os

# If you have custom modules, adjust the path accordingly
# path = '/biggin/b198/orie4254/Documents/scripts/pore_analysis/'
# sys.path.append(path)
# import plot_chap as plt_chap
# import hole_analysis as hole_analysis

# Define cutoff distances for ions (values in nm)
CLA_radial_dist = {
    'water_oxygen': 3.75,  # Å
    'water_hydrogen': 2.75,  # Å
    'protein': 3.0  # Å
}

SOD_radial_dist = {
    'water_oxygen': 3.2,  # Å
    'water_hydrogen': 2.3,  # Å
    'protein': 2.7  # Å
}

POT_radial_dist = {
    'water_oxygen': 3.4,  # Å
    'water_hydrogen': 2.5,  # Å
    'protein': 2.85
}


def PMF_analysis(p, vec=np.arange(-45,46), higher=[], load_data=True,
                 fname='PROD_', title='PMF', fname_higher='', fac=1,
                 prod_dir='david_PROD_CLA', wham_dir='david_WHAM_CLA_unclamp',
                  cutoff_distances=POT_radial_dist):
    f_size = 14  # Font size for plots

    if load_data:
        summary = np.array([])
        summary_add = np.array([])
        print('Number of windows:', len(vec))
        for count, i in enumerate(vec):
            try:
                pullx_file = os.path.join(p, prod_dir, fname + str(i) + '_pullx.xvg')
                d = np.loadtxt(pullx_file, skiprows=17)
                summary = np.append(summary, d[:,1])
                summary_add = np.append(summary_add, d[:,1])
                print(f"Window {i}: Data points = {len(d[:,1])}, Mean = {np.mean(d[:,1]):.3f}, Std = {np.std(d[:,1]):.3f}")
                if abs(np.mean(d[:,1]) - i / 10) > 0.05:
                    print(f"WARNING: Window {i} mean deviation exceeds threshold!")
            except Exception as e:
                print(f"Window {i} not accessible: {e}")

        for i in higher:
            pullx_file = os.path.join(p, prod_dir, fname_higher + str(i) + '_pullx.xvg')
            d = np.loadtxt(pullx_file, skiprows=17)
            summary_add = np.append(summary_add, d[:,1])

        print('Total data points in summary:', len(summary))
        print('Data range:', min(summary), 'to', max(summary))

        # Plot histogram of pullx data
        NBins = 100
        fig, ax = plt.subplots()
        plt.title(title, fontsize=f_size)
        ax.set_ylabel('PDF', fontsize=f_size)
        ax.set_xlabel('COM distance [nm]', fontsize=f_size)
        N3, edges3 = np.histogram(summary, density=True, bins=NBins)
        ax.plot(edges3[1:], N3, '-', label='Pulling Data')
        ax.set_ylim([0, 1.1 * max(N3)])
        ax.legend(loc="best", prop={'size': 12})
        ax.tick_params(axis='both', which='major', labelsize=f_size)
        plt.tight_layout()
        histogram_file = os.path.join(p, wham_dir, f'single_histogram_{title}.png')
        fig.savefig(histogram_file, bbox_inches='tight')
        plt.show()

        # Write files for gmx wham
        os.makedirs(os.path.join(p, wham_dir), exist_ok=True)

        tpr_file_list = os.path.join(p, wham_dir, "files_tpr.dat")
        with open(tpr_file_list, "w") as f:
            for i in vec:
                tpr_file = os.path.join(p, prod_dir, f'PROD_{i}.tpr')
                f.write(tpr_file + '\n')
            for i in higher:
                tpr_file = os.path.join(p, prod_dir, f'{fname_higher}{i}.tpr')
                f.write(tpr_file + '\n')

        pullx_file_list = os.path.join(p, wham_dir, "files_pullx.dat")
        with open(pullx_file_list, "w") as f:
            for i in vec:
                pullx_file = os.path.join(p, prod_dir, f'{fname}{i}_pullx.xvg')
                f.write(pullx_file + '\n')
            for i in higher:
                pullx_file = os.path.join(p, prod_dir, f'{fname_higher}{i}_pullx.xvg')
                f.write(pullx_file + '\n')

    # Plot PMF profiles
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(title, fontsize=f_size)
    for i in range(5, 22, 2):
        pmf_file = os.path.join(p, wham_dir, f'profile_{i}ns.xvg')
        try:
            pmf = np.loadtxt(pmf_file, skiprows=17)
            ax.plot(pmf[:,0], pmf[:,1], '--', label=f'{i} ns per window')
        except Exception as e:
            print(f"Could not read {pmf_file}: {e}")

    ax.set_xlabel("Reaction coordinate (nm)", fontsize=f_size)
    ax.set_ylabel("PMF (kcal/mol)", fontsize=f_size)
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    ax.legend()
    plt.tight_layout()
    pmf_plot_file = os.path.join(p, wham_dir, f"PMF_{title}.png")
    fig.savefig(pmf_plot_file, bbox_inches='tight')

    # Plot PMF for specific time ranges
    digit = 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(title, fontsize=f_size)
    for i in [7, 14, 21]:
        pmf_file = os.path.join(p, wham_dir, f'profile_{i-7}-{i}ns.xvg')
        try:
            pmf = np.loadtxt(pmf_file, skiprows=17)
            offset = pmf[-1,1]
            ind_min = np.argmin(pmf[:,1]-offset)
            ind_max = np.argmax(pmf[:,1]-offset)
            max_str = f" PMF({pmf[ind_max,0]:.{digit}f} nm) = {pmf[ind_max,1]-offset:.{digit}f} kcal/mol"
            label_str = f" PMF({pmf[ind_min,0]:.{digit}f} nm) = {pmf[ind_min,1]-offset:.{digit}f} kcal/mol\n{max_str}"
            ax.plot(pmf[:,0], pmf[:,1]-offset, '-', label=f'{i-7}-{i} ns {label_str}')
        except Exception as e:
            print(f"Could not read {pmf_file}: {e}")

    ax.legend(loc="best")
    ax.set_xlabel("Reaction coordinate (nm)", fontsize=f_size)
    ax.set_ylabel("PMF (kcal/mol)", fontsize=f_size)
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    plt.tight_layout()
    pmf_repeats_file = os.path.join(p, wham_dir, f"PMF_repeats_{title}.png")
    fig.savefig(pmf_repeats_file, bbox_inches='tight')
    plt.show()

    # Average PMF over repeats
    pmf_vec = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(title, fontsize=f_size)
    for i in [7, 14, 21]:
        pmf_file = os.path.join(p, wham_dir, f'profile_{i-7}-{i}ns.xvg')
        try:
            pmf = np.loadtxt(pmf_file, skiprows=17)
            offset = pmf[-1,1]
            pmf_vec.append(pmf[:,1]-offset)
        except Exception as e:
            print(f"Could not read {pmf_file}: {e}")

    if pmf_vec:
        pmf_vec = np.array(pmf_vec)
        pmf_av = np.mean(pmf_vec, axis=0)
        pmf_std = np.std(pmf_vec, axis=0)
        ind = np.argmax(pmf_av)
        label_str = f" PMF({pmf[ind,0]*fac:.{digit}f} nm) = {pmf_av[ind]:.{digit}f} ± {pmf_std[ind]:.{digit}f} kcal/mol"
        ax.plot(pmf[:,0]*fac, pmf_av, label=label_str)
        ax.fill_between(pmf[:,0]*fac, pmf_av - pmf_std, pmf_av + pmf_std, facecolor='grey', alpha=0.3)

    ax.legend(fontsize=12, frameon=False)
    ax.set_xlabel("Position relative to COM (nm)", fontsize=f_size)
    ax.set_ylabel("PMF (kcal/mol)", fontsize=f_size)
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    plt.tight_layout()
    pmf_avg_file = os.path.join(p, wham_dir, f"PMF_repeats_{title}_av_std.png")
    fig.savefig(pmf_avg_file, bbox_inches='tight')
    plt.show()

    

# Define Input Directory
path_GlyR = '/biggin/b232/temp-pritchard/Documents/PMFs/'
p = os.path.join(path_GlyR, '6pm2_GlyR_Pore')

PMF_analysis(p, vec=np.arange(-45,46), higher=[], load_data=False,
             fname='PROD_', title='PMF 6pm2 CLA', fname_higher='',
             prod_dir='PROD_CLA', wham_dir='david_WHAM_CLA_unclamp', 
             cutoff_distances=POT_radial_dist)
