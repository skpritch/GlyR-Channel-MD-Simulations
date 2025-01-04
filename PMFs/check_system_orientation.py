import numpy as np
import MDAnalysis
import MDAnalysis.analysis.rms
import nglview as nv
import matplotlib.pyplot as plt
import pandas
import seaborn as sns
import math
import pickle     # save and load python variables
import time


# load libraries:
import json                             # read in JSON files
import numpy as np                      # manipulate numeric vectors
from matplotlib import pyplot as pl     # plotting facilities
import argparse                         # parse command line arguments

import sys
import os
path='/biggin/b198/orie4254/Documents/scripts/pore_analysis/'
#'/biggin/b198/orie4254/Documents/possible_hydrophobic_gating_problem/'
res2 = os.listdir(path)
#print(res2)
sys.path.append(path)
import plot_chap as plt_chap
import hole_analysis as hole_analysis

from MDAnalysis.analysis import hole2
from MDAnalysis.analysis import align
#%matplotlib inline
hole_exe = '/biggin/b198/orie4254/hole2-ApacheLicense-2.2.005-Linux-x86_64/hole2/exe/hole'
sph_proc = '/biggin/b198/orie4254/hole2-ApacheLicense-2.2.005-Linux-x86_64/hole2/exe/sph_process'
f_size = 18

def write_out_info(system, win, ion, print_resids=False, Drude=False):
    p_GlyR = '/biggin/b232/temp-pritchard/Documents/PMFs/'
    p = p_GlyR + system
    if Drude:
        top = p + 'step3_charmm2omm.psf'
        conf = p + 'PROD/conf'+str(win)+'.pdb'
        u = MDAnalysis.Universe(top, conf, topology_format='psf', format='pdb',)
    else:
        top = p + 'PROD_CLA/PROD_'+str(win)+'.tpr'
        conf = p + 'PROD_CLA/PROD_'+str(win)+'.gro'
        u = MDAnalysis.Universe(top, conf, topology_format='tpr', format='gro',
                                tpr_resid_from_one=True)

    protein = u.select_atoms('protein')
    segIDs = np.unique(protein.segids)
    print('segids of protein', segIDs)
    COM = protein.center_of_mass()
    COG = protein.center_of_geometry()
    print('COM and COG = ',COM, COG)
    
    backbone = u.select_atoms('backbone')
    lipids = u.select_atoms('resname POPC')
    lipids_phosp = u.select_atoms("group acidic and name P", acidic=lipids, updating=True)

    CA_membrane_sys = u.select_atoms('name CA')
    chain1 = u.select_atoms('name CA and segid '+segIDs[0])
    chain2 = u.select_atoms('name CA and segid '+segIDs[1])
    chain3 = u.select_atoms('name CA and segid '+segIDs[2])
    chain4 = u.select_atoms('name CA and segid '+segIDs[3])
    chain5 = u.select_atoms('name CA and segid '+segIDs[4])
    Leu_chain1 = u.select_atoms('resname LEU and segid seg_0_PROA')
    Pro_chain1 = u.select_atoms('resname PRO and segid seg_0_PROA')
    Pm2_gate = []
    L9_gate = []

    for i in range(len(CA_membrane_sys)):
        if i>3 and i<len(CA_membrane_sys)-5:
            r0=CA_membrane_sys[i].resname
            r1=CA_membrane_sys[i+1].resname
            r2=CA_membrane_sys[i+2].resname
            r3=CA_membrane_sys[i+3].resname
            r4=CA_membrane_sys[i+4].resname
            rm1=CA_membrane_sys[i-1].resname
            rm2=CA_membrane_sys[i-2].resname
            rm3=CA_membrane_sys[i-3].resname
            ### test for A-A- P -A-R sequence which is P-2' gate ###
            if rm2=='ALA' and rm1=='ALA' and r0=='PRO' and r1=='ALA' and r2=='ARG':
                Pm2_gate.append(CA_membrane_sys[i].resid)
            ### test for T-T-V- L -T-M-T-T sequence which i sL9' gate ###
            if rm2=='THR' and rm1=='VAL' and r0=='LEU' and r1=='THR' and r2=='MET':
                L9_gate.append(CA_membrane_sys[i].resid)
    if print_resids:
        for i in range(len(chain1)):
            print(chain1[i].resid, chain1[i].resname,
                  chain2[i].resid, chain2[i].resname,
                  chain3[i].resid, chain3[i].resname,
                  chain4[i].resid, chain4[i].resname,
                  chain5[i].resid, chain5[i].resname)
    print('resids of P-2 gate', Pm2_gate)
    print('resids of L9 gate', L9_gate)
    sel_Pm2_gate = 'protein and resname PRO and (resid ' + str(Pm2_gate[0]) 
    sel_L9_gate = 'protein and resname LEU and (resid ' + str(L9_gate[0]) 
    for i in range(1, len(Pm2_gate)):
        sel_Pm2_gate = sel_Pm2_gate + ' or resid ' + str(Pm2_gate[i])
        sel_L9_gate = sel_L9_gate + ' or resid ' + str(L9_gate[i])
    sel_Pm2_gate = sel_Pm2_gate + ' )'
    sel_L9_gate = sel_L9_gate + ' )'
    #print(sel_Pm2_gate, sel_L9_gate)
    
    Pro_gate = u.select_atoms(sel_Pm2_gate)
    Leu_gate = u.select_atoms(sel_L9_gate)
    
    z_Pro_gate = []
    for a in Pro_gate:
        z_Pro_gate.append(a.position[2])
    print('P-2 gate z-cooridnate mean and std', np.mean(z_Pro_gate), np.std(z_Pro_gate))
    z_Leu_gate = []
    for a in Leu_gate:
        z_Leu_gate.append(a.position[2])
    print('L9 gate z-cooridnate mean and std', np.mean(z_Leu_gate), np.std(z_Leu_gate))
    
    CL = u.select_atoms('resname CLA')
    SOD = u.select_atoms('resname SOD')
    #print('CL', len(CL))

    # get transmembrane domain
    extracellular = []
    intracellular = []

    BOUNDARY = COM[2]
    for i in lipids_phosp:
        z = i.position[2]
        if z<BOUNDARY:
            extracellular.append(z)
        else:
            intracellular.append(z)
    com_z = protein.center_of_mass()[2]
    cog_z = protein.center_of_geometry()[2]
    leaflets = np.array([ np.mean(intracellular), np.mean(extracellular) ] )
    leaflet_std = np.array([np.std(intracellular), np.std(extracellular)] )
    inds = np.argsort(leaflets)
    leaflets = leaflets[inds]
    leaflet_std = leaflet_std[inds]
    print('membrane leaflet1', leaflets[0], leaflet_std[0])
    print('membrane leaflet2', leaflets[1], leaflet_std[1])
    print('com_z', com_z)
    print('cog_z', cog_z)
    
    print()
    print('relative to COM (position-COM)')
    print('Leu9 gate',np.mean(z_Leu_gate)-com_z, 'Pro-2 gate',np.mean(z_Pro_gate)-com_z)
    print('membrane boundaries', np.mean(intracellular)-com_z, np.mean(extracellular)-com_z)
    ion_PMF = u.select_atoms('bynum '+str(ion))
    ion_pos = ion_PMF.positions[0][2]
    print('ion_PMF', ion, ion_PMF.names, 'pos_z', ion_pos, 'relative to COM', ion_pos-com_z )



#write_out_info(system='6ply_GlyR_pore/', win=13, ion=32081)
write_out_info(system='6pm2_GlyR_pore/', win=-13, ion=28728) # SP
#write_out_info(system='7m6q_GlyR_pore/', win=13, ion=33617)
#write_out_info(system='7m6s_GlyR_pore/', win=-13, ion=33300)
#write_out_info(system='6ply_GlyR_Drude/', win=13, ion=48892, Drude=True)