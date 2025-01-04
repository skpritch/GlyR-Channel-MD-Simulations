### libraries
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
res2 = os.listdir(path)
#print(res2)
sys.path.append(path)
import plot_chap as plt_chap
import hole_analysis as hole_analysis

folder = '6pm2_GlyR_Pore/'
p = '/biggin/b232/temp-pritchard/Documents/PMFs/'+folder

### last frame of equilibration run 1,000ns ###
top = p + 'cg_output/EQUIL_POT/eq5.tpr'
conf = p + 'cg_output/EQUIL_POT/eq5.gro'
u = MDAnalysis.Universe(top, conf, topology_format='tpr', format='gro', tpr_resid_from_one=True)
print('number of atoms', len(u.atoms))

protein = u.select_atoms('protein')
backbone = u.select_atoms('backbone')
lipids = u.select_atoms('resname POPC')
lipids_phosp = u.select_atoms("group acidic and name P", acidic=lipids, updating=True)
print('number of POPC molecules', len(lipids_phosp))
water = u.select_atoms('resname TIP3')
print('water atoms', len(water), 'moleculs', len(water)/3)

CA_membrane_sys = u.select_atoms('name CA')
CL = u.select_atoms('resname CLA')
K = u.select_atoms('resname POT')
print('CL', len(CL), 'K', len(K))

COM = protein.center_of_mass()
COG = protein.center_of_geometry()
print('COM and COG = ',COM, COG)

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
print('intracellular', np.mean(intracellular), np.std(intracellular))
print('extracellular', np.mean(extracellular), np.std(extracellular))
print('com_z', com_z)
print('cog_z', cog_z)

print('PBC atom')
#print(CA_membrane_sys[30])
#pbc_atom = u.select_atoms('resname LEU and resid 31 and name CA') #
#print(pbc_atom.positions, pbc_atom.types, pbc_atom.indices)
#print()

ion = u.select_atoms('bynum 28728') # 6pm2 CLA SP
# ion = u.select_atoms('bynum 28687') # 6pm2 POT SP
# ion = u.select_atoms('bynum 28687') # 6pm2 SOD SP

print('ion', ion.positions, ion.resnames)
print(CL[0].index)
print(CL[0].position, ion.positions)

Cl_within_pore = u.select_atoms("group acidic  and prop z > 55 and prop z < 93", acidic=CL, updating=True)
#print(Cl_within_pore.indices)
#Cl0 = u.select_atoms('bynum 78170') #78158
#print(Cl0.positions, Cl0.indices)
# bynum = indices + 1
thr = 5
prop_x = 'prop x > '+str(COM[0]-thr)+' and prop x < '+str(COM[0]+thr)+' and '
prop_y = 'prop y > '+str(COM[1]-thr)+' and prop y < '+str(COM[1]+thr)+' and '
prop_z = 'prop z > '+str(COM[2]-50)+' and prop z < '+str(COM[2]+50)+'' #z > 100 
print(prop_x, prop_y, prop_z)
water_reactionCoordinate = u.select_atoms('name OH2 and '+prop_x+prop_y+prop_z, acidic=CL, updating=True)
N = len(water_reactionCoordinate)
print('Number of oxygen atoms on pathway', N)
z_coord = np.zeros(N)
O_index = np.zeros(N)
for i in range(N):
    z_coord[i] = water_reactionCoordinate[i].position[2]
    O_index[i] = water_reactionCoordinate[i].index
ind_sort = np.argsort(z_coord)
com_dist = z_coord[ind_sort] - com_z
#print('sorted z-distances - com_z ', com_dist)
sorted_indices = O_index[ind_sort]

print('distances vec', 10-com_z, 190-com_z)

### visualise reaction coordinate
v = nv.show_mdanalysis(protein)
#v.add_cartoon(component=0)
#v.add_representation('spacefill',  selection='#C', color='blue', radius='1.2')
v.add_representation('cartoon', color='blue')# 

v.add_trajectory(lipids_phosp)
v.add_spacefill(component=1)
v.add_trajectory(Cl_within_pore)
v.add_spacefill(component=2)
#v.add_trajectory(pbc_atom)
#v.add_spacefill(component=3)
v.add_trajectory(water_reactionCoordinate)
v.add_ball_and_stick(component=3, aspectRatio=5)
v.add_trajectory(ion)
v.add_spacefill(component=4)

v._remote_call("setSize", target="Widget", args=["800px", "800px"])
v.center()
#v.render_image()
#v.download_image('glut.jpeg')
v

# replace water by ion
def find_closest_z(d, com_dist):
    for i in range(len(com_dist)):
        if com_dist[i]>d:
            d2 = abs(com_dist[i]-d)
            d1 = abs(com_dist[i-1]-d)
            if d1<d2:
                print('smaller', i-1, com_dist[i-1])
                return i #ind
            else:
                print(i, com_dist[i])
                return i

#folder = '6ply_GlyR_pore/6ply_GlyR_pore/'
#folder = '6pm2_GlyR_pore/6pm2_GlyR_pore/'
#folder = '7m6q_GlyR_pore/7m6q_GlyR_pore/'
#folder = '7m6s_GlyR_pore/7m6s_GlyR_pore/'
#folder = '6pm6_GlyR_pore/'
folder = '6pm2_GlyR_Pore/' # SP
p = '/biggin/b232/temp-pritchard/Documents/PMFs/' + folder

vec = np.arange(-45,46)
ind = 0

ion = u.select_atoms('bynum 28728') # 6pm2 CLA SP
# ion = u.select_atoms('bynum 28687') # 6pm2 POT SP
# ion = u.select_atoms('bynum 28687') # 6pm2 SOD SP

print(CL[0].index)
print(CL[0].position, ion.positions)
print()
com = protein.center_of_mass()

for dist in vec:
    """group1 = u.select_atoms('name *')
    ### find closest z-coordinate in water_reactionCoordinate ###
    ind = find_closest_z(dist, com_dist)
    print('ind', ind)
    index = sorted_indices[ind]
    water_mol = u.select_atoms('index '+str(int(index))+':'+str(int(index+2)))
    print(water_mol.positions,water_mol.types )
    ### Store original positions ###
    ion_pos = ion[0].position.copy()
    water_pos = water_mol[0].position.copy()
    ### Move chloride to water position ###
    ion[0].position = water_mol[0].position
    print(ion.positions)
    ### Compute translation vector for water molecule ###
    translation = ion_pos - water_pos
    water_mol.positions += translation
    ind = ind + 1
     ### Check for clashes ###
    clash_ion = u.select_atoms("(around 1.0 group ion) and not group ion and not group water_mol", ion=ion, water_mol=water_mol)
    clash_water = u.select_atoms("(around 1.0 group water_mol) and not group ion and not group water_mol", ion=ion, water_mol=water_mol)
    new_dist = np.linalg.norm(ion[0].position - COM)
    print('Distance ', dist, 'real dist', new_dist,
          ' - clashes (ion): ', len(clash_ion), ' - clashes (water): ', len(clash_water))

    ### Write out starting configuration - Needs to be altered for other ions###
    if not os.path.exists(p + '/CONFIGS_SOD/'):
        os.makedirs(p + '/CONFIGS_SOD/')
    group1.write(p + '/CONFIGS_SOD/conf' + str(dist) + '.gro')
    print()"""

    group1 = u.select_atoms('name *')
    ### find closest z-coordinate in water_reactionCoordinate ###
    ind = find_closest_z(dist, com_dist)
    print('ind', ind)
    index = sorted_indices[ind]
    water_del = u.select_atoms('index '+str(int(index))+':'+str(int(index+2)))
    print(water_del.positions,water_del.types )
    ### Move chloride to water position ###
    ion[0].position = water_del[0].position
    print(ion.positions)
    ### Delete water molecule ###
    group2 = u.select_atoms('not group id1', id1=water_del)
    print('length of group ', len(group1), len(group2))
    ### Check for clashes ###
    ind = ind + 1
    clash1 = u.select_atoms("(around 1 group basic) and (not group id1)", basic=ion, id1=water_del, updating=True)
    clash2 = u.select_atoms("(around 2 group basic) and (not group id1)", basic=ion, id1=water_del, updating=True)
    clash3 = u.select_atoms("(around 3 group basic) and (not group id1)", basic=ion, id1=water_del, updating=True)
    new_dist = np.linalg.norm(ion.center_of_mass() - com)
    print('distance ', dist, 'real dist', new_dist,
          ' - clashed: ', len(clash1), len(clash2), len(clash3), clash2.indices, clash2.types )
    ### Write out starting configuration ###
    if not os.path.exists(p+'/david_CONFIGS_CLA/'):
        os.makedirs(p+'/david_CONFIGS_CLA/')
    group2.write(p+'/david_CONFIGS_CLA/conf'+str(dist)+'.gro')
    print()