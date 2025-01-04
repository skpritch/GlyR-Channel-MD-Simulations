import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import apply_PBC
import pandas as pd

# Assigns numeric states to ions based on their z-coordinate relative to protein COM
def assign_state(z, z_min, z_max):
    if z > z_max:
        return 1   # 'outside'
    elif z < z_min:
        return -1  # 'inside'
    else:
        return 0   # 'in_channel'

# Tracks ions moving through the channel and calculates net transits
def track_conductance(u, ions, z_min, z_max, xy_cutoff):
    ion_positions = {}
    protein = u.select_atoms('protein')

    for ion in ions:
        ion_positions[ion.index] = []

    # Iterate over every frame in the trajectory
    for ts in u.trajectory:
        box = ts.dimensions
        protein_com = protein.center_of_mass()
        ion_pos = ions.positions
        delta = ion_pos - protein_com
        delta -= box[:3] * np.round(delta / box[:3])

        positions = delta * 0.1  # Convert to nm
        xy_distances = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)

        # Select ions within the x-y cutoff
        within_cutoff = xy_distances <= xy_cutoff
        for i, ion in enumerate(ions):
            ion_index = ion.index
            if within_cutoff[i]:
                z = positions[i, 2]
                state = assign_state(z, z_min, z_max)
                ion_positions[ion_index].append(state)
            else:
                pass  # Ion is outside the x-y cutoff; ignore for transit counting

    # Now, for each ion, process the sequence of states to detect transits with direction
    net_transits = 0
    for ion_index, pos_list in ion_positions.items():
        if len(pos_list) < 3:
            continue  # Not enough data to detect events

        # Create a list of distinct state changes
        events1 = [pos_list[0]]
        for state in pos_list[1:]:
            if state != events1[-1]:
                events1.append(state)

        # Handling PBC errors and avoiding duplicates
        events2 = []
        for count in range(2, len(events1)):
            test = [events1[count-2], events1[count-1], events1[count]]
            if len(set(test)) == 3 and test[1] == 0:
                events2.append(test)

        # Now, count the events with directionality
        for e in events2:
            if e == [1, 0, -1]:
                # Transit from outside to inside
                net_transits += 1
            elif e == [-1, 0, 1]:
                # Transit from inside to outside
                net_transits -= 1
            else:
                # Unexpected transit sequence (possible PBC error)
                print('Unexpected transit sequence for ion', ion_index, ':', e)

    print(f"Net transits: {net_transits}")
    return net_transits, u.trajectory[-1].time

# Calculates conductance based on total charge transported
def calculate_conductance(total_charge, total_time_ps, voltage):
    total_time_sec = total_time_ps * 1e-12  # Convert ps to seconds
    current = total_charge / total_time_sec  # Current in Amperes
    conductance = (current / voltage) * 1e12  # Convert S to pS
    print(f"Conductance: {conductance} pS")
    return conductance

# Processes the simulation for given ion types and calculates conductance
def process_simulation(top, traj, E, ion_types, ion_charges, z_min, z_max, xy_cutoff):
    u = mda.Universe(top, traj)

    # Collect Lz over every 10th frame to compute average voltage
    Lz_list = []
    for ts in u.trajectory[::1]:
        Lz_list.append(ts.dimensions[2])
    Lz_array = np.array(Lz_list)
    avg_Lz = np.mean(Lz_array) * 0.1
    voltage = E * avg_Lz

    u.trajectory.rewind()

    total_charge = 0.0
    total_time_ps = u.trajectory[-1].time

    for ion_type in ion_types:
        ions = u.select_atoms(f'name {ion_type}')
        ion_charge = ion_charges[ion_type]

        transits, time = track_conductance(u, ions, z_min, z_max, xy_cutoff)
        total_charge += transits * ion_charge

    conductance = calculate_conductance(total_charge, total_time_ps, voltage)
    return conductance

# Simulation data
sims = [
    {'pore_label': '6pm2','Cond': 'No NBFIX TIP3P','top': '6pm2_TIP3P/eq5.gro','traj': '6pm2_TIP3P/PROD/PROD_pos_ext_pot.xtc','E': 0.04372456323533785,
        'ion_types': ['CLA', 'POT'],
        'ion_charges': {'CLA': 1.6022e-19,'POT': -1.6022e-19}},
    {'pore_label': '6pm2', 'Cond': 'NBFIX TIP3P', 'top': '6pm2_TIP3P_NBFIX/eq5.gro', 'traj': '6pm2_TIP3P_NBFIX/orig_PROD/PROD_pos_ext_pot.xtc', 'E': 0.04372456323533785,
        'ion_types': ['CLA', 'POT'],
        'ion_charges': {'CLA': 1.6022e-19,'POT': -1.6022e-19}},
    {'pore_label': 'Daemgen', 'Cond': 'No NBFIX TIP3P', 'top': 'Daemgen_TIP3P/eq5.gro', 'traj': 'Daemgen_TIP3P/PROD/PROD_pos_ext_pot.xtc', 'E': 0.04372456323533785,
        'ion_types': ['CLA', 'POT'],
        'ion_charges': {'CLA': 1.6022e-19,'POT': -1.6022e-19}},
    {'pore_label': 'Daemgen', 'Cond': 'NBFIX TIP3P', 'top': 'Daemgen_TIP3P_NBFIX/eq5.gro', 'traj': 'Daemgen_TIP3P_NBFIX/PROD/PROD_pos_ext_pot.xtc', 'E': 0.04372456323533785,
        'ion_types': ['CLA', 'POT'],
        'ion_charges': {'CLA': 1.6022e-19,'POT': -1.6022e-19}}
]

# z_min and z_max relative to protein COM (in nm)
z_min = -1.8  # Lower boundary in nm
z_max = 1.9   # Upper boundary in nm

# x-y cutoff radius in nm
xy_cutoff = 0.7

# Initialize lists to store results
conductances = []
pore_labels = []
conditions = []

for sim in sims:
    cond = process_simulation(
        sim['top'],
        sim['traj'],
        sim['E'],
        sim['ion_types'],
        sim['ion_charges'],
        z_min,
        z_max,
        xy_cutoff)
    conductances.append(cond)
    pore_labels.append(sim['pore_label'])
    conditions.append(sim['Cond'])

# Create a DataFrame with the results
results_df = pd.DataFrame({
    'Pore Label': pore_labels,
    'Condition': conditions,
    'Conductance (pS)': conductances
})

# Pivot the DataFrame to get the desired 2x2 table format
table = results_df.pivot(index='Pore Label', columns='Condition', values='Conductance (pS)')
print(table)