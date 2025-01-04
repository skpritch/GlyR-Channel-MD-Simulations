import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from scipy.optimize import curve_fit
import seaborn as sns

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
    total_transits = 0
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
            total_transits += 1  # Count this transit
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
    print(f"Total transits: {total_transits}")
    return net_transits, total_transits, u.trajectory[-1].time

# Calculates conductance and current based on total charge transported
def calculate_conductance(total_charge, total_time_ps, voltage):
    total_time_sec = total_time_ps * 1e-12  # Convert ps to seconds
    current = total_charge / total_time_sec  # Current in Amperes
    conductance = (current / voltage) * 1e12  # Convert S to pS
    print(f"Conductance: {conductance} pS")
    return conductance, current

# Processes the simulation for given ion types and calculates conductance and current
def process_simulation(top, traj, E, ion_types, ion_charges, z_min, z_max, xy_cutoff):
    u = mda.Universe(top, traj)

    # Collect Lz over every frame to compute average voltage
    Lz_list = []
    for ts in u.trajectory:
        Lz_list.append(ts.dimensions[2])
    Lz_array = np.array(Lz_list)
    avg_Lz = np.mean(Lz_array) * 0.1  # Convert Å to nm
    voltage = E * avg_Lz

    u.trajectory.rewind()

    total_charge = 0.0
    total_time_ps = u.trajectory[-1].time

    total_transits_all = 0
    total_transits_per_ion = {}
    total_transits_CLA = 0

    for ion_type in ion_types:
        ions = u.select_atoms(f'name {ion_type}')
        ion_charge = ion_charges[ion_type]

        transits, total_transits, time = track_conductance(u, ions, z_min, z_max, xy_cutoff)
        total_charge += transits * ion_charge
        total_transits_all += total_transits
        total_transits_per_ion[ion_type] = total_transits
        if ion_type == 'CLA':
            total_transits_CLA += total_transits

    total_transits_others = total_transits_all - total_transits_CLA

    if total_transits_others > 0:
        CLA_ratio = total_transits_CLA / total_transits_others
    else:
        CLA_ratio = np.nan

    conductance, current = calculate_conductance(total_charge, total_time_ps, voltage)

    if total_transits_all > 0:
        error = abs(current) / np.sqrt(total_transits_all)
    else:
        error = np.nan

    return conductance, current, voltage, error, CLA_ratio

# Simulation data
sims = [
    {'pore_label': '6pm2', 'Cond': 'NBFIX TIP3P', 'top': '6pm2_TIP3P_NBFIX/eq5.gro', 'traj': '6pm2_TIP3P_NBFIX/PROD_05/PROD_05_ext_pot_NVT.xtc', 'E': 0.0457611,
     'ion_types': ['CLA', 'POT'],
     'ion_charges': {'CLA': 1.6022e-19, 'POT': -1.6022e-19}},
     {'pore_label': '6pm2', 'Cond': 'NBFIX TIP3P', 'top': '6pm2_TIP3P_NBFIX/eq5.gro', 'traj': '6pm2_TIP3P_NBFIX/PROD_025/PROD_025_ext_pot_NVT.xtc', 'E': 0.022881,
     'ion_types': ['CLA', 'POT'],
     'ion_charges': {'CLA': 1.6022e-19, 'POT': -1.6022e-19}},
    {'pore_label': '6pm2', 'Cond': 'NBFIX TIP3P', 'top': '6pm2_TIP3P_NBFIX/eq5.gro', 'traj': '6pm2_TIP3P_NBFIX/PROD_-025/PROD_-025_ext_pot_NVT.xtc', 'E': -0.022881,
     'ion_types': ['CLA', 'POT'],
     'ion_charges': {'CLA': 1.6022e-19, 'POT': -1.6022e-19}},
    {'pore_label': '6pm2', 'Cond': 'NBFIX TIP3P', 'top': '6pm2_TIP3P_NBFIX/eq5.gro', 'traj': '6pm2_TIP3P_NBFIX/PROD_-05/PROD_-05_ext_pot_NVT.xtc', 'E': -0.0457611,
     'ion_types': ['CLA', 'POT'],
     'ion_charges': {'CLA': 1.6022e-19, 'POT': -1.6022e-19}},
]

# z_min and z_max relative to protein COM (in nm)
z_min = -1.8  # Lower boundary in nm
z_max = 1.9   # Upper boundary in nm

# x-y cutoff radius in nm
xy_cutoff = 0.7

# Initialize lists to store results
conductances = []
currents = []
voltages = []
pore_labels = []
conditions = []
errors = []
CLA_ratios = []

for sim in sims:
    conductance, current, voltage, error, CLA_ratio = process_simulation(
        sim['top'],
        sim['traj'],
        sim['E'],
        sim['ion_types'],
        sim['ion_charges'],
        z_min,
        z_max,
        xy_cutoff)
    conductances.append(conductance)
    currents.append(current)
    voltages.append(voltage)
    pore_labels.append(sim['pore_label'])
    conditions.append(sim['Cond'])
    errors.append(error)
    CLA_ratios.append(CLA_ratio)

# Set pandas display options to show all columns without truncation
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)

# Create a DataFrame with the results
results_df = pd.DataFrame({
    'Pore Label': pore_labels,
    'Condition': conditions,
    'Voltage (V)': voltages,
    'Current (A)': currents,
    'Error (A)': errors,
    'Conductance (pS)': conductances,
    'CLA Ratio': CLA_ratios
})

print(results_df)

# Use a style that is clean and simple
sns.set_theme(style='whitegrid')

# Adjust rcParams for a publication-quality look
mpl.rc('font', size=12)            # Default text font size
mpl.rc('axes', titlesize=14)       # Axes title font size
mpl.rc('axes', labelsize=14)       # Axes label font size
mpl.rc('xtick', labelsize=12)      # X tick label size
mpl.rc('ytick', labelsize=12)      # Y tick label size
mpl.rc('legend', fontsize=12)      # Legend font size
mpl.rc('lines', linewidth=2)       # Line thickness

# Increase figure size
fig, ax = plt.subplots(figsize=(6, 4))

# Plot current vs voltage with error bars
ax.errorbar(voltages, currents, yerr=errors, fmt='o', 
            color='black', capsize=3, elinewidth=1, markeredgewidth=1,
            markersize=6, label='Data')

# Add linear fit through the origin if multiple data points are available
if len(voltages) > 1:
    def linear_func(x, m):
        return m * x
    popt, _ = curve_fit(linear_func, voltages, currents)
    m = popt[0]
    ax.plot(voltages, linear_func(np.array(voltages), m), "g-", linewidth=1.5,
            label=f'Fit Through Origin: y={m:.2e}x')

# Remove top and right spines for a cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set axis labels with units
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Current (A)')

# Use a concise, descriptive title or remove if the figure caption will provide context
ax.set_title('6pm2 Pore: The Effect of Voltage on Current', fontsize=14, pad=10)

# Lighten the grid lines
ax.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.7)
ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

ax.legend(loc='upper left', frameon=False)

# Tight layout for better spacing
plt.tight_layout()

plt.show()