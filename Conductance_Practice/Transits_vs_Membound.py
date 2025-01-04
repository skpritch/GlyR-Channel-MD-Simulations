import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# Assigns states to ions based on their z-coordinate
def assign_state(z, z_min, z_max):
    if z > z_max:
        return 'outside'
    elif z < z_min:
        return 'inside'
    else:
        return 'in_channel'

# Tracks ions moving through the channel and calculates net transits
def track_conductance(u, ions, z_min, z_max):
    transits_inside = 0
    transits_outside = 0
    ion_states = {}

    # Initialize ion states based on initial positions
    for ion in ions:
        z = ion.position[2] * 0.1  # Convert Å to nm
        ion_states[ion.index] = assign_state(z, z_min, z_max)

    # Iterate over trajectory
    for ts in u.trajectory:
        positions = ions.positions * 0.1  # Convert Å to nm
        for i, ion in enumerate(ions):
            z = positions[i, 2]
            ion_index = ion.index

            current_state = assign_state(z, z_min, z_max)
            previous_state = ion_states.get(ion_index, 'unknown')

            # Update ion state and detect transits
            if current_state == 'outside':
                if previous_state == 'towards_outside':
                    transits_outside += 1  # Ion completed transit towards positive z
                ion_states[ion_index] = 'outside'
            elif current_state == 'inside':
                if previous_state == 'towards_inside':
                    transits_inside += 1  # Ion completed transit towards negative z
                ion_states[ion_index] = 'inside'
            elif current_state == 'in_channel':
                if previous_state == 'outside':
                    ion_states[ion_index] = 'towards_inside'
                elif previous_state == 'inside':
                    ion_states[ion_index] = 'towards_outside'
                else:
                    ion_states[ion_index] = previous_state
            else:
                ion_states[ion_index] = 'unknown'

    transits_net = transits_inside - transits_outside
    print(f"Transits Inside: {transits_inside}, Transits Outside: {transits_outside}")
    return transits_net, u.trajectory[-1].time

# Calculates conductance based on total charge transported
def calculate_conductance(total_charge, total_time_ps, voltage):
    total_time_sec = total_time_ps * 1e-12  # Convert ps to s
    current = total_charge / total_time_sec  # Current in Amperes
    conductance = (current / voltage) * 1e12  # Convert S to pS
    print(f"Conductance: {conductance} pS")
    return conductance

# Processes the simulation for given ion types and calculates conductance
def process_simulation(top, traj, E, ion_types, ion_charges, z_min, z_max):
    u = mda.Universe(top, traj)
    voltage = E * u.dimensions[2] * 0.1  # Convert Å to nm

    total_charge = 0.0
    total_time_ps = u.trajectory[-1].time
    
    ion_transits = {}

    for ion_type in ion_types:
        ions = u.select_atoms(f'name {ion_type}')
        ion_charge = ion_charges[ion_type]

        transits_net, time = track_conductance(u, ions, z_min, z_max)
        total_charge += transits_net * ion_charge
        ion_transits[ion_type] = transits_net

    conductance = calculate_conductance(total_charge, total_time_ps, voltage)
    return conductance, ion_transits

# Simulation data
sim = {
    'pore_label': '6pm2', 
    'Cond': 'NBFIX TIP3P', 
    'top': '6pm2_TIP3P_NBFIX/eq5.gro', 
    'traj': '6pm2_TIP3P_NBFIX/PROD_05/PROD_05_ext_pot_NVT.xtc', 
    'E': 0.0457611,
    'ion_types': ['CLA', 'POT'],
    'ion_charges': {'CLA': 1.6022e-19, 'POT': -1.6022e-19}
}

# Original z_min and z_max
z_min = 4.677
z_max = 6.439

# Calculate the center point of the channel
center = (z_min + z_max) / 2

# Define half-widths
half_widths = np.arange(0.1, 5.1, 0.2)

# Initialize lists to store results
conductances = []
cla_transits_list = []
pot_transits_list = []
widths = []

for half_width in half_widths:
    adjusted_z_min = center - half_width
    adjusted_z_max = center + half_width

    print(f"\nCalculating for half-width: {half_width} nm")
    conductance, ion_transits = process_simulation(
        sim['top'],
        sim['traj'],
        sim['E'],
        sim['ion_types'],
        sim['ion_charges'],
        adjusted_z_min,
        adjusted_z_max
    )
    conductances.append(conductance)
    cla_transits_list.append(ion_transits['CLA'])
    pot_transits_list.append(ion_transits['POT'])
    widths.append(half_width)

# ---- Styling Setup ----
sns.set_theme(style='whitegrid')

mpl.rc('font', size=12)           # Default text font size
mpl.rc('axes', titlesize=14)      # Axes title font size
mpl.rc('axes', labelsize=14)      # Axes label font size
mpl.rc('xtick', labelsize=12)     # X tick label size
mpl.rc('ytick', labelsize=12)     # Y tick label size
mpl.rc('legend', fontsize=12)     # Legend font size
mpl.rc('lines', linewidth=2)      # Line thickness

# Plot net transits for cations and anions vs. half-width
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(widths, cla_transits_list, marker='o', color='blue', label='Cl⁻ Net Transits')
ax.plot(widths, pot_transits_list, marker='o', color='red', label='K⁺ Net Transits')
ax.set_xlabel('Half-width of Channel (nm)')
ax.set_ylabel('Net Transits')
ax.set_title('Ion Net Transits vs. Half-width of Channel', fontsize=14, pad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.7)
ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

# Place legend in the center of the plot
ax.legend(loc='center', frameon=False)

plt.tight_layout()
plt.show()

# Plot conductance vs. half-width
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(widths, conductances, marker='o', color='green', label='Conductance')
ax.set_xlabel('Half-width of Channel (nm)')
ax.set_ylabel('Conductance (pS)')
ax.set_title('Conductance vs. Half-width of Channel', fontsize=14, pad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.7)
ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

# Place legend in the center of the plot
ax.legend(loc='center', frameon=False)

plt.tight_layout()
plt.show()