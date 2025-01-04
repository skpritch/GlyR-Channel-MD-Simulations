import numpy as np
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Define cutoff distances for ions (values in Å)
ion_cutoff_distances = {
    'CLA': {
        'water_oxygen': 3.75,
        'water_hydrogen': 2.75,
        'protein': 3.0
    },
    'SOD': {
        'water_oxygen': 3.2,
        'water_hydrogen': 2.3,
        'protein': 2.7
    },
    'POT': {
        'water_oxygen': 3.4,
        'water_hydrogen': 2.5,
        'protein': 2.85
    }
}

def hydration_shell_analysis(path, fname, ion_type, ion_selection, water='TIP3', water_o='OH2',
                             water_H1='H1', water_H2='H2', win_min=-58, win_max=61,
                             win_step=2):

    cutoff_distances = ion_cutoff_distances[ion_type]

    # Prepare list to store data
    data = []

    # Loop over windows
    for win in range(win_min, win_max, win_step):
        try:
            # Construct file paths for trajectory and topology files
            traj_file = os.path.join(path, f"{fname}{win}.xtc")
            top_file = os.path.join(path, f"{fname}{win}.tpr")

            # Check if files exist
            if not os.path.isfile(traj_file) or not os.path.isfile(top_file):
                print(f"Files for window {win} not found. Skipping this window.")
                continue

            # Load the Universe
            u = MDAnalysis.Universe(top_file, traj_file, tpr_resid_from_one=True)
            print(f'Window {win}: Loaded trajectory with {len(u.trajectory)} frames.')

            # Select the ion (assuming ion_selection points to a single ion)
            ion = u.select_atoms(ion_selection)
            print(f'Ion selected: {ion}')
            print(ion.positions)

            # Initialize HydrogenBondAnalysis to guess protein hydrogens as mentor did
            hbonds = HBA(universe=u)

            # Select water hydrogens and oxygens (similar to mentor)
            water_H = u.select_atoms(f"resname {water} and name {water_H1} {water_H2}")
            water_O = u.select_atoms(f"resname {water} and name {water_o}")

            # Select protein hydrogens using guess_hydrogens (like mentor's code)
            protein_H = u.select_atoms(hbonds.guess_hydrogens("protein"))
            print(f'Number of protein hydrogens: {len(protein_H)}')

            # Now select atoms around the ion using the "group id2 and around ... group id1" syntax:
            # This matches the mentor's approach, ensuring the counting is done similarly.
            waterH_near_ion = u.select_atoms(
                f"group id2 and around {cutoff_distances['water_hydrogen']} group id1",
                id1=ion, id2=water_H, updating=True
            )
            waterO_near_ion = u.select_atoms(
                f"group id2 and around {cutoff_distances['water_oxygen']} group id1",
                id1=ion, id2=water_O, updating=True
            )
            proteinH_near_ion = u.select_atoms(
                f"group id2 and around {cutoff_distances['protein']} group id1",
                id1=ion, id2=protein_H, updating=True
            )

            print(f'Initial counts - Water H: {len(waterH_near_ion)}, Water O: {len(waterO_near_ion)}, Protein H: {len(proteinH_near_ion)}')

            num_H = []
            num_O = []
            num_Hprotein = []
            ion_z_coords = []

            for ts in u.trajectory:
                num_H.append(len(waterH_near_ion))
                num_O.append(len(waterO_near_ion))
                num_Hprotein.append(len(proteinH_near_ion))
                ion_z_coords.append(ion.positions[0][2])  # Assuming ion is a single atom

            # Calculate means and standard deviations
            mean_num_H = np.mean(num_H)
            std_num_H = np.std(num_H)
            mean_num_O = np.mean(num_O)
            std_num_O = np.std(num_O)
            mean_num_Hprotein = np.mean(num_Hprotein)
            std_num_Hprotein = np.std(num_Hprotein)
            mean_ion_z = np.mean(ion_z_coords)
            std_ion_z = np.std(ion_z_coords)

            print(f'Window {win}:')
            print(f'  Water H: {mean_num_H:.2f} ± {std_num_H:.2f}')
            print(f'  Water O: {mean_num_O:.2f} ± {std_num_O:.2f}')
            print(f'  Protein H: {mean_num_Hprotein:.2f} ± {std_num_Hprotein:.2f}')
            print(f'  Ion Z-coordinate: {mean_ion_z:.2f} ± {std_ion_z:.2f}')

            # Store data
            data.append({
                'window': win,
                'reaction_coordinate': win / 10,  # Use 'win' as the reaction coordinate
                'mean_num_H': mean_num_H,
                'std_num_H': std_num_H,
                'mean_num_O': mean_num_O,
                'std_num_O': std_num_O,
                'mean_num_Hprotein': mean_num_Hprotein,
                'std_num_Hprotein': std_num_Hprotein,
            })

        except Exception as e:
            print(f'ERROR with window {win}: {e}')

    # Convert data to DataFrame
    df = pd.DataFrame(data)

    # Sort DataFrame by reaction coordinate
    df = df.sort_values(by='reaction_coordinate')

    # Plot the interactions versus reaction coordinate
    plt.figure(figsize=(5, 5))
    plt.plot(df['reaction_coordinate'], df['mean_num_H'], 'b-', label='n(Cl-H) water')  # Blue line
    plt.plot(df['reaction_coordinate'], df['mean_num_O'], 'm-', label='n(Cl-O) water')  # Magenta line
    plt.plot(df['reaction_coordinate'], df['mean_num_Hprotein'], 'r-', label='n(Cl-*) protein')  # Red line

    # Add shading to represent standard deviations
    plt.fill_between(df['reaction_coordinate'], df['mean_num_H'] - df['std_num_H'], df['mean_num_H'] + df['std_num_H'], color='blue', alpha=0.3)
    plt.fill_between(df['reaction_coordinate'], df['mean_num_O'] - df['std_num_O'], df['mean_num_O'] + df['std_num_O'], color='magenta', alpha=0.3)
    plt.fill_between(df['reaction_coordinate'], df['mean_num_Hprotein'] - df['std_num_Hprotein'], df['mean_num_Hprotein'] + df['std_num_Hprotein'], color='red', alpha=0.3)

    # Axis labels and title
    plt.xlabel('Reaction Coordinate (Å)', fontsize=12)
    plt.ylabel('Coordination Number', fontsize=12)
    plt.title('Coordination Number (CHARMM36m)', fontsize=14)

    # Improved legend formatting
    plt.legend(fontsize=10, frameon=True, loc='upper center', framealpha=1, borderpad=1)

    # Tight layout and save the figure
    plt.tight_layout()
    plt.savefig(os.path.join(path, f'Hydration_Shell_{ion_type}_styled.png'), dpi=300)
    plt.show()

# Example usage of the function
if __name__ == "__main__":
    hydration_shell_analysis(
        path='/biggin/b232/temp-pritchard/Documents/PMFs/6pm2_GlyR_Pore/PROD/PROD_CLA',
        fname='PROD_',
        ion_type='CLA',
        ion_selection='bynum 28728',  # Adapt this to the correct ion selection
        water='TIP3',
        water_o='OH2',
        water_H1='H1',
        water_H2='H2',
        win_min=-48,
        win_max=49,
        win_step=1,
    )
