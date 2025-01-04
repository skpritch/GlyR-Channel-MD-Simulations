# Investigating Anion Selectivity of the GlyR Channel Using Imposed Voltage and Umbrella Sampling

**Author:** Simon Pritchard  
**Affiliation:** Stanford University, Department of Biology  
**Contact:** skpritch@stanford.edu  

## Abstract

This research investigates the anion selectivity mechanisms of the glycine receptor (GlyR) channel using molecular dynamics simulations and umbrella sampling techniques. Building upon previous work, this study examines the conductance properties and ion permeation characteristics of the transmembrane domain under various applied potentials, validating and extending the findings of Dr. Seiferth regarding GlyR channel properties. 

Potential of mean force (PMF) calculations reveal substantially higher free energy barriers for cations compared to chloride ions, providing a quantitative basis for the channel's anion selectivity. Additionally, a significant interaction between chloride ions and the K292 residue at the extracellular end of the truncated transmembrane domain was identified using hydration shell plots, highlighting limitations in using simplified channel models for MD simulations. 

These findings demonstrate both the utility and limitations of computational approaches in understanding ion channel selectivity mechanisms.

## PMFs:
- **generate_starting_configurations_PMF.py**: takes a .gro and .xtc file as inputs and generates .gro files for each z-axis ion starting position
- **write_index.sh**: For a given starting CONFIG directory generates the requisite index files
- **write_mdp.pl**: Generates the necessary .mdp files for each run (need to alter the harmonic positional restraint z-axis coordinate for each coordinate)
- **script_wham.sh**: Runs a WHAM analysis from GROMACS across various time-length windows to prepare files for convergence analysis
- **PMF_analysis.py**: Accepts a PROD and WHAM directory and outputs PMF profile graphs
- **hydration_shell.py**: Accepts a PROD directory and calculates interaction numbers for each time window
- **check_system_orientation.py**: evaluates the results of production runs and checks ion position against key pore residue positions

## Conductance Practice
- **6pm2**: contains the forcefield, mdp, topology, and .gro files for the 6pm2 pore
- **Daemgen**: contains the forcefield, mdp, topology, and .gro files for the Daemgen pore
- **cond_practice.py**: calculates the conductance given a .gro and .xtc
- **current_vs_voltage.py**: given a spread of applied voltages calculates the respective currents and plots them
- **transits_vs_membound.py**: given the outputs of a production run plots the conductance estimates for different assigned channel z-axis boundaries

## Figures
- Contains the figures used in the report
