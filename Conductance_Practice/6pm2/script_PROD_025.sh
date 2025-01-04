#!/bin/bash
#SBATCH --nodes=1
## The use of ntasks-per-socket instead of ntasks-per-node is *very* important
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J 025_6pm2_TIP3P_NBFIX
#SBATCH -p gpu-khalid

module purge
module load all/GROMACS/2021.5-CUDA-11.7.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo ${SLURM_CPUS_PER_TASK}

# Directory for production run
mkdir -p PROD_025
prev=eq5
fname=PROD_025_ext_pot_NVT
gmx grompp -f MDP/PROD_025_ext_pot_NVT.mdp -o PROD_025/PROD_025_ext_pot_NVT.tpr -c eq5.gro -r eq5.gro -p topol.top -n index.ndx -maxwarn 1
cd PROD_025
gmx mdrun -v -deffnm PROD_025_ext_pot_NVT -ntomp ${SLURM_CPUS_PER_TASK} -update gpu -cpi
cd ..


