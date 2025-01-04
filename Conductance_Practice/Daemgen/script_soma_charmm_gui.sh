#!/bin/bash
#SBATCH --nodes=1
## The use of ntasks-per-socket instead of ntasks-per-node is *very* important
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J Daemgen_TIP3P_NBFIX
#SBATCH -p gpu-khalid

module purge
module load all/GROMACS/2021.5-CUDA-11.7.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo ${SLURM_CPUS_PER_TASK}

# Directory for production run
mkdir -p PROD

prev=eq5
fname=PROD_pos_ext_pot
gmx grompp -f MDP/PROD_pos_ext_pot.mdp -o PROD/PROD_pos_ext_pot.tpr -c eq5.gro -r eq5.gro -p topol.top -n index.ndx
cd PROD
gmx mdrun -v -deffnm PROD_pos_ext_pot -ntomp ${SLURM_CPUS_PER_TASK} -update gpu
cd ..
