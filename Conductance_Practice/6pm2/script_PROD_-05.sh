#!/bin/bash
#SBATCH --nodes=1
## The use of ntasks-per-socket instead of ntasks-per-node is *very* important
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J neg05_6pm2_TIP3P_NBFIX
#SBATCH -p gpu-khalid

module purge
module load all/GROMACS/2021.5-CUDA-11.7.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo ${SLURM_CPUS_PER_TASK}

# Directory for production run
mkdir -p PROD_-05
prev=eq5
fname=PROD_-05_ext_pot_NVT
gmx grompp -f MDP/PROD_-05_ext_pot_NVT.mdp -o PROD_-05/PROD_-05_ext_pot_NVT.tpr -c eq5.gro -r eq5.gro -p topol.top -n index.ndx -maxwarn 1
cd PROD_-05
gmx mdrun -v -deffnm PROD_-05_ext_pot_NVT -ntomp ${SLURM_CPUS_PER_TASK} -update gpu -cpi
cd ..


