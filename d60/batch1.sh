#! /bin/bash -l
####! /bin/tcsh
#PBS -l walltime=07:30:00
#PBS -l nodes=1:ppn=12
#SBATCH --time=24:00:00
#SBATCH --mem=256gb
#     #SBATCH --ntasks-per-node=1
df -h 
echo $PATH
module list 

dir /home/mat236/miniconda/bin

cd /home/mat236/dcfp/matear-nb/d60
conda activate pangeo

python bgc2.py 
