#!/bin/bash

#SBATCH --job-name=GMD-Loop

#SBATCH --partition=norm

#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=09:00:00
#SBATCH -o example1.out

source /home/$USER/.bashrc
hostname -s

mamba activate FNLGMD

WRK="/mnt/projects/RAS-CompChem/static/Mayukh/FNLGMD"
echo $WRK
cd $WRK

CODE="/mnt/projects/RAS-CompChem/static/Mayukh/FNLGMD/source/main.py"
CONF="/mnt/projects/RAS-CompChem/static/Mayukh/FNLGMD/examples/LogP_JTVAE/config.yaml"
echo "Before"
python $CODE -config $CONF
echo "After"
