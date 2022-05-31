#!/bin/bash

#SBATCH --job-name=gridSearch
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/user/home/uw20204/scratch/CScapeModels/conservation

dataset1=("$1")
dataset2=("$2")
featureGroup=("$3")

python gridSearchFinal.py ${dataset1} ${dataset2} ${featureGroup}
