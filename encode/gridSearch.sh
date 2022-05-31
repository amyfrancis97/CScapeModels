#!/bin/bash

#SBATCH --job-name=gridSearch
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/user/home/uw20204/scratch/CScapeModels/encode

feature="$1"
python gridSearchFinal.py "${feature}"
