#!/bin/bash

#SBATCH --job-name=LOCO-CV-cons
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/user/home/uw20204/scratch/CScapeModels/conservation

python loco_final.py
