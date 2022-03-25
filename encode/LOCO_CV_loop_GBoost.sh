#!/bin/bash

#SBATCH --job-name=LOCO-CV_GBOOST
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/user/home/uw20204/scratch/CScapeModels

cd /user/home/uw20204/scratch/CScapeModels/encode
python LOCO_CV_loop_GBoost.py
