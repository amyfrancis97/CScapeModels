#!/bin/bash

#SBATCH --job-name=BigWigBedGraphTest
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/user/home/uw20204/mrcieu_data/encode/public/Histone+ChIP-seq/released/2021-10-21/data/training/coding

python LOCO_CV_loop.py
