#!/bin/bash

#SBATCH --job-name=gridSearch
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/user/home/uw20204/mrcieu_data/ucsc/public/k24.Bismap.MultiTrackMappability/released/2021-11-01/

##wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k24.Bismap.MultiTrackMappability.bw

bigWigToBedGraph /user/home/uw20204/mrcieu_data/ucsc/public/k24.Bismap.MultiTrackMappability/released/2021-11-01/k24.Bismap.MultiTrackMappability.bw /user/home/uw20204/mrcieu_data/ucsc/public/k24.Bismap.MultiTrackMappability/released/2021-11-01/k24.Bismap.MultiTrackMappability.bedGraph
