#!/bin/bash

featureGroups=("Histone+ChIP-seq" "TF+ChIP-seq" "eCLIP" "ChIA-PET")
for featureGroup in ${featureGroups[@]}; do sbatch gridSearch.sh ${featureGroup}; done
