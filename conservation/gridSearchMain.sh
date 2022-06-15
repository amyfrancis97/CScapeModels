#!/bin/bash

dataset1Array=('hg38.phastCons4waytrainingVariantsCoding_web.CSV' 'hg38.phastCons17waytrainingVariantsCoding_web.CSV' 'hg38.phastCons20waytrainingVariantsCoding_web.CSV' 'hg38.phastCons30way.trainingVariantsCoding_web.CSV' 'k24.Bismap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV' 'k36.Bismap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV' 'k50.Bismap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV' 'k100.Bismap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV')
dataset2Array=('hg38.phyloP4waytrainingVariantsCoding_web.CSV' 'hg38.phastCons17waytrainingVariantsCoding_web.CSV' 'hg38.phyloP20waytrainingVariantsCoding_web.CSV' 'hg38.phyloP30way.trainingVariantsCoding_web.CSV' 'k24.Umap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV' 'k36.Umap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV' 'k50.Umap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV' 'k100.Umap.MultiTrackMappabilitytrainingVariantsCoding_web.CSV')
featureGroupArray=('4WayCons' '17WayCons' '20WayCons' '30WayCons' '24MultiTrackMappability' '36MultiTrackMappability' '50MultiTrackMappability' '100MultiTrackMappability')

for i in {0..8}; do sbatch gridSearch.sh "${dataset1Array[$i]}" "${dataset2Array[$i]}" "${featureGroupArray[$i]}"; echo ${featureGroupArray[$i]}; done

