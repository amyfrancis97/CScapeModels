import pandas as pd 

# reading in the feature group in CSV format
dataset = pd.read_csv("/user/home/uw20204/mrcieu_data/ucsc/public/hg38.phastCons30way/released/2021-10-27/" + "hg38.phastCons30way.trainingVariantsCoding_web.CSV")
dataset2 = pd.read_csv("/user/home/uw20204/mrcieu_data/ucsc/public/hg38.phyloP30way/released/2021-10-27/" + "hg38.phyloP30way.trainingVariantsCoding_web.CSV")
featureGroup = "30WayCons"
