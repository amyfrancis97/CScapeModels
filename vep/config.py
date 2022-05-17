import pandas as pd 


consequences = pd.read_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/VEP_training_coding_csv_transcriptNoIncl_web.txt", sep = "\t")
aminoAcids = pd.read_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/VEP_web_aminoAcid.txt", sep = "\t")
