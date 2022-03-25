import numpy as np
import pandas as pd
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.utils import shuffle
import sys
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report

variants = pd.read_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/ensembl-vep/variants_VEP.out", comment='#', on_bad_lines='skip', sep = "\t", header = None)
variants = variants.dropna()

variantsOrig = pd.read_csv("/user/home/uw20204/CScapeVariants/trainingVariantsCoding.txt", comment='#', on_bad_lines='skip', sep = "\t")
variantsOrig.pos = pd.to_numeric(variantsOrig.pos, downcast='integer')

df = pd.DataFrame(columns=["chrom", "driver_stat", "splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant",
"stop_lost","start_lost","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant",
"protein_altering_variant","splice_region_variant","incomplete_terminal_codon_variant","start_retained_variant",
"stop_retained_variant","synonymous_variant","coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant",
"3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","NMD_transcript_variant",
"non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification",
"TF_binding_site_variant","regulatory_region_ablation","regulatory_region_amplification","feature_elongation",
"regulatory_region_variant","feature_truncation","intergenic_variant"])
for i in range(0, len(variants[1].unique())):
    myList = variants[variants[1] == variants[1].unique()[i]][6].unique()
    consReturned = []
    for y in [x.split(',') for x in myList]:
        for z in y:
            consReturned.append(z)
    consReturned = np.unique(consReturned)
    consNotReturned = np.unique(df.columns[range(2, len(df.columns))].difference(consReturned))
    driver_stat = variantsOrig.loc[(variantsOrig.chrom == variants[variants[1] == variants[1].unique()[i]].reset_index(drop = True)[1][0].split(":")[0]) 
    & (variantsOrig.pos == int(variants[variants[1] == variants[1].unique()[i]].reset_index(drop = True)[1][0].split(":")[1]))].reset_index(drop = True).driver_status[0]

    chrom = variantsOrig.loc[(variantsOrig.chrom == variants[variants[1] == variants[1].unique()[i]].reset_index(drop = True)[1][0].split(":")[0]) 
    & (variantsOrig.pos == int(variants[variants[1] == variants[1].unique()[i]].reset_index(drop = True)[1][0].split(":")[1]))].reset_index(drop = True).chrom[0]

    df.loc[i, consReturned] = 1
    df.loc[i, "driver_stat"] = driver_stat
    df.loc[i, consNotReturned] = 0
    df.loc[i, "chrom"] = chrom.split(":")[0]

#there were some mislabelled data in the training set
df.loc[df.splice_donor_region_variant == 1, "splice_donor_variant"] = 1
df.loc[df.splice_polypyrimidine_tract_variant == 1, "splice_region_variant"] = 1
df.loc[df.splice_donor_5th_base_variant == 1, "splice_donor_variant"] = 1
df.drop(['splice_donor_region_variant', 'splice_polypyrimidine_tract_variant', 'splice_donor_5th_base_variant'], inplace=True, axis=1)

df.to_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/VEP_training_coding_csv.txt",  sep = "\t", index = None)

