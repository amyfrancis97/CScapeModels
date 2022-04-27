import pandas as pd

# reading in the feature group in CSV format

dataset = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/Histone+ChIP-seq/released/2021-10-21/data/Histone+ChIP-seqSigVals/CSV.txt", header=None)

# reading in the corresponding variants to get genomic positions/ which chromosome the variants is on

variants = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/Histone+ChIP-seq/released/2021-10-21/data/Histone+ChIP-seqSigVals/" +
"variants_incl.txt", names = ["chrom", "pos", "driverStat", "refAllele", "altAllele"], header=None, sep = "\t")
