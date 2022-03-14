import numpy as np
import pandas as pd
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

dataset = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/TF+ChIP-seq/released/2021-10-21/data/training/coding/" +
"CSV.txt", header=None)

variants = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/TF+ChIP-seq/released/2021-10-21/data/training/coding/" +
"variants_incl.txt", names = ["chrom", "pos", "driverStat", "refAllele", "altAllele"], header=None, sep = "\t")

variants = variants.iloc[range(0, len(dataset)), :]

dataset_variants = pd.concat([variants, dataset], axis=1)
rows_with_nan = [index for index, row in dataset_variants.iterrows() if row.isnull().any()]
dataset_variants = dataset_variants.drop(labels=rows_with_nan, axis=0)

dataset = dataset_variants.rename(columns={0: "class"})

dfWeightedAvPrec = []
for chrom in range(1, 22):
    heldOutChrom = "chr" + str(chrom)
    # LOCO-CV for loop
    datasetTest = dataset[dataset["chrom"] == heldOutChrom]
    X_test = datasetTest.drop(["class", "altAllele", "refAllele", "chrom", "pos", "driverStat"], axis=1)
    y_test = datasetTest["class"]

    datasetTrain = dataset[dataset["chrom"] != heldOutChrom]
    print("the number of samples in the training set are:" + str(len(datasetTrain)))
    datasetTrain = pd.concat([datasetTrain[datasetTrain["class"] == 1].sample(3000), datasetTrain[datasetTrain["class"] == -1].sample(3000)])
    datasetTrain = shuffle(datasetTrain)
    datasetTrain = datasetTrain.reset_index(drop = True)
    X_train = datasetTrain.drop(["class", "altAllele", "refAllele", "chrom", "pos", "driverStat"], axis=1)
    y_train = datasetTrain["class"]

    svclassifier = SVC(kernel='rbf', C = 10, gamma = 0.0001)
    svclassifier.fit(X_train, y_train)
    y_pred = svclassifier.predict(X_test)
    df = (classification_report(y_test, y_pred, output_dict=True))
    df = df.get('weighted avg').get('precision')
    dfWeightedAvPrec.append(df)
def Average(lst):
    return sum(lst) / len(lst)
print(round(Average(dfWeightedAvPrec), 3))
