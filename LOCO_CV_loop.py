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

# reading in the feature group in CSV format

dataset = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/TF+ChIP-seq/released/2021-10-21/data/training/coding/" +
"CSV.txt", header=None)

# reading in the corresponding variants to get genomic positions/ which chromosome the variants is on

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

    # select stated chrom to hold out as the test set
    datasetTest = dataset[dataset["chrom"] == heldOutChrom]

    X_test = datasetTest.drop(["class", "altAllele", "refAllele", "chrom", "pos", "driverStat"], axis=1) # test dataset X are the signal values ONLY
    y_test = datasetTest["class"] # test dataset y, or labels, are the driver status ONLY

    # hold out stated chrom to generate training dataset
    datasetTrain = dataset[dataset["chrom"] != heldOutChrom]

    print("the number of samples in the training set are:" + str(len(datasetTrain)))

    # randomly sample 3000 +ive and 3000 -ive examples from the training dataset to carry forward
    datasetTrain = pd.concat([datasetTrain[datasetTrain["class"] == 1].sample(3000), datasetTrain[datasetTrain["class"] == -1].sample(3000)])
    datasetTrain = shuffle(datasetTrain)
    datasetTrain = datasetTrain.reset_index(drop = True)

    X_train = datasetTrain.drop(["class", "altAllele", "refAllele", "chrom", "pos", "driverStat"], axis=1) # train dataset X are the signal values ONLY
    y_train = datasetTrain["class"] # train dataset y, or labels, are the driver status ONLY

    # specify the model parameters for training, based on previous grid search
    svclassifier = SVC(kernel='rbf', C = 10, gamma = 0.0001)
    svclassifier.fit(X_train, y_train) # fit the model
    y_pred = svclassifier.predict(X_test) # validate the model using all of the available data for the left out chromosome
    df = (classification_report(y_test, y_pred, output_dict=True)) # generate a classification report
    df = df.get('weighted avg').get('precision') # get the weighted average from the classification report
    dfWeightedAvPrec.append(df)

# get the average of the weighted averages, for each chromosome in LOCO-CV
def Average(lst):
    return sum(lst) / len(lst)
print(round(Average(dfWeightedAvPrec), 3))
