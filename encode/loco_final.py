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
import ast
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import GradientBoostingClassifier
from config import *

# SVM classifier
def LOCO_SVM(dataset, kernel, C, gamma):
    rows_with_nan = [index for index, row in dataset.iterrows() if row.isnull().any()]
    dataset = dataset.drop(labels=rows_with_nan, axis=0)

    dfWeightedAvPrec2 = []
    for i in range(1, 30):
        dfWeightedAvPrec1 = []
        for chrom in range(1, 22):
            
            heldOutChrom = "chr" + str(chrom)

            # LOCO-CV for loop

            # select stated chrom to hold out as the test set
            datasetTest = dataset[dataset["chrom"] == heldOutChrom]

            X_test = datasetTest.drop(["class", "chrom"], axis=1) # test dataset X are the signal values ONLY
            y_test = datasetTest["class"] # test dataset y, or labels, are the driver status ONLY

            # hold out stated chrom to generate training dataset
            datasetTrain = dataset[dataset["chrom"] != heldOutChrom]

            # randomly sample 3000 +ive and 3000 -ive examples from the training dataset to carry forward
            datasetTrain = pd.concat([datasetTrain[datasetTrain["class"] == 1].sample(2000), datasetTrain[datasetTrain["class"] == -1].sample(2000)])
            datasetTrain = shuffle(datasetTrain)
            datasetTrain = datasetTrain.reset_index(drop = True)

            X_train = datasetTrain.drop(["class", "chrom"], axis=1) # train dataset X are the signal values ONLY
            y_train = datasetTrain["class"] # train dataset y, or labels, are the driver status ONLY

            # specify the model parameters for training, based on previous grid search
            svclassifier = SVC(kernel=kernel, C = C, gamma = gamma)
            svclassifier.fit(X_train, y_train) # fit the model
            y_pred = svclassifier.predict(X_test) # validate the model using all of the available data for the left out chromosome
            df = (classification_report(y_test, y_pred, output_dict=True)) # generate a classification report
            df = df.get('weighted avg').get('precision') # get the weighted average from the classification report
            dfWeightedAvPrec1.append(df)

        # get the average of the weighted averages, for each chromosome in LOCO-CV
        def Average(lst):
            return sum(lst) / len(lst)
        dfWeightedAvPrec2.append(Average(dfWeightedAvPrec1))

    # get the average of the weighted averages, for each round of LOCO (30 rounds)
    def Average(lst):
        return sum(lst) / len(lst)
    weightedAvFinal = round(Average(dfWeightedAvPrec2), 4)
    return weightedAvFinal

# GB classifier
# SVM classifier
def LOCO_GB(dataset, learning_rate, n_estimators, max_depth):
    rows_with_nan = [index for index, row in dataset.iterrows() if row.isnull().any()]
    dataset = dataset.drop(labels=rows_with_nan, axis=0)

    dfWeightedAvPrec2 = []
    for i in range(1, 30):
        dfWeightedAvPrec1 = []
        for chrom in range(1, 22):
            
            heldOutChrom = "chr" + str(chrom)

            # LOCO-CV for loop

            # select stated chrom to hold out as the test set
            datasetTest = dataset[dataset["chrom"] == heldOutChrom]

            X_test = datasetTest.drop(["class", "chrom"], axis=1) # test dataset X are the signal values ONLY
            y_test = datasetTest["class"] # test dataset y, or labels, are the driver status ONLY

            # hold out stated chrom to generate training dataset
            datasetTrain = dataset[dataset["chrom"] != heldOutChrom]

            # randomly sample 3000 +ive and 3000 -ive examples from the training dataset to carry forward
            datasetTrain = pd.concat([datasetTrain[datasetTrain["class"] == 1].sample(2000), datasetTrain[datasetTrain["class"] == -1].sample(2000)])
            datasetTrain = shuffle(datasetTrain)
            datasetTrain = datasetTrain.reset_index(drop = True)

            X_train = datasetTrain.drop(["class", "chrom"], axis=1) # train dataset X are the signal values ONLY
            y_train = datasetTrain["class"] # train dataset y, or labels, are the driver status ONLY

        # specify the model parameters for training, based on previous grid search
        gb_clf = GradientBoostingClassifier(learning_rate=learning_rate, n_estimators=n_estimators,max_depth=max_depth, min_samples_split=2, min_samples_leaf=1, subsample=1)
        gb_clf.fit(X_train, y_train) # fit the model
        y_pred = gb_clf.predict(X_test) # validate the model using all of the available data for the left out chromosome
        df = (classification_report(y_test, y_pred, output_dict=True)) # generate a classification report
        df = df.get('weighted avg').get('precision') # get the weighted average from the classification report
        dfWeightedAvPrec1.append(df)

        # get the average of the weighted averages, for each chromosome in LOCO-CV
        def Average(lst):
            return sum(lst) / len(lst)
        dfWeightedAvPrec2.append(Average(dfWeightedAvPrec1))

    # get the average of the weighted averages, for each round of LOCO (30 rounds)
    def Average(lst):
        return sum(lst) / len(lst)
    weightedAvFinal = round(Average(dfWeightedAvPrec2), 4)
    return weightedAvFinal

# reading in best params file
df = pd.read_csv("bestparamsnew2.txt")

# carry out SVM models
for featureGroup in df["featureGroup"].unique():
    # reading in the feature group in CSV format
    dataset = pd.read_csv(location + featureGroup + "_csv.txt", header=None)

    # reading in the corresponding variants to get genomic positions/ which chromosome the variants is on
    variants = pd.read_csv(location + featureGroup + "_variants.txt", names = ["chrom", "pos", "driverStat", "refAllele", "altAllele"], header=None, sep = "\t")

    # reformat to match criteria for grid search function
    variants = variants.iloc[range(0, len(dataset)), :]
    dataset_variants = pd.concat([variants, dataset], axis=1)
    dataset_variants = dataset_variants.drop(columns = ["pos", "refAllele", "altAllele", 0])
    dataset = dataset_variants.rename(columns={"driverStat": "class"})

    # carrying out function for LOCO-CV SVM and writing into parameter file.
    # pull out list of best params for each model
    myDict = ast.literal_eval(df.loc[(df["featureGroup"] == featureGroup) & (df["model"] == "SVM"), "bestParams"].reset_index(drop = True)[0])
    myDict = ast.literal_eval(df.loc[(df["featureGroup"] == featureGroup) & (df["model"] == "SVM"), "bestParams"].reset_index(drop = True)[0])
    if len(myDict.keys()) == 3:
        consRes = LOCO_SVM(dataset = dataset, kernel = list(myDict.values())[2], C = list(myDict.values())[0], gamma = list(myDict.values())[1])
        df.loc[(df["featureGroup"] == featureGroup) & (df["model"] == "SVM"), "LOCO_weightedAv"] = consRes
    else:
        print(list(myDict.values())[0], list(myDict.values())[1])
        consRes = LOCO_SVM(dataset = dataset, kernel = list(myDict.values())[1], C = list(myDict.values())[0], gamma = "auto")
        df.loc[(df["featureGroup"] == featureGroup) & (df["model"] == "SVM"), "LOCO_weightedAv"] = consRes

    # carry out GB models
    myDict = ast.literal_eval(df.loc[(df["featureGroup"] == featureGroup) & (df["model"] == "GB"), "bestParams"].reset_index(drop = True)[0])
    consRes = LOCO_GB(dataset, list(myDict.values())[1], list(myDict.values())[2], list(myDict.values())[0])
    df.loc[(df["featureGroup"] == featureGroup) & (df["model"] == "GB"), "LOCO_weightedAv"] = consRes

    # write back to the CSV
    df.to_csv("/user/home/uw20204/scratch/CScapeModels/encode/bestparamsnew2.txt", index = False)
