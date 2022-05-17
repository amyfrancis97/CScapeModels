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
df = pd.read_csv("bestparamsnew.txt")

consequences = consequences.rename(columns={"driver_stat": "class"})
aminoAcids = aminoAcids.rename(columns={"driver_stat": "class"})

# carrying out function for LOCO-CV SVM and writing into parameter file.
# pull out list of best params for each model

# carry out SVM models
myDict = ast.literal_eval(df.loc[(df["featureGroup"] == "consequences") & (df["model"] == "SVM"), "bestParams"].reset_index(drop = True)[0])
consRes = LOCO_SVM(consequences, list(myDict.values())[2], list(myDict.values())[0], list(myDict.values())[1])
df.loc[(df["featureGroup"] == "consequences") & (df["model"] == "SVM"), "LOCO_weightedAv"] = consRes

myDict = ast.literal_eval(df.loc[(df["featureGroup"] == "aminoAcids") & (df["model"] == "SVM"), "bestParams"].reset_index(drop = True)[0])
AARes = LOCO_SVM(aminoAcids, list(myDict.values())[2], list(myDict.values())[0], list(myDict.values())[1])
df.loc[(df["featureGroup"] == "aminoAcids") & (df["model"] == "SVM"), "LOCO_weightedAv"] = AARes

aminoAcids2 = aminoAcids.iloc[:, 3:]
aminoAcidsConsMerged = pd.concat([consequences, aminoAcids2], axis=1)
myDict = ast.literal_eval(df.loc[(df["featureGroup"] == "aminoAcidsConsMerged") & (df["model"] == "SVM"), "bestParams"].reset_index(drop = True)[0])
mergedRes = LOCO_SVM(aminoAcidsConsMerged, list(myDict.values())[1], list(myDict.values())[0], 'auto')
df.loc[(df["featureGroup"] == "aminoAcidsConsMerged") & (df["model"] == "SVM"), "LOCO_weightedAv"] = AARes

# carry out GB models
myDict = ast.literal_eval(df.loc[(df["featureGroup"] == "consequences") & (df["model"] == "GB"), "bestParams"].reset_index(drop = True)[0])
consRes = LOCO_GB(consequences, list(myDict.values())[1], list(myDict.values())[2], list(myDict.values())[0])
df.loc[(df["featureGroup"] == "consequences") & (df["model"] == "GB"), "LOCO_weightedAv"] = consRes

myDict = ast.literal_eval(df.loc[(df["featureGroup"] == "aminoAcids") & (df["model"] == "GB"), "bestParams"].reset_index(drop = True)[0])
AARes = LOCO_GB(aminoAcids, list(myDict.values())[1], list(myDict.values())[2], list(myDict.values())[0])
df.loc[(df["featureGroup"] == "aminoAcids") & (df["model"] == "GB"), "LOCO_weightedAv"] = AARes

myDict = ast.literal_eval(df.loc[(df["featureGroup"] == "aminoAcidsConsMerged") & (df["model"] == "GB"), "bestParams"].reset_index(drop = True)[0])
mergedRes = LOCO_GB(aminoAcidsConsMerged, list(myDict.values())[1], list(myDict.values())[2], list(myDict.values())[0])
df.loc[(df["featureGroup"] == "aminoAcidsConsMerged") & (df["model"] == "GB"), "LOCO_weightedAv"] = AARes

# write back to the CSV
df.to_csv("/user/home/uw20204/scratch/CScapeModels/vep/bestparamsnewv3.txt", index = False)
