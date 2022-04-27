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
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import GridSearchCV
from config import *

# reading in the feature group in CSV format
dataset = dataset.reset_index(drop = True)
merged_data= dataset2.merge(dataset, on=["chrom","pos"])
merged_data = merged_data.iloc[:, [0,1,2,11,6]]

dataset = merged_data
rows_with_nan = [index for index, row in dataset.iterrows() if row.isnull().any()]
dataset = dataset.drop(labels=rows_with_nan, axis=0)

############## SVM CLASSIFIER #################

dfWeightedAvPrec2 = []
dfBalancedAcc2 = []
for i in range(1, 30):
    dfWeightedAvPrec1 = []
    dfBalancedAcc1 = []
    for chrom in range(1, 22):
        
        heldOutChrom = "chr" + str(chrom)

        # LOCO-CV for loop

        # select stated chrom to hold out as the test set
        datasetTest = dataset[dataset["chrom"] == heldOutChrom]
        
        X_test = datasetTest.drop(["driver_status_x", "chrom", "pos"], axis=1) # test dataset X are the signal values ONLY
        y_test = datasetTest["driver_status_x"] # test dataset y, or labels, are the driver status ONLY

        # hold out stated chrom to generate training dataset
        datasetTrain = dataset[dataset["chrom"] != heldOutChrom]

        # randomly sample 3000 +ive and 3000 -ive examples from the training dataset to carry forward
        datasetTrain = pd.concat([datasetTrain[datasetTrain["driver_status_x"] == 1].sample(1000), datasetTrain[datasetTrain["driver_status_x"] == -1].sample(1000)])
        datasetTrain = shuffle(datasetTrain)
        datasetTrain = datasetTrain.reset_index(drop = True)

        X_train = datasetTrain.drop(["driver_status_x", "chrom", "pos"], axis=1) # train dataset X are the signal values ONLY
        y_train = datasetTrain["driver_status_x"] # train dataset y, or labels, are the driver status ONLY

        # specify the model parameters for training, based on previous grid search
        #svclassifier = SVC(kernel='rbf', C = 100, gamma = 0.001)
        svclassifier = SVC(kernel='linear', C = 100000)
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
print("SVM weighted average:")
print(featureGroup)
print(round(Average(dfWeightedAvPrec2), 3))
weightedAvFinal = Average(dfWeightedAvPrec2)

