# grid search SVM
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.utils import shuffle
import sys
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import GradientBoostingClassifier
from config import *

# reading in the feature group in CSV format
dataset = dataset.reset_index(drop = True)
merged_data= dataset2.merge(dataset, on=["chrom","pos"])
merged_data = merged_data.iloc[:, [0,1,2,11,6]]
dataset = merged_data

# reformat to match criteria for grid search function
dataset = dataset.drop(columns = "pos")
dataset = dataset.rename(columns={"driver_status_x": "class"})

# function to perform multiple grid searches
def gridSearch(dataset):
    rows_with_nan = [index for index, row in dataset.iterrows() if row.isnull().any()]
    dataset = dataset.drop(labels=rows_with_nan, axis=0)

    result = pd.concat([dataset[dataset["class"] == 1].sample(1000), dataset[dataset["class"] == -1].sample(1000)])
    dataset = shuffle(result)
    dataset = dataset.reset_index(drop = True)

    X = dataset.drop(["class", "chrom"], axis=1)
    y = dataset["class"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

    # Set the parameters by cross-validation
    tuned_parameters = [
        {"kernel": ["rbf"], "gamma": [1e-3, 1e-4], "C": [1, 10, 100, 1000]},
        {"kernel": ["linear"], "C": [0.001, 0.01, 0.1, 1,  10, 100, 1000, 10000, 100000]},
    ]

    # svm grid search
    print("svm grid search results: ")
    print("# Tuning hyper-parameters for %s" % "balanced_accuracy")
    print()
    clf = GridSearchCV(SVC(), tuned_parameters, scoring="balanced_accuracy")
    clf.fit(X_train, y_train)
    print("Best parameters set found on development set:")
    print(clf.best_params_, clf.best_score_)


# call the grid search function for the merged dataset
file_path = '/user/home/uw20204/scratch/CScapeModels/conservation/bestParamsSVM.txt'
sys.stdout = open(file_path, "w")

print("phyloP & phastcons 30-way grid search: ")
gridSearch(dataset)
