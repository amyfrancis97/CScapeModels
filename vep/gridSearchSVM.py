import numpy as np
import pandas as pd
import sklearn
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
    return clf.best_params_, clf.best_score_

consequences = consequences.rename(columns={"driver_stat": "class"})
aminoAcids = aminoAcids.rename(columns={"driver_stat": "class"})

file_path = '/user/home/uw20204/scratch/CScapeModels/vep/SVMbestParams.txt'
sys.stdout = open(file_path, "w")

print("consequences grid search: ")
gridSearch(consequences)

print("aminoAcids grid search: ")
gridSearch(aminoAcids)

print("merged grid search: ")
aminoAcids = aminoAcids.iloc[:, 3:]
merged = pd.concat([consequences, aminoAcids], axis=1)
gridSearch(merged)
