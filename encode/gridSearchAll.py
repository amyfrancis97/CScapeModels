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

    # gradient boosting grid search for learning rate & no. estimators
    # Set the parameters by cross-validation
    p_test3 = {'learning_rate':[0.15,0.1,0.05,0.01,0.005,0.001], 'n_estimators':[100,250,500,750,1000,1250,1500,1750]}

    tuning = GridSearchCV(estimator =GradientBoostingClassifier(max_depth=4, min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
                param_grid = p_test3, scoring='accuracy',n_jobs=4, cv=5)
    tuning.fit(X_train,y_train)
    print("gradient boosting grid search results:")
    print(tuning.best_params_, tuning.best_score_)
    learningRateNoEstParams = list(tuning.best_params_.values())

    # gradient boosting grid search for max depth
    # Set the parameters by cross-validation
    p_test2 = {'max_depth':[2,3,4,5,6,7]}

    tuning = GridSearchCV(estimator =GradientBoostingClassifier(learning_rate=learningRateNoEstParams[0],n_estimators=learningRateNoEstParams[1], min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
            param_grid = p_test2, scoring='accuracy',n_jobs=4,iid=False, cv=5)
    tuning.fit(X_train,y_train)
    print("optimum maxdepth gboost")
    print(tuning.best_params_, tuning.best_score_)

# call the grid search function for the merged dataset
# reading in the feature group in CSV format


# reformat to match criteria for grid search function
variants = variants.iloc[range(0, len(dataset)), :]
dataset_variants = pd.concat([variants, dataset], axis=1)
dataset_variants = dataset_variants.drop(columns = ["pos", "refAllele", "altAllele", 0])
dataset = dataset_variants.rename(columns={"driverStat": "class"})

file_path = '/user/home/uw20204/scratch/CScapeModels/encode/bestParams.txt'
sys.stdout = open(file_path, "w")


print("encode grid search: ")
gridSearch(dataset)
