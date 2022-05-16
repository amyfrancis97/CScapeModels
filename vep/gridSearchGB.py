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

    # gradient boosting grid search for learning rate & no. estimators
    # Set the parameters by cross-validation
    p_test3 = {'learning_rate':[0.15,0.1,0.05,0.01,0.005,0.001], 'n_estimators':[100,250,500,750,1000,1250,1500,1750]}

    tuning = GridSearchCV(estimator =GradientBoostingClassifier(max_depth=4, min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
                param_grid = p_test3, scoring='accuracy',n_jobs=4,iid=False, cv=5)
    tuning.fit(X_train,y_train)
    print("gradient boosting grid search results:")
    print(tuning.best_params_, tuning.best_score_)
    learningRateNoEstParams = list(tuning.best_params_.values())

    # gradient boosting grid search for max depth
    # Set the parameters by cross-validation
    p_test2 = {'max_depth':[2,3,4,5,6,7]}

    tuning = GridSearchCV(estimator =GradientBoostingClassifier(learning_rate=learningRateNoEstParams[0],n_estimators=learningRateNoEstParams[1], min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
            param_grid = p_test2, scoring='accuracy',n_jobs=4, cv=5)
    tuning.fit(X_train,y_train)
    print("optimum maxdepth gboost")
    print(tuning.best_params_, tuning.best_score_)

consequences = consequences.rename(columns={"driver_stat": "class"})
aminoAcids = aminoAcids.rename(columns={"driver_stat": "class"})

file_path = '/user/home/uw20204/scratch/CScapeModels/vep/GBbestParams.txt'
sys.stdout = open(file_path, "w")

print("consequences grid search: ")
gridSearch(consequences)

print("aminoAcids grid search: ")
gridSearch(aminoAcids)

print("merged grid search: ")
aminoAcids = aminoAcids.iloc[:, 3:]
merged = pd.concat([consequences, aminoAcids], axis=1)
gridSearch(merged)
