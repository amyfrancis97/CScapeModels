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

def gridSearchSVM(dataset):
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
    clf = GridSearchCV(SVC(), tuned_parameters, scoring="balanced_accuracy")
    clf.fit(X_train, y_train)
    return clf.best_params_, clf.best_score_

def gridSearchGB(dataset):
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
    best_params_learning_est = tuning.best_params_
    learningRateNoEstParams = list(tuning.best_params_.values())

    # gradient boosting grid search for max depth
    # Set the parameters by cross-validation
    p_test2 = {'max_depth':[2,3,4,5,6,7]}

    tuning = GridSearchCV(estimator =GradientBoostingClassifier(learning_rate=learningRateNoEstParams[0],n_estimators=learningRateNoEstParams[1], min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
            param_grid = p_test2, scoring='accuracy',n_jobs=4, cv=5)
    tuning.fit(X_train,y_train)
    bestParams = tuning.best_params_.copy()
    bestParams.update(best_params_learning_est)
    return bestParams, tuning.best_score_


# SVM grid search
dfFin = pd.DataFrame(columns = ["overallFeatureCat", "featureGroup", "model", "bestParams", "PA_bestScore","LOCO_weightedAv"])
consequences = consequences.rename(columns={"driver_stat": "class"})
aminoAcids = aminoAcids.rename(columns={"driver_stat": "class"})
VEPConsqRes = gridSearchSVM(consequences)
df = pd.DataFrame({"overallFeatureCat": "vep", "featureGroup": ["consequences"], "model": ["SVM"], "bestParams":[VEPConsqRes[0]], "PA_bestScore": [round(VEPConsqRes[1], 4)],"LOCO_weightedAv": ["N/A"]})
dfFin = dfFin.append(df)
VEPaaRes = gridSearchSVM(aminoAcids)
df = pd.DataFrame({"overallFeatureCat": "vep", "featureGroup": ["aminoAcids"], "model": ["SVM"], "bestParams":[VEPaaRes[0]], "PA_bestScore": [round(VEPaaRes[1], 4)],"LOCO_weightedAv": ["N/A"]})
dfFin = dfFin.append(df)
aminoAcids2 = aminoAcids.iloc[:, 3:]
merged = pd.concat([consequences, aminoAcids2], axis=1)
VEPaaConsqRes = gridSearchSVM(merged)
df = pd.DataFrame({"overallFeatureCat": "vep", "featureGroup": ["aminoAcidsConsMerged"], "model": ["SVM"], "bestParams":[VEPaaConsqRes[0]], "PA_bestScore": [round(VEPaaConsqRes[1], 4)],"LOCO_weightedAv": ["N/A"]})
dfFin = dfFin.append(df)

# gradient boosting grid search
dfFin = dfFin
VEPConsqRes = gridSearchGB(consequences)
df = pd.DataFrame({"overallFeatureCat": "vep", "featureGroup": ["consequences"], "model": ["GB"], "bestParams":[VEPConsqRes[0]], "PA_bestScore": [round(VEPConsqRes[1], 4)],"LOCO_weightedAv": ["N/A"]})
dfFin = dfFin.append(df)
VEPaaRes = gridSearchGB(aminoAcids)
df = pd.DataFrame({"overallFeatureCat": "vep", "featureGroup": ["aminoAcids"], "model": ["GB"], "bestParams":[VEPaaRes[0]], "PA_bestScore": [round(VEPaaRes[1], 4)],"LOCO_weightedAv": ["N/A"]})
dfFin = dfFin.append(df)
merged = pd.concat([consequences, aminoAcids2], axis=1)
VEPaaConsqRes = gridSearchGB(merged)
df = pd.DataFrame({"overallFeatureCat": "vep", "featureGroup": ["aminoAcidsConsMerged"], "model": ["GB"], "bestParams":[VEPaaConsqRes[0]], "PA_bestScore": [round(VEPaaConsqRes[1], 4)],"LOCO_weightedAv": ["N/A"]})
dfFin = dfFin.append(df)

dfFin.to_csv("/user/home/uw20204/scratch/CScapeModels/vep/bestparamsnew.txt", index = False)
