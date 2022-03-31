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
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import GradientBoostingClassifier

# reading in the feature group in CSV format
dataset = pd.read_csv("/user/home/uw20204/mrcieu_data/ucsc/public/hg38.phastCons30way/released/2021-10-27/" + "hg38.phastCons30way.trainingVariantsCoding.CSV")
dataset = dataset.reset_index(drop = True)

dataset2 = pd.read_csv("/user/home/uw20204/mrcieu_data/ucsc/public/hg38.phyloP30way/released/2021-10-27/" + "hg38.phyloP30way.trainingVariantsCoding.CSV")

merged_data= dataset2.merge(dataset, on=["chrom","pos"])
merged_data = merged_data.iloc[:, [0,1,2,11,6]]

dataset = merged_data
rows_with_nan = [index for index, row in dataset.iterrows() if row.isnull().any()]
dataset = dataset.drop(labels=rows_with_nan, axis=0)

result = pd.concat([dataset[dataset["driver_status_x"] == 1].sample(1000), dataset[dataset["driver_status_x"] == -1].sample(1000)])
dataset = shuffle(result)
dataset = dataset.reset_index(drop = True)

X = dataset.drop(["driver_status_x", "chrom", "pos"], axis=1) 
y = dataset["driver_status_x"] 

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

# Set the parameters by cross-validation
p_test3 = {'learning_rate':[0.15,0.1,0.05,0.01,0.005,0.001], 'n_estimators':[100,250,500,750,1000,1250,1500,1750]}

for i in range(1, 10):
    tuning = GridSearchCV(estimator =GradientBoostingClassifier(max_depth=4, min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
            param_grid = p_test3, scoring='accuracy',n_jobs=4,iid=False, cv=5)
    tuning.fit(X_train,y_train)
    print(tuning.cv_results_, tuning.best_params_, tuning.best_score_)
