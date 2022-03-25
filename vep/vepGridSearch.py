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
from sklearn.metrics import classification_report

df = pd.read_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/trainingVEP2.txt", sep = "\t")
print(df.head())
dataset = df.rename(columns={"driver_stat": "class"})
result = pd.concat([dataset[dataset["class"] == 1].sample(1000), dataset[dataset["class"] == -1].sample(1000)])
dataset = shuffle(result)
dataset = dataset.reset_index(drop = True)
print(dataset.head())

#X = dataset.drop(["class", "chrom"], axis=1)

X = dataset.drop(["class"], axis=1)
y = dataset["class"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

# Set the parameters by cross-validation
tuned_parameters = [
    {"kernel": ["rbf"], "gamma": [1e-3, 1e-4], "C": [1, 10, 100, 1000]},
    {"kernel": ["linear"], "C": [0.001, 0.01, 0.1, 1,  10, 100, 1000, 10000, 100000]},
]

for i in range(1, 10):
    print("# Tuning hyper-parameters for %s" % "balanced_accuracy")
    print()

    clf = GridSearchCV(SVC(), tuned_parameters, scoring="balanced_accuracy")
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
