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
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import GradientBoostingClassifier

df = pd.read_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/VEP_training_coding_csv_transcriptNoIncl.txt", sep = "\t")
print(df.head())
dataset = df.rename(columns={"driver_stat": "class"})
result = pd.concat([dataset[dataset["class"] == 1].sample(5000), dataset[dataset["class"] == -1].sample(5000)])
dataset = shuffle(result)
dataset = dataset.reset_index(drop = True)
print(dataset.head())

X = dataset.drop(["class", "chrom"], axis=1)
y = dataset["class"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

cs = [10, 20, 30, 40, 50, 70, 100, 200]
for c in cs:

    gb_clf = GradientBoostingClassifier(learning_rate=1, n_estimators=c,max_depth=3, min_samples_split=2, min_samples_leaf=1, subsample=1)
    gb_clf.fit(X_train, y_train) # fit the model
    y_pred = gb_clf.predict(X_test) # validate the model using all of the available data for the left out chromosome

    #Import scikit-learn metrics module for accuracy calculation
    from sklearn import metrics

    # Model Accuracy: how often is the classifier correct?
    print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

    # Model Precision: what percentage of positive tuples are labeled as such?
    print("Precision:",metrics.precision_score(y_test, y_pred))

    # Model Recall: what percentage of positive tuples are labelled as such?
    print("Recall:",metrics.recall_score(y_test, y_pred))