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

df = pd.read_csv("/user/home/uw20204/mrcieu_data/ensembl/public/vep/VEP_training_coding_csv_transcriptNoIncl.txt", sep = "\t")
print(df.head())
dataset = df.rename(columns={"driver_stat": "class"})
result = pd.concat([dataset[dataset["class"] == 1].sample(3000), dataset[dataset["class"] == -1].sample(3000)])
dataset = shuffle(result)
dataset = dataset.reset_index(drop = True)
print(dataset.head())

X = dataset.drop(["class", "chrom"], axis=1)
y = dataset["class"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

# Set the parameters by cross-validation
p_test3 = {'learning_rate':[0.15,0.1,0.05,0.01,0.005,0.001], 'n_estimators':[100,250,500,750,1000,1250,1500,1750]}

for i in range(1, 10):
    tuning = GridSearchCV(estimator =GradientBoostingClassifier(max_depth=4, min_samples_split=2, min_samples_leaf=1, subsample=1,max_features='sqrt', random_state=10), 
            param_grid = p_test3, scoring='accuracy',n_jobs=4,iid=False, cv=5)
    tuning.fit(X_train,y_train)
    print("gBoost params")
    print(tuning.cv_results_, tuning.best_params_, tuning.best_score_)
