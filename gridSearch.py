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

# reading in the feature group in CSV format

dataset = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/TF+ChIP-seq/released/2021-10-21/data/training/coding/" +
"CSV.txt", header=None)

# reading in the corresponding variants to get genomic positions/ which chromosome the variants is on

variants = pd.read_csv("/user/home/uw20204/mrcieu_data/encode/public/TF+ChIP-seq/released/2021-10-21/data/training/coding/" +
"variants_incl.txt", names = ["chrom", "pos", "driverStat", "refAllele", "altAllele"], header=None, sep = "\t")

variants = variants.iloc[range(0, len(dataset)), :]

dataset_variants = pd.concat([variants, dataset], axis=1)
rows_with_nan = [index for index, row in dataset_variants.iterrows() if row.isnull().any()]
dataset_variants = dataset_variants.drop(labels=rows_with_nan, axis=0)

dataset = dataset_variants.rename(columns={0: "class"})

result = pd.concat([dataset[dataset["class"] == 1].sample(1000), dataset[dataset["class"] == -1].sample(1000)])
dataset = shuffle(result)
dataset = dataset.reset_index(drop = True)
print(dataset.head())

X = dataset.drop(["class", "altAllele", "refAllele", "chrom", "pos", "driverStat"], axis=1)
y = dataset["class"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

# Set the parameters by cross-validation
tuned_parameters = [
    {"kernel": ["rbf"], "gamma": [1e-3, 1e-4], "C": [1, 10, 100, 1000]},
    {"kernel": ["linear"], "C": [0.001, 0.01, 0.1, 1,  10, 100, 1000, 10000, 100000]},
]

scores = ["precision", "recall"]

for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

    clf = GridSearchCV(SVC(), tuned_parameters, scoring="%s_macro" % score)
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_["mean_test_score"]
    stds = clf.cv_results_["std_test_score"]
    for mean, std, params in zip(means, stds, clf.cv_results_["params"]):
        print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))
    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print()
