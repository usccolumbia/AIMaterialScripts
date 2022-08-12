# This code trains the data from Ion_Data_Featurized.csv using the random forest model

#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from crystal_feature import CrystalFeature
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LinearRegression
from numpy import mean
from numpy import absolute
from numpy import sqrt

from sklearn.preprocessing import StandardScaler


# In[2]:


df = pd.read_csv("Ion_Data_Featurized.csv")
print(df.columns)
y = df["Ionic_Conductivity"]
for index, v in enumerate(y):
    if v >= 10**-4:
        df.iloc[index] = 1
    else:
        df.iloc[index] = 0
print(y)
        
a = ['aav', 'sdlc', 'sdli','lbi', 'lnc', 'llb', 'sbi', 'snc', 'afc', 'aasd', 'vpa', 'lasd', 'llsd', 
     'ens', 'pf', 'spf', 'slpw', 'slpe', 'rbi', 'rnc']
X = df[a]
print(X)
print(y)


# In[4]:


clf = RandomForestClassifier(max_depth=2, random_state=0)
# clf.fit(X, y)
# print(clf.predict(X.iloc[0:]))

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)
clf.fit(X_train, y_train)

loocv = LeaveOneOut()
scores = cross_val_score(clf, X, y, cv=loocv)
print(f"Accuracy: {mean(scores.mean()*100.0)}%")

# sc = StandardScaler()
# X_train = sc.fit_transform(X_train)
# X_test = sc.transform(X_test)


# In[5]:


import joblib
joblib.dump(clf, "./random_forest.joblib")


# In[ ]:




