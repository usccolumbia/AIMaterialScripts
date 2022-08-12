#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
df = pd.read_csv("IonDatasets/ion_conductivity.csv")
df


# In[2]:


from pymatgen.io.cif import CifParser
from crystal_feature import CrystalFeature


# In[3]:


crystal_features = []
all_crystalfeatures = []
count = 0
for file in df['Filename']:
    crystal = CrystalFeature(f'IonDatasets/{file}')
    crystal_features.append(crystal)
    all_crystalfeatures.append(crystal.get_all_features())
    print(f'{count} done')
    count += 1


# In[4]:


features = crystal_features[0].get_feature_names()
features


# In[5]:


for feature in features:
    feature_col = []
    for crystal in all_crystalfeatures:
        feature_col.append(crystal[feature])
    df[feature] = feature_col


# In[6]:


df


# In[7]:


df.to_csv("Ion_Data_Featurized.csv",index=None)


# In[ ]:




