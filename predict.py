# This code takes the a cif file and predicts the ion conductivity

from pymatgen.io.cif import CifParser # get cifs
from crystal_feature import CrystalFeature # feature python code
import sys # sys.argv[1]
import pandas as pd 

if len(sys.argv)<=1:
     print(f"python {sys.argv[0]} your.cif")
     exit(0)

crystal = CrystalFeature(f'{sys.argv[1]}')
dictionary = crystal.get_all_features() # this only gives the values of the features, not the feature names
df = pd.DataFrame()
for key, val in dictionary.items():
	df[key] = [val]
print(df)

import joblib
# load, no need to initialize the loaded_rf
loaded_rf = joblib.load("./random_forest.joblib") # if i run this on terminal I am in the folder 

print(loaded_rf.predict(df))
