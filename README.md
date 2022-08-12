# AIMaterialScripts

## Machine learning model for ion conductivity prediction

   2022. Summer

### Jessica Chen and Jeffrey Hu from DFHS

working at <a href="http://mleg.cse.sc.edu" target="_blank">Machine Learning and Evolution Laboratory</a>, University of South Carolina


### Description of algorithm

**Structural Features calculated:**

More information about calculated features as well as the training dataset can be found in [Paper](https://doi.org/10.1039/C6EE02697D). 

1. Average atomic volume, AAV 
2. Standard deviation in Li neighbor count, SDLC
3. Standard deviation in Li bonding ionicity, SDLI
4. Average Li bond ionicity, LBI
5. Average Li neighbor count, LNC
6. Average Li-Li bonds per Li, LLB
7. Average bond ionicity of sublattice, SBI
8. Average sublattice neighbor count, SNC
9. Anion framework coordination, AFC
10. Average shortest anion-anion separation distance, AASD
11. Volume per anion, VPA
12. Average shortest Li-anion separation distance, LASD
13. Average shortest Li-Li separation distance, LLSD
14. Average electronegativity of sublattice, ENS
15. Packing fraction of full crystal, PF
16. Packing fraction of sublattice, SPF
17. Average straight-line path width, SLPW
18. Average straight-line path electronegativity, SLPE
19. Ratio of average Li bond ionicity to average sublattice bond ionicity, RBI
20. Ratio of average Li neighbor count to average sublattice neighbor count, RNC


**Machine learning model:**

The machine learning model used to train the data (from [Paper](https://doi.org/10.1039/C6EE02697D)) was a random forest classifier. 
// I don't know what else to add


### Datasets for training

ion conductivity dataset (32 samples)
//this one we didnt train with icsd we trained with the 30 things



All above datasets can be downloaded from [Figshare](https://figshare.com/articles/dataset/MT_dataset/20122796)
reference to Stanford paper.. (where we get the data)

### Model performance

report MAE R2 performance of model..


### How to train with your own dataset

#### Installation
1. Create your own conda or other enviroment. // is this required for ours?
2. install basic packages
```
pip install sklearn
pip install pymatgen

pip install pandas
pip install numpy
```
3. Install `pytorch` from [pytorch web](https://pytorch.org/get-started/previous-versions/) given your python & cuda version // is this required for ours?


#### Data preparation
Download IonDatasets folder from files above.
If predicting more than one file, create a csv in the `IonDatasets` folder  with the Filename and Composition organized as:
| Filename (cif) | Composition |
| ----------- | ----------- |
| LiI.cif | LiI |


#### Training 
An example is to train a MT-GPT model on the Hybrid-mix dataset. 
```
python train.py
```
The training for other models is similar to MT-GPT.

### Predicting ion conductivity

```
python predict.py test.cif
```


