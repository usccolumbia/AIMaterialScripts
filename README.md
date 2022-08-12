# AIMaterialScripts

## Machine learning model for ion conductivity prediction

   2022. Summer

### Jessica Chen and Jeffrey Hu

by <a href="http://mleg.cse.sc.edu" target="_blank">Machine Learning and Evolution Laboratory</a>, University of South Carolina


### Description of algorithm

Structural Features calculated: 

Calculated features and datasets were based on research from [Paper](https://doi.org/10.1039/C6EE02697D)

1. Feature1...
2. 


Machine learning model: 


### Datasets for training

ICSD-mix dataset (52317 samples)



All above datasets can be downloaded from [Figshare](https://figshare.com/articles/dataset/MT_dataset/20122796)
reference to Stanford paper.. (where we get the data)

### Model performance

report MAE R2 performance of model..


### How to train with your own dataset

#### Installation
1. Create your own conda or other enviroment.
2. install basic packages
```
pip install -r requirements.txt
```
3. Install `pytorch` from [pytorch web](https://pytorch.org/get-started/previous-versions/) given your python & cuda version


#### Data preparation
Download datasets from the above link, then unzip it under `MT_dataset` folder.
After the above, the directory should be:


```

```
#### Training 
An example is to train a MT-GPT model on the Hybrid-mix dataset. 
```
python ./MT_model/MT_GPT/train_GPT.py  --tokenizer ./MT_model/tokenizer/   --train_data  ./MT_Dataset/hy_mix/train.txt  --valid_data ./MT_Dataset/hy_mix/valid.txt  --output_dir ./output
```
The training for other models is similar to MT-GPT.

### Predicting ion conductivity

```
python ./MT_model/MT_GPT/train_GPT.py  --tokenizer ./MT_model/tokenizer/   --train_data  ./MT_Dataset/hy_mix/train.txt  --valid_data ./MT_Dataset/hy_mix/valid.txt  --output_dir ./output
```


