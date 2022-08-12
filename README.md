# AIMaterialScripts


Ciation: Nihang Fu, Lai Wei, Yuqi Song, Qinyang Li, Rui Xin, Sadman Sadeed Omee, Rongzhi Dong, Edirisuriya M. Dilanga Siriwardane, Jianjun Hu.  Materials Transformer Language Models for Generative Materials Design: a Benchmark Study. Arxiv 2022

### Jessica Chen and Jeffrey Hu

by <a href="http://mleg.cse.sc.edu" target="_blank">Machine Learning and Evolution Laboratory</a>, University of South Carolina

2022. Summer.


### Benchmark Datasets for training inorganic materials composition transformers

ICSD-mix dataset (52317 samples)



All above datasets can be downloaded from [Figshare](https://figshare.com/articles/dataset/MT_dataset/20122796)

### Trained Materials Transformer Models



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

