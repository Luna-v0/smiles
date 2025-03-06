#!/bin/bash

cd data

datasets=(
    "tox21"
    "muv"
    "hiv"
    "bace"
    "bbbp"
    "sider"
    "clintox"
    "delaney"
    "lipo"
    "sampl"
)

for dataset in "${datasets[@]}"; do
    wget -O "${dataset}.csv.gz" "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/${dataset}.csv.gz"
    gunzip "${dataset}.csv.gz" 
    rm "${dataset}.csv.gz"
done

echo "All datasets downloaded and extracted!"
