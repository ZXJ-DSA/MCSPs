#!/bin/bash
data_path=$1
for strategy in wave ave pri
do
   CTR 256 1 $strategy DC
   USA 256 1 $strategy DC
   Orkut 1024 1 $strategy DC
done