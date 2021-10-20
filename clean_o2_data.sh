#!/bin/bash

# Note that for this to work, especially the last step, you should use the login node to run this. 
cd /n/data1/bch/hemonc/camargo/li/DATA/o2_data
rm -r -f .git
git init
git remote add origin git@github.com:ShouWenWang/o2_data.git
git add -A
git commit -m 'Refresh git'
git push -f --set-upstream origin master