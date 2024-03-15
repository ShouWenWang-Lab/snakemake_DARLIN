#!/bin/bash
template=$1
root_dir=$2

if $use_cCARLIN
then
        rsync -avP ${root_dir}/Custom_CARLIN/ ${root_dir}/cDARLIN/
        cp ${root_dir}/cDARLIN/CARLIN_def_cCARLIN.m ${root_dir}/cDARLIN/CARLIN_def.m
else
        rsync -avP ${root_dir}/Custom_CARLIN/ ${root_dir}/Tigre_DARLIN/
        cp ${root_dir}/Tigre_DARLIN/CARLIN_def_Tigre.m ${root_dir}/Tigre_DARLIN/CARLIN_def.m
fi
