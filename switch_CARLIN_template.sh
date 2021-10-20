#!/bin/bash
template=$1
root_dir=$2

if $use_cCARLIN
then
        rsync -avP ${root_dir}/Custom_CARLIN/ ${root_dir}/cCARLIN/
        cp ${root_dir}/cCARLIN/CARLIN_def_cCARLIN.m ${root_dir}/cCARLIN/CARLIN_def.m
else
        rsync -avP ${root_dir}/Custom_CARLIN/ ${root_dir}/Tigre_CARLIN/
        cp ${root_dir}/Tigre_CARLIN/CARLIN_def_Tigre.m ${root_dir}/Tigre_CARLIN/CARLIN_def.m
fi