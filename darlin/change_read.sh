#!/bin/bash
root_dir=$1
DataList_DNA=$2
ini_cut_length=$3
end_cut_length=$4
UMI_length=$5

#DataList_DNA='CC_DNA_poor,CC_RNA_poor'
#cut_length=40

cd $root_dir
Field_Separator=$IFS
# set comma as internal field separator for the string list
IFS=,

command="\$ln++; s/^(.{$UMI_length})(.{$ini_cut_length})(.*)(.{$end_cut_length})\$/\$1\$3/ if (\$ln % 2 == 0)"
#command="\$ln++; s/(.{40})\$// if (\$ln % 2 == 0)"
echo $command

for val in $DataList_DNA;
do
echo $val
#cat ${val}.trimmed.pear.assembled.fastq | perl -plane '$ln++; s/(.{40})$// if ($ln % 2 == 0)' > ${val}_mo.trimmed.pear.assembled.fastq 
cat ${val}.trimmed.pear.assembled.fastq | perl -plane "$command" > ${val}_mo.trimmed.pear.assembled.fastq 

done

