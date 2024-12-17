#!/bin/bash

sample=$1
data_path=$2

cd $data_path

echo "---------sample: $sample -----------"
echo '' > $sample.split_locus.tsv
X0=`wc -l  $sample.trimmed.pear.assembled.fastq`
locus='all'
echo "$locus,$X0,$sample"
echo "$locus,$X0,$sample" >> $sample.split_locus.tsv


primer='GCAACTAGAAGGCACAGTCG' # 3'
#primer='CCGCTTACTTGTACAGCT' # 5'
locus='cCARLIN'
grep -A 2 -B 1 $primer $sample.trimmed.pear.assembled.fastq  | sed '/^--$/d' > $sample.$locus.trimmed.pear.assembled.fastq
echo "cCARLIN"
X0=`wc -l $sample.$locus.trimmed.pear.assembled.fastq`
echo "$locus,$X0,$sample"
echo "$locus,$X0,$sample" >> $sample.split_locus.tsv

primer='GCAACTAGAAGGCACCGACAA' #3'
#primer='GGCGATTCGCGAGGTACC' #5'
locus='Tigre'
grep -A 2 -B 1 $primer $sample.trimmed.pear.assembled.fastq  | sed '/^--$/d' > $sample.$locus.trimmed.pear.assembled.fastq
echo "Tigre"
X0=`wc -l $sample.$locus.trimmed.pear.assembled.fastq`
echo "$locus,$X0,$sample"
echo "$locus,$X0,$sample" >> $sample.split_locus.tsv


primer='GCAACTAGAAGGCACACAGCA' # 3'
#primer='GGCGCGGCCGCTTTACTT' #5'
locus='Rosa'
#grep $primer $sample.trimmed.pear.assembled.fastq > $sample.$locus.trimmed.pear.assembled.fastq
grep -A 2 -B 1 $primer $sample.trimmed.pear.assembled.fastq  | sed '/^--$/d' > $sample.$locus.trimmed.pear.assembled.fastq
echo "Rosa"
X0=`wc -l $sample.$locus.trimmed.pear.assembled.fastq`
echo "$locus,$X0,$sample"
echo "$locus,$X0,$sample" >> $sample.split_locus.tsv