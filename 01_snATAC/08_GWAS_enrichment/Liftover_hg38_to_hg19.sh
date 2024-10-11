#!/bin/sh
# liftOver hg38 to hg19

liftOver=./04_Softwares/liftOver/liftOver
chain=./04_Softwares/liftOver/hg38ToHg19.over.chain.gz

for samp in `ls *bed`
do
$liftOver $samp $chain ${samp%.*}".hg19.bed" ${samp%.*}".hg19.unlifted.bed"
done