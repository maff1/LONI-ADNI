#!/usr/bin/env bash

cd /t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI

module add R-cbrg/current

export RSCRIPT=/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI/CODE
export METABOLITES=/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI/DATA/metabolites.list

for met in $(cat $METABOLITES)
do
sbatch -p batch -c 4 --job-name lmx_${met} -o %j.out -e %j.err --mail-user=mfernand --mail-type=ALL \
--wrap="Rscript --vanilla $RSCRIPT/f_011_bootstrapping.R \"MMSE\" \"${met}\" 4 1000"
done
