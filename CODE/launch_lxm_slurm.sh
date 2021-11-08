#!/usr/bin/env bash

cd /t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI

module add R-cbrg/current

export RSCRIPT=/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI/CODE
export PHENOTYPES=/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI/DATA/pheno.list

for pheno in $(cat $PHENOTYPES)
do
sbatch -p batch -c 8 --job-name lmx_${pheno} -o %j.out -e %j.err --mail-user=mfernand --mail-type=ALL \
--wrap="Rscript --vanilla $RSCRIPT/f_008_groupedKcv_v2.R \"${pheno}\" 8"
done
