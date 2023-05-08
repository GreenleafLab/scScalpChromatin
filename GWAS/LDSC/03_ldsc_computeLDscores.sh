#!/bin/bash

############################################
# Compute LD scores with annot files
############################################

# Stop execution if we get an error along the way
#set -o errexit 

# Source config file
source ${HOME}/git_clones/scScalpChromatin/GWAS/LDSC/00_ldsc_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/boberrey/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

# Compute LD scores if .ldscore files are not present
for bf in ${cts_peak_dir}/*.bed
do
    prefix=$(basename ${bf} | sed 's/_specific_peaks.bed//')

    for i in {1..22} # Sex chromosomes cannot be used in ldsc presently
    do
        ldsfile=${ldscore_dir}/${peakClass}.${prefix}.${i}.l2.ldscore.gz
        annotfile=${ldscore_dir}/${peakClass}.${prefix}.${i}.annot.gz

        # Make sure annot file was previously generated
        if [ ! -f ${annotfile} ]; then echo "ERROR!: annot file ${annotfile} does not exist!"; exit 1; fi

        if [ ! -f ${ldsfile} ]
        then
            echo "Calculating LD scores for chr${i} of ${prefix}..."
            logfile=${ldscore_dir}/${prefix}.${i}.ldsc.out
            sbatch -p wjg,biochem,sfgf -t 01:00:00 --mem=15G --cpus-per-task=1 \
            --job-name=calc_LD --output=${logfile} --error=${logfile} \
            --wrap "${HOME}/git_clones/ldsc/ldsc.py \
            --l2 --bfile ${plink_dir}/1000G.EUR.hg38.${i} \
            --ld-wind-cm 1 --annot ${ldscore_dir}/${peakClass}.${prefix}.${i}.annot.gz \
            --thin-annot --out ${ldscore_dir}/${peakClass}.${prefix}.${i} \
            --print-snps ${snp_dir}/list.txt" 
            sleep 0.25s
        else
            echo "LD score file ${ldsfile} already exists. Skipping..."
        fi
    done
done