#!/bin/bash

############################################
# Create annot files for each set of peaks:
############################################

# Stop execution if we get an error along the way
#set -o errexit 

# Source config file
source ${HOME}/git_clones/scScalpChromatin/GWAS_scripts/LDSC/00_ldsc_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/boberrey/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

# Make ldscore directory if it doesn't already exist
if [ ! -d ${ldscore_dir} ]; then mkdir ${ldscore_dir}; fi

# Delete previous contrast file if it already exists
if [ -f ${ldcts_file} ]; then rm ${ldcts_file}; fi
    
# Generate annot files for each peak set for each chromosome
for bf in ${cts_peak_dir}/*.bed
do
    prefix=$(basename ${bf} | sed 's/_specific_peaks.bed//')

    for i in {1..22} # Sex chromosomes cannot be used in ldsc
    do
        annotfile=${ldscore_dir}/${peakClass}.${prefix}.${i}.annot.gz

        if [ ! -f ${annotfile} ]
        then
            echo "Making annot file for chr${i} of ${prefix}..."
            logfile=${ldscore_dir}/${prefix}.${i}.annot.out
            sbatch -p wjg,biochem,sfgf -t 00:30:00 --mem=10G --cpus-per-task=1 \
            --job-name=ldsc_annot --output=${logfile} --error=${logfile} \
            --wrap "${HOME}/git_clones/ldsc/make_annot_BOR.py \
            --bed-file ${bf} --bimfile ${plink_dir}/1000G.EUR.hg38.${i}.bim \
            --annot-file ${annotfile}"
            sleep 0.25s
        else
            echo "Annot file ${annotfile} alreay exists. Skipping..."
        fi
    done
    # Append to ldcts contrast file
    # This file will define what contrasts to run in ldsc --h2-cts (i.e. subgroup vs allPeaks)
    if [ ${prefix} != "allPeaks" ]
    then
        echo -e "${prefix}\t${ldscore_dir}/${peakClass}.${prefix}.,${ldscore_dir}/${peakClass}.allPeaks." >> ${ldcts_file}
    fi
done


############################################
