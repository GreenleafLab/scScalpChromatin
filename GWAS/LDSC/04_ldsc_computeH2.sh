#!/bin/bash

# Stop execution if we get an error along the way
#set -o errexit 

# Source config file
source ${HOME}/git_clones/scScalpChromatin/GWAS_scripts/LDSC/00_ldsc_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/boberrey/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

############################################
# Compute Partitioned Heritability
############################################

# Perform regression on all .sumstat.gz files:

if [ ! -d ${h2_result_dir} ]; then mkdir ${h2_result_dir}; fi

for sumstat_f in ${sumstat_dir}/*.sumstats
do
    # Get the prefix for output from sumstats file
    trait=$(basename ${sumstat_f} | sed 's/.sumstats//')

    for bf in ${cts_peak_dir}/*.bed
    do
    prefix=$(basename ${bf} | sed 's/_specific_peaks.bed//')
    resultfile=${h2_result_dir}/${prefix}.${trait}.results
    if [ ! -f ${resultfile} ]
    then
        echo "Computing partitioned heritability for ${trait} in ${prefix}..."
        logfile=${h2_result_dir}/${prefix}.${trait}.h2.out

        # Calculate partitioned heritability
        sbatch -p wjg,biochem,sfgf -t 01:00:00 --mem=15G --cpus-per-task=1 \
        --job-name=${prefix}_h2 --output=${logfile} --error=${logfile} \
        --wrap "${HOME}/git_clones/ldsc/ldsc.py --h2 ${sumstat_f} \
        --ref-ld-chr ${ldscore_dir}/${peakClass}.${prefix}.,${baseline_dir}/baselineLD. \
        --w-ld-chr ${weights_dir}/weights.hm3_noMHC. --frqfile-chr ${frq_dir}/1000G.EUR.QC. \
        --out ${h2_result_dir}/${prefix}.${trait} --overlap-annot" 
        sleep 0.25s
    else
        echo "h2 results file ${resultfile} already exists. Skipping..."
    fi
    done
done
