#!/bin/bash

###########################################################
# Perform all steps of LD score regression for ArchR scATAC
###########################################################

# Stop execution if we get an error along the way
#set -o errexit 


# Source config file
source ${HOME}/git_clones/scScalpChromatin/GWAS_scripts/LDSC/00_ldsc_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/boberrey/miniconda3/etc/profile.d/conda.sh
conda activate ldsc


############################################
# Perform cell type-specific analysis
############################################

# Perform regression on all .sumstat.gz files:
if [ ! -d ${cts_result_dir} ]; then mkdir ${cts_result_dir}; fi


for sumstat_f in ${sumstat_dir}/*.sumstats
do
    # Get the prefix for output from sumstats file
    trait=$(basename $sumstat_f | sed 's/.sumstats//')
    resultfile=${cts_result_dir}/${trait}.cell_type_results.txt
    if [ ! -f ${resultfile} ]
    then
        echo "Computing cell-type-specific partitioned heritability for ${trait}..."
        logfile=${cts_result_dir}/${trait}.h2-cts.out

        #Calculate cell-type-specific partitioned heritability
        sbatch -p wjg,biochem,sfgf -t 00:30:00 --mem=15G --cpus-per-task=1 \
        --job-name=${trait}_h2 --output=${logfile} --error=${logfile} \
        --wrap "${HOME}/git_clones/ldsc/ldsc.py --h2-cts ${sumstat_f} \
        --ref-ld-chr ${baseline_dir}/baselineLD. --out ${cts_result_dir}/${trait} \
        --ref-ld-chr-cts ${ldcts_file} --w-ld-chr ${weights_dir}/weights.hm3_noMHC."
        sleep 0.25s
    else
        echo "h2-cts results file ${resultfile} already exists. Skipping..."
    fi
done
