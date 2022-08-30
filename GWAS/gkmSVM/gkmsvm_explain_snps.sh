#!/bin/bash

# Stop execution if we get an error along the way
#set -o errexit 

############################################
# Run gkmexplain on ref and alt SNPs 
############################################

# https://github.com/kundajelab/lsgkm-svr
# https://github.com/kundajelab/gkmexplain

# Load required modules
ml system gcc

# Add executibles to path
export PATH=$HOME/git_clones/lsgkm-svr/bin:$PATH

# Directory containing fit models
model_dir=/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/full_models_1000bp 

# Directory for prediction results
result_dir=/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/snp_explanations
if [ ! -d ${result_dir} ]; then mkdir ${result_dir}; fi

# Fasta files to predict
fasta_dir=/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/snp_fastas
snp_fastas=(${fasta_dir}/*snp_seqs*.fasta)

# All previously fit models
model_files=(${model_dir}/*.model.txt)
file_basenames=("${model_files[@]##*/}")
model_headers=("${file_basenames[@]/.model.txt/}")

#############################################################
# Predict accessibility for all SNP fastas using all models
#############################################################

for model in "${model_headers[@]}"
do
    # Prepare filepaths
    model_file=${model_dir}/${model}.model.txt
    job_header=gkm_250exp
    outpre=${result_dir}/${model}

    for fasta in "${snp_fastas[@]}"
    do
        # Prepare output file
        fasta_basename="${fasta##*/}"
        fasta_basename="${fasta_basename/.fasta/}"
        result_file=${outpre}.${fasta_basename}.gkmexplain.txt
        log_file=${outpre}.${fasta_basename}.log
        if [ ! -f ${result_file} ]
        then
            # Explain model importance
            echo "Using model ${model} to calculate gkmexplain scores for ${fasta_basename}..."
            sbatch -p wjg,biochem,sfgf -t 23:59:00 --mem=12G --cpus-per-task=1 \
            --job-name=${job_header} --output=${log_file} --error=${log_file} \
            --wrap "gkmexplain ${fasta} ${model_file} ${result_file}"
        else
            echo "File ${result_file} already exists! Skipping..."
        fi
    done
done



