#!/bin/bash

# Stop execution if we get an error along the way
#set -o errexit 

########################################################################################
# Train gkmSVM models on previously generated peak sequences and corresponding null seqs
########################################################################################

# See:
# https://github.com/Dongwon-Lee/lsgkm/
# https://github.com/kundajelab/lsgkm-svr

# Load required modules
ml system gcc

# Add executibles to path
export PATH=$HOME/git_clones/lsgkm-svr/bin:$PATH

# Get output directory ready
model_dir=/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/fit_models_1000bp
full_model_dir=/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/full_models_1000bp
if [ ! -d ${model_dir} ]; then mkdir ${model_dir}; fi
if [ ! -d ${full_model_dir} ]; then mkdir ${full_model_dir}; fi

# Find all peak categories and fastas
fasta_dir=/oak/stanford/groups/wjg/boberrey/hairATAC/scratch_copy/scratch/analyses/GWAS/gkmSVM/fastas_1000bp_randOnly
chr1_true_seq_files=(${fasta_dir}/*chr1_true_seqs.fasta)
file_basenames=("${chr1_true_seq_files[@]##*/}")
model_headers=("${file_basenames[@]/_chr1_true_seqs.fasta/}")

# Model training parameters
yval="0"        #    
tval="5"        # Kernel (5 = wgkmrbf)
lval="11"       # word length
kval="7"        # number of informative positions
dval="3"        # number of mismatches to consider
gval="2.0"      # Recommended to be set from 1.0 to 2.0 for RBF kernels (https://github.com/Dongwon-Lee/lsgkm)
cval="10.0"     # Recommended to be set from 1.0 to 10.0 for RBF kernels
Mval="50"       # initial value for exponential decay function
Hval="50"       # half-life parameter
eval="0.001"    # precision
mval="20000.0"  # Cache memory size in MB

# Define training/testing folds
all_chrs=($(seq 1 1 22))
all_chrs+=("X")
all_chrs=("${all_chrs[@]/#/chr}")

# 10-fold cross-validation
fold0=(chr1)
fold1=(chr2 chr19)
fold2=(chr3 chr20)
fold3=(chr6 chr13 chr22)
fold4=(chr5 chr16)
fold5=(chr4 chr15 chr21)
fold6=(chr7 chr14 chr18)
fold7=(chr11 chr17 chrX)
fold8=(chr9 chr12)
fold9=(chr8 chr10)
n_folds=($(seq 0 1 9))

all_folds=("${n_folds[@]/#/fold}")

# Directory to hold training data
train_dir=${fasta_dir}/training_fastas
if [ ! -d ${train_dir} ]; then mkdir ${train_dir}; fi

############################################
# Submit jobs
############################################

job_header=gkm_train

for model in "${model_headers[@]}"
do
    
    # First, fit cross-validation models for assessing model performance

    for fold in "${all_folds[@]}"
    do
        # Define fold chromosomes
        fchr=${fold}[@]
        fchr=("${!fchr}")

        # Get chromosomes to use for training
        train_chr=( $(printf "%s\n" "${all_chrs[@]}" "${fchr[@]}" | sort | uniq -u) )

        echo "Training chromosomes for ${fold}:"

        # Get true and null seq fasta files for training
        true_seq_files=( $(printf "${fasta_dir}/${model}_%s_true_seqs.fasta " "${train_chr[@]}") )
        null_seq_files=( $(printf "${fasta_dir}/${model}_%s_null_seqs.fasta " "${train_chr[@]}") )

        outpre=${model_dir}/${model}.${fold}
        resultfile=${outpre}.model.txt

        if [ ! -f ${resultfile} ]
        then
            # Create files for combined training data
            echo "Combining training fastas..."
            true_seqs=${train_dir}/${model}.${fold}.true_seqs.fasta
            null_seqs=${train_dir}/${model}.${fold}.null_seqs.fasta
            cat "${true_seq_files[@]}" > ${true_seqs}
            cat "${null_seq_files[@]}" > ${null_seqs}

            # Train models
            echo "Training model for ${model}..."
            sbatch -p wjg,biochem,sfgf -t 72:00:00 --mem=40G --cpus-per-task=16 \
            --job-name=${job_header} --output=${outpre}.log --error=${outpre}.log \
            --wrap "gkmtrain -y ${yval} -t ${tval} -l ${lval} -k ${kval} -d ${dval} -g ${gval} -c ${cval} -M ${Mval} -H ${Hval} -e ${eval} \
            -m ${mval} -r 123 -T 16 ${true_seqs} ${null_seqs} ${outpre}"
        else
            echo "File ${resultfile} already exists! Skipping..."
        fi
    done

    # Next fit 'full models' that are trained on all training data
    
    # Get true and null seq fasta files for training
    true_seq_files=( $(printf "${fasta_dir}/${model}_%s_true_seqs.fasta " "${all_chrs[@]}") )
    null_seq_files=( $(printf "${fasta_dir}/${model}_%s_null_seqs.fasta " "${all_chrs[@]}") )

    outpre=${full_model_dir}/${model}.full
    resultfile=${outpre}.model.txt

    if [ ! -f ${resultfile} ]
    then
        # Create files for combined training data
        echo "Combining training fastas..."
        true_seqs=${train_dir}/${model}.full.true_seqs.fasta
        null_seqs=${train_dir}/${model}.full.null_seqs.fasta
        cat "${true_seq_files[@]}" > ${true_seqs}
        cat "${null_seq_files[@]}" > ${null_seqs}

        # Train models
        echo "Training model for ${model}..."
        sbatch -p wjg,biochem,sfgf -t 72:00:00 --mem=40G --cpus-per-task=16 \
        --job-name=${job_header} --output=${outpre}.log --error=${outpre}.log \
        --wrap "gkmtrain -y ${yval} -t ${tval} -l ${lval} -k ${kval} -d ${dval} -g ${gval} -c ${cval} -M ${Mval} -H ${Hval} -e ${eval} \
        -m ${mval} -r 123 -T 16 ${true_seqs} ${null_seqs} ${outpre}"
    else
        echo "File ${resultfile} already exists! Skipping..."
    fi

done
