#!/bin/bash

###########################################################
# Download required data and prepare peakset for ldsc
###########################################################

# Stop execution if we get an error along the way
set -o errexit 

# Source config file
source ${HOME}/git_clones/scScalpChromatin/GWAS/LDSC/00_ldsc_config.conf

# Make results directory if it doesn't already exist
if [ ! -d ${wd} ]; then mkdir ${wd}; fi

# Download baseline LD scores:
if [ ! -d "${baseline_dir}" ]
then
    echo "Downloading baseline ldscores..."
    wget ${baselinesource} -P ${wd}
    tar -xvzf ${wd}/${baselinefile}
    rm ${wd}/${baselinefile}
    if [ ! -d "${baseline_dir}" ]
    then
        mkdir ${baseline_dir}
        mv baselineLD.* ${baseline_dir}
    fi
else
    echo "baseline ldscores already present..."
fi

# Get HapMap3 weights (without HLA loci)
if [ ! -d "${weights_dir}" ]
then
    echo "Downloading HapMap3 weights..."
    wget ${weightssource} -P ${wd}
    tar -xvzf ${wd}/${weightsfile}
    rm ${wd}/${weightsfile}
else
    echo "HapMap3 weights already present..."
fi

# Download 1000G Phase3 plink files:
if [ ! -d "${plink_dir}" ]
then
    echo "Downloading 1000G Phase3 plink files..."
    wget ${plinksource} -P ${wd}
    tar -xvzf ${wd}/${plinkfile}
    rm ${wd}/${plinkfile}
else
    echo "1000G Phase3 plink files already present..."
fi

# Download HapMap3 SNPs:
if [ ! -d "${snp_dir}" ]
then
    echo "Downloading HapMap3 SNPs..."
    wget ${snpsource} -P ${wd}
    mkdir ${snp_dir}
    mv ${snpfile} ${snp_dir}
    # tar -xvzf ${wd}/${snpfile}
    # rm ${wd}/${snpfile}
else
    echo "HapMap3 SNPs already present..."
fi

# Download frequency files, if needed:
if [ ! -z "${frqfile}" ]
then
    if [ ! -d "${frq_dir}" ]
    then
        echo "Downloading frequency files..."
        wget ${frqsource} -P ${wd}
        tar -xvzf ${wd}/${frqfile}
        rm ${wd}/${frqfile}
    else
        echo "Frq files already present..."
    fi
fi


############################################
# Identify cell type-specific peaks:
############################################

# Load required modules
ml R/4.0.2 system 

if [ ! -d "${cts_peak_dir}" ]
then
    echo "Identifying cell type-specific peaks..."
    mkdir ${cts_peak_dir}
    logfile=${cts_peak_dir}/get_cell_type_peaks.out
    sbatch -p wjg,biochem,sfgf -t 00:20:00 --mem=50G --cpus-per-task=8 \
    --job-name=ctspeaks --output=${logfile} --error=${logfile} \
    --wrap "Rscript ${HOME}/git_clones/scScalpChromatin/GWAS/LDSC/get_cell_type_peaks.R ${peakClass} ${archrPath} ${wd}"
else
    echo "cell type-specific peaks already present..."
fi

############################################