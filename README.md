# scScalpChromatin

Repository to host code for ["Integrated single-cell chromatin and transcriptomic analyses of human scalp identify gene-regulatory programs and critical cell types for hair and skin diseases" by Ober-Reynolds et al. Nature Genetics. 2023.](https://doi.org/10.1038/s41588-023-01445-4)

Scripts for performing initial scATAC data processing and subclustering are found in scATAC directory. Scripts for performing initial scRNA data processing and subclustering are found in the scRNA directory. Numbers preceeding these scripts indicate the order in which scripts should be used.

Fine-mapped SNP analyses (LDSR, intersection with peak-to-gene linkages, and GkmSVM models) are found in the GWAS directory.

Processed and annotated Seurat objects for scRNA-seq data can be found [here](https://drive.google.com/drive/folders/1klScH010LvYxdU-TZiSkGVip9IrilZN2?usp=sharing). This includes the objects for the subclustered keratinocytes, fibroblasts, endothelial, lymphoid, and myeloid groups. (We may change the location of these files in the future if we find a better place to host them, but will update this page with the location.)
