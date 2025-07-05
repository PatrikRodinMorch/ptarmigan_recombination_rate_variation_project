#!/bin/bash -l

vcf=/crex/proj/snic2020-6-18/mapping/picard_out/duplMarked_bam/Haplotypecaller/gVCF_Patrik/mapped_to_RP_ref/merged_vcf/merged_genotyped/filtered_vcf/indel_clean/analysis_ready_vcf/analysis/ReLERNN/ptarm_data/${1}_biallelic.vcf

directory=/crex/proj/snic2020-6-18/mapping/picard_out/duplMarked_bam/Haplotypecaller/gVCF_Patrik/mapped_to_RP_ref/merged_vcf/merged_genotyped/filtered_vcf/indel_clean/analysis_ready_vcf/analysis/ReLERNN/ptarm_data

conda activate ReLERNN

./ReLERNN/ReLERNN_PREDICT --vcf $vcf \
--projectDir $directory/${1} \
--seed 42 \
--unphased

conda deactivate
