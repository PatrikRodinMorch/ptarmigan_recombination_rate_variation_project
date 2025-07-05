#!/bin/bash -l


conda activate ReLERNN

./ReLERNN/ReLERNN_SIMULATE --vcf ptarm_data/${1}_biallelic.vcf \
--genome ptarm_data/bLagMut_sorted_autosomes_final.bed \
--mask ptarm_data/unmapped_regions_unfolded_new.bed \
--projectDir ptarm_data/${1} \
--demographicHistory ptarm_data/${1}.csv \
--assumedMu 0.3e-8 \
--assumedGenTime 2 \
--seed 42 \
--unphased

conda deactivate
