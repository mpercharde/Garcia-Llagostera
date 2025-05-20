#!/bin/bash

# load modules
module load openjdk

# run ChromHMM OverlapEnrichment for all states over genome features
java -mx8000M -jar /home/programs/ChromHMM/ChromHMM.jar OverlapEnrichment segments.bed  \
  /home/mnt/network/DATA/projects/isgs/analysis/chromHMM/ptm_tar/oe_analysis/regions/genome_features/ model_14-oe_genome_features

# run ChromHMM OverlapEnrichment for all states over gene sets
java -mx8000M -jar /home/programs/ChromHMM/ChromHMM.jar OverlapEnrichment -uniformscale segments.bed \
/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/ptm_tar/oe_analysis/regions/gene_sets/ model_14-oe_gene_sets

echo "done"
