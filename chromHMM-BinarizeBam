#!/bin/bash
#SBATCH --job-name=chromHMM_BinarizeBam.mm39
#SBATCH --output=/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/reports/chromHMM_BinarizeBam.mm39_%A.out
#SBATCH --time=48:00:00
#SBATCH --partition=hmem
#SBATCH --cpus-per-task=20
#SBATCH --qos=qos_batch

######################################################
# README

# this script runs ChromHMM BinarizeBed for files in a specified directory, sent to a specified queue
# sbatch /path/to/script/

######################################################
# MODULES + PATHS + REPORTING

## Load modules
module load openjdk

echo "~ChromHMM parameters and info~"
echo "genome: GRCm39/mm39, annotation: Ensembl v110 (gtf)"
echo "running ChromHMM BinarizeBam"
echo "binsize is 200 (default) and shiftsize is 100 (default)"

## paths and information
inBam_tar="/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/in_bam/targets/"
inBam_control="/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/in_bam/input/"
out_tar="/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/ptm_tar/binarized_bam/targets"
out_control="/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/ptm_tar/binarized_bam/input"
out_signal="/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/ptm_tar/signal_out"

info_tab="/home/mnt/network/DATA/projects/isgs/analysis/chromHMM/sample_info"
chrom_sz="/home/genomes/mm39/GRCm39_primary_assembly.chrom_sizes"

######################################################
# RUN

java --version

echo "Running BinarizeBam"
java -mx4000M -jar /home/programs/ChromHMM/ChromHMM.jar BinarizeBam -b 200 -c ${inBam_control} -n 100 -o ${out_control} -t ${out_signal} -mixed ${chrom_sz} ${inBam_tar} $info_tab/cellmarkfiletable.txt ${out_tar}

echo "Run completed"


