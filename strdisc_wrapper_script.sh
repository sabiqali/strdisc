#!/bin/sh

#dependencies
strdisc_conda_dir = "path/to/conda/dir"

#required parameters
indel_file_location = "path/to/indel/file"
bam_file_location = "path/to/bam/file"
output_bed_file = "output/file/name/with/path"

#optional parameters
extra_options = "--upper_length 3" #please refer to strdisc help menu for extra options

#using the SV analysis pipeline to generate the required bed files. 
snakemake --jobs 500 --rerun-incomplete -s ./nanopore-SV-analysis/Snakefile --keep-going --latency-wait 120 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -q all.q -P simpsonlab -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_rt={params.run_time} -b y"

#use the bed file generated with strdisc to generate bed file with up to top 3 candidates from each strand
$strdic_conda_dir/python ./scripts/strdisc.py --indel_file $indel_file_location --bam $bam_file_location --bed $output_bed_file $extra_options 