#!/bin/sh

#make sure to edit the SV pipeline config and then edit this

#dependencies
strdisc_conda_dir="path/to/conda/envs/strdisc/bin"
strdisc_dir="path/to/strdisc"
output_dir="path/to/output/dir" #output directory where you want to store SV pipeline output along with the strdisc output
sample_names="sample_name" #enter same sample name here as the config for the SV calling pipeline

#required parameters
bam_file_location="path/to/bam/file"
output_bed_file="output/bed/file/name"

#indel file locations
indel_file_locations="$output_dir/$sample_names/structural_variants/$sample_names.insertions.bed $output_dir/$sample_names/structural_variants/$sample_names.duplications.bed"

#strdisc python from conda directory
strdisc_python="$strdisc_conda_dir/python"

#snakefile location
snakefile_loc="$strdisc_dir/nanopore-SV-analysis/Snakefile"

#strdisc executable
strdisc_exec="$strdisc_dir/scripts/strdisc.py"

#change to output directory and create subdirectories for the samples
mkdir $output_dir
cp ./config.taml $output_dir
cd $output_dir
for sample_name in $sample_names; do
    mkdir $sample_name
done

#optional parameters
extra_options="--upper_length 3" #please refer to strdisc help menu for extra options

#using the SV analysis pipeline to generate the required bed files. 
snakemake --jobs 500 --rerun-incomplete -s $snakefile_loc --keep-going --latency-wait 120 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -q all.q -P simpsonlab -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_rt={params.run_time} -b y"

#use the bed file generated with strdisc to generate bed file with up to top 3 candidates from each strand
for indel_file_location in $indel_file_locations; do
    qsub -cwd -V -N strdisc -l h_vmem=32G -l h_stack=32M -l h_rt=5:0:0:0 -P simpsonlab -b y "$strdisc_python $strdisc_exec --indel_file $indel_file_location --bam $bam_file_location --bed $output_bed_file $extra_options"
    wait 
done
