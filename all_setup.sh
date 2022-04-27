#download and install the packages and dependencies needed by the SV analysis pipeline
conda env create -n nanopore-SV-analysis -f ./nanopore-SV-analysis/nanopore-SV-analysis.yml

#download and install the packages and dependencies needed by strdisc
conda env create -n strdisc -f strdisc.yml
conda activate strdisc

#get references needed by SV analysis pipeline
chmod u+x ./download_reference.sh
./download_reference.sh
