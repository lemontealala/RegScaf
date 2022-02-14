##!/bin/bash

script_path=/home/limengtian/RegCodes_final/scripts
output_dir=$2
genome=$1

if [ -d "$output_dir" ]; then
        mv $output_dir ${output_dir}_old;
fi
mkdir $output_dir
cp $1 $output_dir/genome.fasta

##Prepare contigs:
if [ -d "$output_dir/CONTIG" ]; then
	rm -rf $output_dir/CONTIG;
fi
mkdir $output_dir/CONTIG
cd $output_dir
#$script_path/CutContigAtN genome.fasta ./CONTIG/ ./detail.txt ./evidence.txt.cutcontigatn
python $script_path/CutScaffoldAtN_assess.py
ln -s initial_contig.fa reference.fasta

##Mapping
cp ../config/Mapping.sh .
bash ./Mapping.sh > Mapping.log 2>&1

##Filter high coverage region
if [ -e "automap.bed" ]; then
	mv automap.bed automap.bed.old
fi

if [ -e "rautomap.bed" ]; then
        mv rautomap.bed rautomap.bed.old
fi

if [ -e "highdepth.bed" ]; then
        mv highdepth.bed highdepth.bed.old
fi

cp ../config/Bed.sh ../config/name_list.txt .
bash Bed.sh name_list.txt 6 > Bed.log 2>&1

##Scaffolding by RegScaf
if [ -d "RegScaf" ]; then
	mv RegScaf RegScaf_old
fi
mkdir RegScaf

cp ../config/MakeTab.sh .;
bash ./MakeTab.sh >MakeTab.log 2>&1
cp ../config/RegLib.txt ../config/ScafAssess.sh RegScaf/;
cd RegScaf/;cp ../detail.txt .; bash ./ScafAssess.sh > ScafAssess.log 2>&1
