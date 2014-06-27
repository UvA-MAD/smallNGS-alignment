#!/bin/bash
#### Alignment to spikes
#### 2014-05-12 Inez Terpstra
#### Spikes and reads are trimmed with Trimmomatic0.30
#### Trimmed reads are aligned to spikes with Bowtie2
#### - alignment to Forward strand only
#### - 10% mismatches allowed
#### conversion of sam output to bam with samtools
#### Alignments are filtered in R
#### - aligned reads should be 80% of spike length
#----------------------------------------------------------------------------------------------------

base_dir=$1
exp=$2
filename=$3
chip=$4
barcode=$5
db_dir=$6
db_in=$7
crop=$8
seqnr=$9

echo "Alignment to spikes"

#### set environment
src_dir=/mad/MAD-RBAB/10_Protocols/DryLab/smallNGS_alignment
map_dir=$base_dir/Scratch/$exp/Data/spike_aligned_data
trim_dir=$base_dir/Scratch/$exp/Data/trimmed_data
log_dir=$base_dir/Scratch/$exp/Data/log_data

fbase=$chip"_"$barcode"_trm"$seqnr

echo "File base name is $fbase"

#### trim data
java -jar /mad/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar SE -threads 8 -phred33 -trimlog $log_dir/$fbase.log $filename $trim_dir/$fbase.fastq CROP:$crop

#### alignment
f_map=$fbase"_"$db_in

cd $db_dir
## check if database exists
if [ ! -f "$db_in.1.bt2" ]; then
  echo "Creating database index"
  /mad/bin/bowtie2-build -f $db_in".fa" $db_in
fi

/mad/bin/bowtie2 -L 6 -i S,0,0.5 --ignore-quals --norc --score-min L,-1,-0.6 -D 20 -t -p 8 -q -x $db_in -U $trim_dir/$fbase.fastq -S $map_dir/$f_map.sa
    
cd $map_dir
samtools view -bS -o $f_map.bam $f_map.sa
samtools sort $f_map.bam sorted_$f_map
samtools index sorted_$f_map.bam

#### process spike alignments in R
$src_dir/process_spikes.R $base_dir $exp $filename $chip $barcode $f_map $db_dir $db_in



