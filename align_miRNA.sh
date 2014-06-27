#!/bin/bash
#### Alignment to miRNAs
#### 2014-05-12 Inez Terpstra
#### Reads are trimmed with Trimmomatic0.30
#### Trimmed reads are aligned to pre-miRNA with SHRiMP2
#### - database projected with seeds: 
####   111011110111,111011011111,11111110111,111101110111
####   seeds can be a bit longer for dre alignments:
####   1101110110111,11101100101111,11101111010011,1110110110111
#### - local and global alignment
#### - alignment to Forward strand only
#### conversion of sam output to bam with samtools
#### Alignments are filtered in R
#### - 10% mismatches are allowed
#### - each InDel counts as 1 mismatch
#### - a 3p overhang of 3 nucleotides is allowed
#### - aligned reads are assigned to mature miRNA if the read extends the mature sequence with 3 or less nucleotides
#----------------------------------------------------------------------------------------------------

base_dir=$1
exp=$2
filename=$3
chip=$4
barcode=$5
org=$6
db_dir=$7
db_in=$8
crop=$9
minlen=${10}
seqnr=${11}
full_gap_penalty=${12}

echo "Alignment to miRNAs"

#### set environment
src_dir=/mad/MAD-RBAB/10_Protocols/DryLab/smallNGS_alignment
map_dir=$base_dir/Scratch/$exp/Data/miRNA_aligned_data
trim_dir=$base_dir/Scratch/$exp/Data/trimmed_data
log_dir=$base_dir/Scratch/$exp/Data/log_data
export SHRIMP_FOLDER=/mad/software/src/SHRiMP_2_2_3

fbase=$chip"_"$barcode"_trm"$seqnr
echo "File base name is $fbase"

#### trim data
java -jar /mad/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar SE -threads 8 -phred33 -trimlog $log_dir/$fbase.log $filename $trim_dir/$fbase.fastq CROP:$crop MINLEN:$minlen

#### alignment
f_map=$fbase"_"$db_in

cd $db_dir
## check if database exists
if [ ! -f "$db_in-ls.genome" ]; then
  echo "Creating database projection"
  if [ $org = dre ] ; then $SHRIMP_FOLDER/utils/project-db.py --h-flag --seed 1101110110111,11101100101111,11101111010011,1110110110111 --shrimp-mode ls "$db_dir/$db_in.fa"
  else 
    $SHRIMP_FOLDER/utils/project-db.py --h-flag --seed 111011110111,111011011111,11111110111,111101110111 --shrimp-mode ls "$db_dir/$db_in.fa"
fi

#align to forward strand only
$SHRIMP_FOLDER/bin/gmapper-ls -L $db_dir/$db_in"-ls" $trim_dir/$fbase.fastq -F -H -w 100% -g -10 -q -10 -e -1 -f -1 -n 1 -o 1 --strata -N 8 --qv-offset 33 --single-best-mapping >$map_dir/$f_map"_"global.sa
$SHRIMP_FOLDER/bin/gmapper-ls -L $db_dir/$db_in"-ls" $trim_dir/$fbase.fastq -F -H -w 100% -g -10 -q -10 -e -1 -f -1 -n 1 -o 1 --strata --local -N 8 --qv-offset 33 --single-best-mapping >$map_dir/$f_map"_"local.sa

cd $map_dir
samtools view -bS -o $f_map"_"global.bam $f_map"_"global.sa
samtools sort $f_map"_"global.bam sorted_$f_map"_"global
samtools index sorted_$f_map"_"global.bam
samtools view -bS -o $f_map"_"local.bam $f_map"_"local.sa
samtools sort $f_map"_"local.bam sorted_$f_map"_"local
samtools index sorted_$f_map"_"local.bam

#### process miRNA alignments in R
$src_dir/process_miRNAs.R $base_dir $exp $filename $chip $barcode $org $f_map $db_dir $db_in $full_gap_penalty


