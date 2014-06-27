#!/bin/bash
#### Create experiment directories if they don't exist
#### 2014-05-12 Inez Terpstra
#----------------------------------------------------------------------------------------------------
base_dir=$1
exp=$2

if [ ! -d "$base_dir/Results/$exp" ]
then
  echo -n "Directory $base_dir/Results/$exp does not exist, should it? y/n: "
  read proceed
  case $proceed in
    "n") 
      exit 1
    ;;
    "") 
      exit 1
    ;;
  esac
  mkdir $base_dir/Images/$exp 
  mkdir $base_dir/Images/$exp/Box 
  mkdir $base_dir/Images/$exp/MA 
  mkdir $base_dir/Images/$exp/PCA 
  if [ ! -d "$base_dir/Scripts/$exp" ]
  then
    mkdir $base_dir/Scripts/$exp 
  fi
  mkdir $base_dir/Scratch/$exp 
  mkdir $base_dir/Scratch/$exp/Fastq
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_spikes
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_miRNA
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_mature
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_piRNA
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_tRNA
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_rRNA
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_cDNA
  mkdir $base_dir/Scratch/$exp/Fastq/unmapped_genome
  mkdir $base_dir/Scratch/$exp/Data 
  mkdir $base_dir/Scratch/$exp/Data/spike_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/miRNA_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/mature_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/piRNA_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/tRNA_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/rRNA_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/cDNA_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/genome_aligned_data
  mkdir $base_dir/Scratch/$exp/Data/log_data
  mkdir $base_dir/Scratch/$exp/Data/trimmed_data
  mkdir $base_dir/Results/$exp
  mkdir $base_dir/Results/$exp/Output
  mkdir $base_dir/Results/$exp/DESeq
  mkdir $base_dir/Results/$exp/DESeq/Info
  mkdir $base_dir/Results/$exp/DESeq/Results
fi

