#!/bin/bash

#### Alignment analysis per Barcode
#### 2014-05-12 Inez Terpstra
#### all scripts collected in one procedure to align smallRNA sequencing data to spikes, miRNA, piRNA, rRNA
#### spikes - version 1 with RNA size (8,10,16,19,22,25,28,34,40,50,60,70 nt) and 19 normalization spikes
#### miRNA - miRBase version 20  (http://www.mirbase.org/)
#### piRNA - piRNABank april 2014 (http://pirnabank.ibab.ac.in)
#### environment: genseq hpc cluster MAD-RBAB
#### Start in fastq directory (Scratch) and get project environment from current directory
#### Start with sff files in future

#----------------------------------------------------------------------------------------------------

src_dir=/mad/MAD-RBAB/10_Protocols/DryLab/smallNGS_alignment

fq_dir=$PWD

base_dir=${fq_dir%%/Scratch*}
client_type=${fq_dir##*/MAD-RBAB/}
client_type=${client_type%%/MAD[0-9]*}
experiment=${base_dir##*/}
client_nr=${experiment:0:7}
project_nr=${experiment##$client_nr-}
project_nr=${project_nr:0:4}
experiment_nr=${experiment##$client_nr-$project_nr-}
experiment_nr=${experiment_nr:0:9}

echo -n "What is the run environment? (bioinf experiment name): "
read exp
echo -n "Which files to process? (a(ll) or pattern: eg 001 or 005*RID0035): "
read selection
echo -n "Select organism (ath,cel,dme,dre,hsa,mmu,rno): "
read org
echo -n "Align to spikes? y/n: "
read spikes
echo -n "Align to miRNA? y/n: "
read mir
echo -n "Align to piRNA? y/n: "
read pir
echo -n "Align to rRNA? y/n: "
read rrna
echo -n "Align to tRNA? y/n: "
read trna
echo -n "Normalize data? y/n: "
read norm

if [ ${selection:0:1} = a ] 
then 
  list_of_files=*.fastq
else
  list_of_files=IonXpressRNA_$selection*fastq
fi

echo client_type is $client_type
echo client nr is $client_nr
echo project nr is $project_nr
echo experiment nr is $experiment_nr
echo experiment name is ${experiment##$client_nr-$project_nr-$experiment_nr_}
echo The run environment is $exp
echo The files to process are ${list_of_files##"$fq_dir/"}
echo The selected organism is $org
echo Alignment to:
if [ $spikes = "y" ]; then echo "spikes"; fi
if [ $mir = "y" ]; then echo "miRNA"; fi
if [ $pir = "y" ]; then echo "piRNA"; fi
if [ $rrna = "y" ]; then echo "rRNA"; fi
if [ $trna = "y" ]; then echo "tRNA"; fi
if [ $norm = "y" ]; then echo "normalize data"; fi
echo -n "Proceed? y/n: "
read proceed
case $proceed in
  "n") 
    exit 1
  ;;
  "") 
    exit 1
  ;;
esac

#### create environment if necessary
$src_dir/create_experiment_directory.sh $base_dir $exp 

## reference databases
db_dir_base=/mad/MAD-RBAB/05_Reference-db

db_dir_spikes=$db_dir_base/RBAB/spikes
db_spikes=spikes40
crop_spikes=40

db_dir_mir=$db_dir_base/external/$org/RNA/miRNA
db_mir=hairpin20
case $org in
  "ath")
    minlen_mir=15
    crop_mir=30
  ;;
  "cel")
    minlen_mir=15
    crop_mir=40
  ;;
  "dme")
    minlen_mir=15
    crop_mir=40
  ;;
  "dre")
    minlen_mir=15
    crop_mir=30
  ;;
  "hsa")
    minlen_mir=12
    crop_mir=40
  ;;
  "mmu")
    minlen_mir=12
    crop_mir=40
  ;;
  "rno")
    minlen_mir=12
    crop_mir=40
  ;;
esac

db_dir_pir=$db_dir_base/external/$org/RNA/piRNA
db_pir=piRNA
minlen_pir=12
crop_pir=40

db_dir_rrna=$db_dir_base/external/$org/RNA/rRNA
case $org in
  "ath")
    db_rrna=unknown
  ;;
  "cel")
    db_rrna=unknown
  ;;
  "dme")
    db_rrna=unknown
  ;;
  "dre")
    db_rrna=rRNA
  ;;
  "hsa")
    db_rrna=rRNA
  ;;
  "mmu")
    db_rrna=unknown
  ;;
  "rno")
    db_rrna=unknown
  ;;
esac
minlen_rrna=12
crop_rrna=40

db_dir_trna=$db_dir_base/external/$org/RNA/tRNA
case $org in
  "ath")
    db_rrna=unknown
  ;;
  "cel")
    db_rrna=unknown
  ;;
  "dme")
    db_rrna=unknown
  ;;
  "dre")
    db_rrna=rRNA
  ;;
  "hsa")
    db_rrna=rRNA
  ;;
  "mmu")
    db_rrna=unknown
  ;;
  "rno")
    db_rrna=unknown
  ;;
esac
minlen_trna=12
crop_trna=40

#---------------------------------------------------------------------------------------------------
## perform alignments
#---------------------------------------------------------------------------------------------------
## process files
#---------------------------------------------------------------------------------------------------

for filename in $list_of_files ; do

  chip=${filename%%.fastq}
  chip=${chip##*-}
  echo $chip
  barcode=${filename##*RNA_}
  barcode=${barcode%%_*}
  echo $barcode

  new_filename=$base_dir/Scratch/Fastq/$filename

  ### spikes
  if [ $spikes = "y" ] 
  then
    echo ${new_filename##$base_dir/Scratch/}
    $src_dir/align_spikes.sh $base_dir $exp $filename $chip $barcode $db_dir_spikes $db_spikes $crop_spikes "1"
    new_filename=$base_dir/Scratch/$exp/Fastq/unmapped_spikes/$chip"_"$barcode"_unm.fastq"
  fi

  ### miRNA
  if [ $mir = "y" ] 
  then
    ## check if database is present
    if [ ! -f "$db_dir_mir/$db_mir.fa" ]; then
    echo "miRNA data is not present for this organism; skipping miRNA alignment"
    else 
      echo ${new_filename##$base_dir/Scratch/}
      full_gap_penalty="FALSE"
      $src_dir/align_miRNA.sh $base_dir $exp $new_filename $chip $barcode $org  $db_dir_mir $db_mir $crop_mir $minlen_mir "2" $full_gap_penalty
      new_filename=$base_dir/Scratch/$exp/Fastq/unmapped_miRNA/$chip"_"$barcode"_unm.fastq"
    fi
  fi

  ### piRNA
  if [ $pir =  "y" ]
  then
    ## check if database is present
    if [ ! -f "$db_dir_pir/$db_pir.fa" ]; then
    echo "piRNA data is not present for this organism; skipping piRNA alignment"
    else 
      echo ${new_filename##$base_dir/Scratch/}
      full_gap_penalty="TRUE"
      $src_dir/align_piRNA.sh $base_dir $exp $new_filename $chip $barcode $org  $db_dir_pir $db_pir $crop_pir $minlen_pir "3" $full_gap_penalty
      new_filename=$base_dir/Scratch/$exp/Fastq/unmapped_piRNA/$chip"_"$barcode"_unm.fastq"
    fi
  fi

  ### rRNA
  if [ $rrna =  "y" ]
  then
    ## check if database is present
    if [ ! -f "$db_dir_rrna/$db_rrna.fa" ]; then
    echo "rRNA data is not present for this organism; skipping rRNA alignment"
    else 
      echo ${new_filename##$base_dir/Scratch/}
      full_gap_penalty="TRUE"
      $src_dir/align_smallRNA.sh $base_dir $exp $new_filename $chip $barcode $org "rRNA" $db_dir_rrna $db_rrna $crop_rrna $minlen_rrna "4" $full_gap_penalty
      new_filename=$base_dir/Scratch/$exp/Fastq/unmapped_rRNA/$chip"_"$barcode"_unm.fastq"
    fi
  fi

  ### tRNA
  if [ $trna =  "y" ]
  then
    ## check if database is present
    if [ ! -f "$db_dir_trna/$db_trna.fa" ]; then
    echo "tRNA data is not present for this organism; skipping tRNA alignment"
    else 
      echo ${new_filename##$base_dir/Scratch/}
      full_gap_penalty="TRUE"
      $src_dir/align_smallRNA.sh $base_dir $exp $new_filename $chip $barcode $org "tRNA" $db_dir_trna $db_trna $crop_trna $minlen_trna "5" $full_gap_penalty
      new_filename=$base_dir/Scratch/$exp/Fastq/unmapped_tRNA/$chip"_"$barcode"_unm.fastq"
    fi
  fi

done

#---------------------------------------------------------------------------------------------------
## collect results of all aligned reads
#---------------------------------------------------------------------------------------------------
if [ $spikes = "y" ]; then $src_dir/collect_results.R $base_dir $exp "spike"; fi
if [ $mir = "y" ]
then
  $src_dir/collect_results.R $base_dir $exp "miRNA"
  $src_dir/collect_results_miRNA.R $base_dir $exp
fi
if [ $pir = "y" ]; then $src_dir/collect_results.R $base_dir $exp "piRNA"; fi
if [ $rrna = "y" ]; then $src_dir/collect_results.R $base_dir $exp "rRNA"; fi
if [ $trna = "y" ]; then $src_dir/collect_results.R $base_dir $exp "tRNA"; fi

#---------------------------------------------------------------------------------------------------
## normalizes results of all aligned reads
#---------------------------------------------------------------------------------------------------
if [ $norm = "y" ]; 
then 
  $src_dir/normalize_data.R $base_dir $exp $org; fi

#---------------------------------------------------------------------------------------------------
## plot results of all aligments
#---------------------------------------------------------------------------------------------------
$src_dir/plot_results_general.R $base_dir $exp 
if [ $spikes = "y" ]; then $src_dir/plot_results_spikes.R $base_dir $exp $db_dir_spikes $db_spikes; fi
if [ $mir = "y" ]; then $src_dir/plot_results.R $base_dir $exp "miRNA" $spikes; fi
if [ $pir = "y" ]; then $src_dir/plot_results.R $base_dir $exp "piRNA" $spikes; fi
if [ $rrna = "y" ]; then $src_dir/plot_results.R $base_dir $exp "rRNA" $spikes; fi
if [ $trna = "y" ]; then $src_dir/plot_results.R $base_dir $exp "tRNA" $spikes; fi

exit


