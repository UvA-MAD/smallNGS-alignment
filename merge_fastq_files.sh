#!/bin/bash

#### Merge fastq files from different chips
#### 2014-05-22 Inez Terpstra

#----------------------------------------------------------------------------------------------------

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

echo -n "First chip (RIDxxxx): "
read chip1
echo -n "Second chip (RIDxxxx): "
read chip2
echo -n "Which files to process? (a(ll) or pattern: eg 001 or 005*RID0035): "
read selection

if [ ${selection:0:1} = a ] 
then 
  list_of_files=*.fastq
else
  list_of_files=IonXpressRNA_$selection*fastq
fi

#echo "The experiment environment is: $base_dir"
echo client_type is $client_type
echo client nr is $client_nr
echo project nr is $project_nr
echo experiment nr is $experiment_nr
echo experiment name is ${experiment##$client_nr-$project_nr-$experiment_nr_}
echo The run environment is $exp
echo The files to process are ${list_of_files##"$fq_dir/"}
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

#---------------------------------------------------------------------------------------------------
## merge files
#---------------------------------------------------------------------------------------------------
## process files
#---------------------------------------------------------------------------------------------------

for filename in $fq_dir/$chip1/$list_of_files ; do
  echo $filename 
  barcode=${filename##*RNA_}
  barcode=${barcode%%_R*}
  echo $barcode
  cat $filename > "IonXpressRNA_"$barcode"_merged-"$chip1.fastq
done

for filename in $fq_dir/$chip2/$list_of_files ; do
  echo $filename 
  barcode=${filename##*RNA_}
  barcode=${barcode%%_R*}
  echo $barcode
  cat $filename >> "IonXpressRNA_"$barcode"_merged-"$chip1.fastq
done

exit

