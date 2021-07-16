#!/usr/bin/bash

#####	02/01/2019	#####
#####	Downloading the data from PacBio tar files (from Shana) without directory structure	#####

set -e 

if [ "$#" -ne 1 ]; then
    printf "Illegal number of parameters:\t	$#	\n";
    printf "\tPlease RUN as : sh extract_pacbiorun_tar.sh FilePattern\n";
    printf "\t\tExample:	sh extract_pacbiorun_tar.sh /PATH/TO/FILES/\n";				#sh extract_pacbiorun_tar.sh $PWD
    exit;
fi

TAR_PATHING=$1"/*.tar"
echo $TAR_PATHING

checking=$(ls $TAR_PATHING | wc -l)
if [ "$checking" -eq 0 ]; then
    printf ".tar files with pattern $checking dont exist\n\tPLEASE CHECK for $1*.tar\n";
    exit
fi

for file1 in $TAR_PATHING; do
  file2=$(echo $file1 | sed -e 's/\.tar$//'); 
  printf "Starting $file1\n";
  mkdir -p $file2;

  stripping=$(less $file1 | cat | awk '{print $NF}' | sed -e 's/\//\t/g' | awk -F'\t' '{print NF-1}' | sort -un)
  depth=$(less $file1 | cat | awk '{print $NF}' | sed -e 's/\//\t/g' | awk -F'\t' '{print NF-1}' | sort -un | wc -l)

  if [ "$depth" -eq "0" ]; then
    printf "Something wrong with folder depth; Check if there are other subfolders with files in the following output\n";
    less $file1 | cat | awk '{print $NF}';
    exit;
  fi

  printf "\tRunning:\t tar -xvf $file1 -C $file2/ --strip=$stripping\n";
  tar -xvf $file1 -C $file2/ --strip=$stripping

###  tar -xf $file1 -C $file2/ --transform='s/.*\///';				### Transforming Without stripping
###  printf "\tRemoving excess directories for $file1\n";
###  find $file2 -type d -empty -delete 					### To Delete empty directory, but we shouldn't have subdirectories

  printf "\tProcessed   $file1 \n"; 
done	


