#!/bin/bash
# Okt 2019

echo "---------- CHANGING FOLDER NAME FROM BEARCAVE-master TO BEARCAVE -----------"

cd ..
mv BEARCAVE-master/ BEARCAVE/
cd BEARCAVE/

echo "---------- CREATING SUBFOLDER STRUCTURE -----------"


mkdir rawdata/old_metadata
mkdir software
mkdir trimdata
mkdir trimdata/trimlogs
mkdir refgenomes

softdir=./software/miniconda3/


if [[ "$OSTYPE" == "linux"* ]]; then

  echo "---------- DOWNLOADING MINICONDA3 4.5.11 for Linux -----------"
  wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh
  chmod u+x Miniconda3-4.5.11-Linux-x86_64.sh

  echo "---------- INSTALLING MINICONDA3 for Linux -----------"
  ./Miniconda3-4.5.11-Linux-x86_64.sh -f -p $PWD/software/miniconda3;

elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo "---------- DOWNLOADING MINICONDA3 4.5.11 for macOS -----------"
  curl https://repo.continuum.io/miniconda/Miniconda3-4.5.11-MacOSX-x86_64.sh -o Miniconda3-4.5.11-MacOSX-x86_64.sh
  chmod u+x Miniconda3-4.5.11-MacOSX-x86_64.sh

  echo "---------- INSTALLING MINICONDA3 for macOS -----------"
  ./Miniconda3-4.5.11-MacOSX-x86_64.sh -f -p $softdir;
else
  echo "---------- Operatingsystem not supported -----------"
  exit
fi

echo "------Setting up Python 3.5.6 enviroment------"
$softdir/bin/conda create -f --name py35 python=3.5.6
sleep 2
echo "---------- INSTALLING CUTADAPT 1.12 -----------"
$softdir/bin/conda install -y -f -n py35 -c bioconda cutadapt=1.12 
sleep 2
echo "---------- INSTALLING FLASH 1.2.11 -----------"
$softdir/bin/conda install -y -f -n py35 -c bioconda flash=1.2.11
sleep 2
echo "---------- INSTALLING BWA 0.7.15 -----------"
$softdir/bin/conda install -y -f -n py35 -c bioconda bwa=0.7.15
sleep 2
echo "---------- INSTALLING SAMTOOLS 1.3.1 -----------"
$softdir/bin/conda install -y -f -n py35 -c bioconda samtools=1.3.1


echo "---------- SETTING PERMISSIONS -----------"
chmod 770 .

chmod 440 install.sh

chmod -R a-w software/
chmod -R o-rx software/
chmod -R 550 scripts/

chmod 770 rawdata/
chmod 440 rawdata/metadata.txt
chmod 770 rawdata/old_metadata/


chmod -R 770 trimdata/

chmod 770 refgenomes/

echo "---------- BEARCAVE INSTALLATION COMPLETE -----------"
echo -e '\n'
echo "You may want to update your working directory by navigating one folder up and back into BEARCAVE/"
echo -e '\n'
