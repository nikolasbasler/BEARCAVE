#!/bin/bash
# December 2018


mkdir rawdata/old_metadata
mkdir software
mkdir trimdata
mkdir trimdata/trimlogs
mkdir refgenomes

cd software

echo "---------- DOWNLOADING MINICONDA3 4.5.11 -----------"
wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh
chmod u+x Miniconda3-4.5.11-Linux-x86_64.sh

echo "---------- INSTALLING MINICONDA3 -----------"
./Miniconda3-4.5.11-Linux-x86_64.sh -f -p $PWD/miniconda3


echo "---------- INSTALLING CUTADAPT 1.12 -----------"
./miniconda3/bin/conda install -y -f -c bioconda cutadapt=1.12


echo "---------- INSTALLING FLASH 1.2.11 -----------"
./miniconda3/bin/conda install -y -f -c bioconda flash=1.2.11

echo "---------- INSTALLING BWA 0.7.15 -----------"
./miniconda3/bin/conda install -y -f -c bioconda bwa=0.7.15

echo "---------- INSTALLING SAMTOOLS 1.3.1 -----------"
./miniconda3/bin/conda install -y -f -c bioconda samtools=1.3.1

cd ..

echo "---------- SETTING PERMISSIONS -----------"
chmod 770 .

chmod 440 install.sh
chmod 440 version_info.txt

chmod -R a-w software/
chmod -R o-rx software/
chmod -R 550 scripts/

chmod 770 rawdata/
chmod 440 rawdata/metadata.txt
chmod 770 rawdata/old_metadata/


chmod -R 770 trimdata/

chmod 770 refgenomes/

echo "---------- CHANGING FOLDER NAME FROM BEARCAVE-master TO BEARCAVE -----------"

cd ..
mv BEARCAVE-master/ BEARCAVE/
cd BEARCAVE/


echo "---------- BEARCAVE INSTALLATION COMPLETE -----------"
echo -e '\n'
echo "You may want to update your working directory by navigating one folder up and back into BEARCAVE/"
echo -e '\n'


