#!/bin/bash

basedir="/home/mounika/bioinformatics/test_project"

wget -P $basedir "file url"

./guppy_basecaller -i $basedir -o "output folder location" -c "config file"

pycoQC -f "folder containing sequencing_summary.txt from guppy basecalling" -o "output folder to report.html"

minimap2 -a "fastq file from guppy" "fasta file" > output sam file

 
