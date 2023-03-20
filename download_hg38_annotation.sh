#!/bin/bash

mkdir -p data/ref

echo "Downloading annotation hg38 version 106...."

wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz -O data/ref/annotation.gtf.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O data/ref/cdna.fa.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O data/ref/genome.fa.gz
wget https://github.com/datngu/Circall/releases/download/v1.0.0/Homo_sapiens.GRCh38.106_BSJ_sequences.fa.gz -O data/ref/bsj.fa.gz


echo "Unzip...!"

gunzip data/ref/annotation.gtf.gz
gunzip data/ref/cdna.fa.gz
gunzip data/ref/genome.fa.gz
gunzip data/ref/bsj.fa.gz

echo "DONE!"



