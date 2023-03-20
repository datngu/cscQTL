#!/bin/bash

cd data/ref

echo "Downloading annotation hg38 version 106...."

wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz -O annotation.gtf.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O cdna.fa.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O genome.fa.gz
wget https://github.com/datngu/Circall/releases/download/v1.0.0/Homo_sapiens.GRCh38.106_BSJ_sequences.fa.gz -O bsj.fa.gz


echo "Unzip...!"

gunzip annotation.gtf.gz
gunzip cdna.fa.gz
gunzip genome.fa.gz
gunzip bsj.fa.gz

echo "DONE!"



