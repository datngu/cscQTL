#!/bin/bash

#input parameters are:
    #txIdx -- quasi-index of txFasta
    #bsjIdx -- quasi-index of BSJ reference fasta file
    #read1 -- input read1
    #read2 -- input read2
    #p -- number of thread
    #o -- output folder
#

#export LD_LIBRARY_PATH=/path/to/linux/lib:$LD_LIBRARY_PATH
#export PATH=/path/to/linux/bin:$PATH

syntax="./Circall_recount.sh \n\t-txIdx [quasi-index of txFasta] \n\t-bsjIdx [quasi-index of BSJ fasta file] \n\t-read1 [read1 fastq.gz file] \n\t-read2 [read2 fastq.gz file] \n\t-thread [number of thread] \n\t-o [output prefix]"

#add path
#export LC_ALL=C
#export LD_LIBRARY_PATH=/path/to/linux/lib:$LD_LIBRARY_PATH
#export PATH=/path/to/linux/bin:$PATH


while [[ $# -gt 1 ]]
do
key="$1"

case $key in


     -txIdx|--txIdx) 
     IndexTranscriptome=$(readlink -f $2)
     shift
     ;;

     -bsjIdx|--bsjIdx) 
     IndexBSJ=$(readlink -f $2)
     shift
     ;;


     -read1|--read1) 
     inRead1=$(readlink -f $2)
     shift
     ;;

     -read2|--read2) 
     inRead2=$(readlink -f $2)
     shift
     ;;

     -p|--thread) 
     CPUNUM=$2
     shift
     ;; 


     -o|--out)
     outFile=$2
     shift
     ;;

     *)

esac
shift
done



if [[ -z "$IndexTranscriptome" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no index directory of the transcript sequences (cDNA) !"
   echo ""
   exit
fi

if [[ -z "$IndexBSJ" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no index directory of back-splicing-junction (BSJ) reference !"
   echo ""
   exit
fi


if [[ -z "$inRead1" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no input read1 file !"
   echo ""
   exit
fi


if [[ -z "$inRead2" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no input read2 file !"
   echo ""
   exit
fi



if [[ -z "$CPUNUM" ]]; then
   CPUNUM=4
   echo "No specific setting for thread number, so use the default setting (-p 4)"
fi


if [[ -z "$outFile" ]]; then
   outFile="circRNA_recount.txt"
   echo "No specific setting for thread number, so use the default setting: <circRNA_recount.txt>"
fi

used_agrs=("IndexTranscriptome" "IndexBSJ" "inRead1" "inRead2" "CPUNUM" "outFile")
inFileFormat=${inRead1##*.}

echo -e "
---------------------------------------------------------------------------
                     Checking arguments successfully!!
---------------------------------------------------------------------------
Your input arguments list:
"

for (( i=0; i< ${#in_agr[@]} ; i++ ));
do
  echo "${in_agr[$i]} was specified as: ${!used_agrs[$i]}"
done

echo ""


## for testing
## docker run -it --rm -v /sigma4:/sigma4 --name rna ndatth/rna-tools:v0.0.0
# cd /sigma4/data/genome_ref_hg38/test
# IndexTranscriptome=/sigma4/data/genome_ref_hg38/IndexTranscriptome
# IndexBSJ=/sigma4/projects/rnaQTL_hg38/trace_dir/recount_consensus_idx/index_pseudo_seqs
# inRead1=/sigma4/data/40T-cell-blueprint/RNA_seq_data/total_RNA/EGAR00001193012_1.fastq.gz
# inRead2=/sigma4/data/40T-cell-blueprint/RNA_seq_data/total_RNA/EGAR00001193012_2.fastq.gz
# CPUNUM=56
# outFile=sample12.txt


working_dir=${outFile}_tem_dir
mkdir ${working_dir}
# cd ${working_dir}

Circall_log=${outFile}.log

# ####### Extract upmapped reads
Circall_wt -i $IndexTranscriptome -1 <(gunzip -c $inRead1) -2 <(gunzip -c $inRead2) -o ${working_dir}/Circall_wt -p $CPUNUM >>$Circall_log 2>&1

# Merge the outputs from multiple cores
cat ${working_dir}/Circall_wt/UN_read1_* >  ${working_dir}/Unmapped_read_1.fa
cat ${working_dir}/Circall_wt/UN_read2_* >  ${working_dir}/Unmapped_read_2.fa
cat ${working_dir}/Circall_wt/fragmentInfo.txt > ${working_dir}/Circall_wt_fragmentInfo.txt

# clear space
rm -rf ${working_dir}/Circall_wt

## mapping to BSJ reference
Circall_bsj -i $IndexBSJ -1 ${working_dir}/Unmapped_read_1.fa -2 ${working_dir}/Unmapped_read_2.fa -o ${working_dir}/Circall_bsj -p $CPUNUM >>$Circall_log 2>&1
cat ${working_dir}/Circall_bsj/UN_read* > ${working_dir}/BSJ_mapped_read.txt

# clear space
rm -rf ${working_dir}/Circall_bsj
rm ${working_dir}/Unmapped_read_1.fa ${working_dir}/Unmapped_read_2.fa

## moving output out
mv ${working_dir}/Circall_wt_fragmentInfo.txt ${outFile}_Circall_wt_fragmentInfo.txt
mv ${working_dir}/BSJ_mapped_read.txt ${outFile}_BSJ_mapped_read.txt
# clear space
rm -rf ${working_dir}


#recount_BSJ.py --bsj ${working_dir}/BSJ_mapped_read.txt --info ${working_dir}/Circall_wt_fragmentInfo.txt --out $outFile

## bash Circall_recount.sh -p 56 -txIdx /sigma4/data/genome_ref_hg38/IndexTranscriptome -bsjIdx /sigma4/projects/rnaQTL_hg38/trace_dir/recount_consensus_idx/index_pseudo_seqs -read1 /sigma4/data/40T-cell-blueprint/RNA_seq_data/total_RNA/EGAR00001193012_1.fastq.gz -read2 /sigma4/data/40T-cell-blueprint/RNA_seq_data/total_RNA/EGAR00001193012_2.fastq.gz

