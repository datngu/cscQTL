# cscQTL: An integrative framework for circular RNA quantitative trait locus discovery

  


## 1. Introduction


cscQTL is an integrative framework for circular RNA quantitative trait locus (circQTL) discovery. The motivation for this framework is current circQTL studies relies hevily on a single circRNA calling method to annotate and quantify circRNAs. Regardless the efforts of state-of-the-art circRNA calling methods, circRNA detection is still suffer from a certain amount of false positive and circRNA detection exhibits little agreement between calling tools that cause divergence results in circQTL downstream analyses. cscQTL resolves the problem by combining inputs of serveal circRNA calling algorithm and re-quantify circRNA expression with Quasi-mapping by construction the pseudo circRNA references. cscQTL is implemented as an ready to use pipeline based on Nextflow. It also intergrate an automatic procedure to perform collocation with COLOC. For comparison purpose, the tool also allow users perform single method circQTL with Circall, CIRI2, CIRCExplorer2.



## 2. Installation
  

### 2.0 Hardware requirement

 

The pipeline requires at least a linux system with 32 CPUs and 64GB of RAM to run.

 

### 2.1 Dependency

  

Since cscQTL is implemted with Nextflow (DSL2), you would need Nextflow to run it. Further information to install Nextflow can be found on its home page:[https://www.nextflow.io/](https://www.nextflow.io/)

  

For tool dependency, I have placed all dependencies in their corresponding container, so you would need to run with Docker or Singularity.

  

The default desgins of cscQTL is running with either Docker using a single machine (**-profile standard**) or with sigularity using slurm HPC(**-profile cluster**), you may need to change the **nextflow.config** to adapt to the your system.

  

Further customization, I recommend to consult Nextflow homepage: [https://www.nextflow.io/docs/latest/config.html](https://www.nextflow.io/docs/latest/config.html)

  

  

### 2.2 Clone the pipeline

  

```sh

 
git clone  https://github.com/datngu/cscQTL.git


```

 

## 3. Parameters

  
| Parameters          | Description                                                                                                         | Default setting                        |
| :------------------ | :------------------------------------------------------------------------------------------------------------------ | :------------------------------------- |
| genome              | genome in fasta format                                                                                              | $baseDir/data/ref/genome.fa            |
| cdna                | transcripts (cDNA) in fasta format                                                                                  | $baseDir/data/ref/cdna.fa              |
| bsj                 | back-splice junction database of Circall in fastq format                                                            | $baseDir/data/ref/bsj.fa               |
| annotation          | ensembl annotation in gtf format                                                                                    | $baseDir/data/ref/annotation.gtf       |
| reads               | ribo minus RNA seq reads in fastq.gz format                                                                         | $baseDir/data/reads/\*\_{1,2}.fastq.gz |
| genotype            | genotype in vcf.gz format                                                                                           | $baseDir/data/genotype.vcf.gz          |
| meta                | meta data in csv format - using for matching covariates and rna file names                                          | $baseDir/data/meta.csv                 |
| sumstat\_files      | summary statistics for COLOC analyses                                                                               | $baseDir/data/sumstat/\*.gz            |
| meta\_sumstat       | summary statistics information - for matching file names and detail sample sizes                                    | $baseDir/data/meta\_sumstat.csv        |
| trace\_dir          | directory for tracing output - all intermediate files in the analyses will be soft-linked here - used for debugging | $baseDir/trace\_dir                    |
| outdir              | directory of output                                                                                                 | $baseDir/results                       |
| Running parameters  |                                                                                                                     |                                        |
| chrom               | chromsome range used for QTL mapping                                                                                | 1\..22                                 |
| peer                | PEER range used for optimize the number of PEER factors                                                             | 1\..20                                 |
| genotype\_PCs       | number of genotype principle components used for QTL mapping                                                        | 4                                      |
| bsj\_filter         | BSJ read cutoff for circRNA candidate filtering                                                                     | 2                                      |
| consensus           | number of circRNA calling methods requires to accept a circRNA candidate for consensus filtering                    | 1                                      |
| exp\_prop           | population expression cutoff for filtering circRNAs                                                                 | 0\.3                                   |
| maf                 | MAF cutoff for genotype filtering                                                                                   | 0\.05                                  |
| fdr                 | q-value cutoff                                                                                                      | 0\.05                                  |
| fastqtl\_window     | FastQTL window size                                                                                                 | 1\.00E+06                              |
| Pipeline parameters |                                                                                                                     |                                        |
| cscqtl              | running consensus-based circQTL                                                                                     | true                                   |
| coloc               | running collocation analysis - summary statistics required                                                          | false                                  |
| circall             | running circQTL with Circall                                                                                        | false                                  |
| ciri2               | running circQTL with CIRI2                                                                                          | false                                  |
| circexplorer2       | running circQTL with CIRCExplorer2                                                                                  | false                                  |



## 4. Input data explaination


We prepared a script to download the needed annotation of human genome hg38, annotation v106. Note: the script will automaticly places the files in **data/ref**

  

```sh


bash download_hg38_annotation.sh
  

```



You can also download the example sumary statistics for collocation analyses. Note: the script will automaticly places the files in **data/sumstat**


```sh


bash dowload_sumstats.sh


```



For other input files, we prepare example input data for your convenient to figure how the tool works:

**It is noted that these data are not the real data - they are for explaination purpose only**


- example of reads data : data/reads

The input argument to the nextflow pipeline should be: data/reads/\*\_{1,2}.fastq.gz for this directory. It is noted that, all RNA seq sample must be named in the same format, the "\*" is the wildcard of all sample ID.

- example of genotype data : data/genotype.vcf.gz

A standard vcf.gz file is accepted. rsID annotation is required - for collocation testing).


- example of meta.csv: data/meta.csv

The **meta.csv** file is very important and need to pay careful attention. Its need to includes at least 2 colums: **rna_id, genotype_id**. The column names must be exacted, and separator must be ",". The purpose of this file is to correctly mapping genotype ID and RNA ID. You can also add orther columns, they will be included as numeric covariates. Categorical variable such as Gender must be one-hot encoded in this file.

- example of meta_sumstat.csv: data/meta_sumstat.csv

This file is to provide the tools file names of summary statistics and its case-control ratio. The column names must be exact: file_name, case_prop.


## 5. Real data analysis with 40 t-cell dataset

Bellow is the real scripts I used to apply the tool for the 40 t-cell dataset.

```sh
### select the working directory - you may need to custome yourself
cd /sigma4/projects/

# clonning the pipeline
git clone  https://github.com/datngu/cscQTL.git

# change to the pipeline directory
cd cscQTL

# download the annotation - the files will be placed to matched the default setting
bash download_hg38_annotation.sh

# download the summary statistics - the files will be placed to matched the default setting
bash dowload_sumstats.sh

# the real RNA seqs of 40 T-cell ribo minus: can be download from https://www.synapse.org/#!Synapse:syn22250947
# you would need to rename the RNA seq file names to be the similar names at data/reads

reads="/sigma4/data/40T-cell-blueprint/RNA_seq_data/total_RNA/*_{1,2}.fastq.gz"

# you would need to change "/sigma4/data/40T-cell-blueprint/RNA_seq_data/total_RNA" to your own path

# the genotype data, I provide the hg38 genome (converted to hg38 with crossmap)
wget https://github.com/datngu/data/releases/download/v.0.0.2/tcel_hg38.vcf.gz -O tcel_hg38.vcf.gz
genotype="$PWD/tcel_hg38.vcf.gz"

```



If you run in a local computer with Docker:

```sh

nextflow run main.nf -resume --reads $reads --outdir "TEST_WITH_LOCAL" --genotype $genotype --consensus 3 --coloc true --circall true --ciri2 true --circexplorer2 true -with-report -profile standard


```

If you run in a HCP server with Singularity and Slurm:

```sh

nextflow run main.nf -resume --reads $reads --outdir "TEST_WITH_HPC" --genotype $genotype --consensus 3 --coloc true --circall true --ciri2 true --circexplorer2 true -with-report -profile cluster


```



## 6. License
  
cscQTL uses GNU General Public License GPL-3.


## 7. References

Nguyen, D. T. et al. Circall: fast and accurate methodology for discovery of circular rnas from paired-end rna-sequencing data. BMC bioinformatics 22, 1–18 (2021).

Gao, Y., Zhang, J. & Zhao, F. Circular rna identification based on multiple seed matching. Briefings in bioinformatics 19, 803–810 (2018).

Zhang, X.-O. et al. Diverse alternative back-splicing and alternative splicing landscape of circular rnas. Genome research 26, 1277–1287 (2016).

Ongen, H., Buil, A., Brown, A. A., Dermitzakis, E. T. & Delaneau, O. Fast and efficient qtl
mapper for thousands of molecular phenotypes. Bioinformatics 32, 1479–1485 (2016).

Stegle, O., Parts, L., Piipari, M., Winn, J. & Durbin, R. Using probabilistic estimation of expression residuals (peer) to obtain increased power and interpretability of gene expression analyses. Nature protocols 7, 500–507 (2012).


Hormozdiari, F. et al. Colocalization of gwas and eqtl signals detects target genes. The
American Journal of Human Genetics 99, 1245–1260 (2016).

Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nature biotechnology 35, 316–319 (2017).


Chang, C. C. et al. Second-generation plink: rising to the challenge of larger and richer
datasets. Gigascience 4, s13742–015 (2015).

Storey, J. D. & Tibshirani, R. Statistical significance for genomewide studies. Proceedings of the National Academy of Sciences 100, 9440–9445 (2003).





















 