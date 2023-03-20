mkdir data/sumstat

# Publication: Multiple common variants for celiac disease influencing immune gene expression. [Dubois et al., 2010] nature genetics
# celiac disease [CEL] 

wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST000001-GCST001000/GCST000612/harmonised/20190752-GCST000612-EFO_0001060-Build37.f.tsv.gz -O data/sumstat/CEL.tsv.gz

# Publication: Association analyses identify 38 susceptibility loci for inflammatory bowel disease and highlight shared genetic risk across populations. [nature genetics 2015]
# Ulcerative colitis (UC)
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003045/harmonised/26192919-GCST003045-EFO_0000729-Build37.f.tsv.gz -O data/sumstat/UC.tsv.gz

# Crohn's disease (CD)
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003044/harmonised/26192919-GCST003044-EFO_0000384-Build37.f.tsv.gz -O data/sumstat/CD.tsv.gz

# Inflammatory bowel disease (IBD)
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003043/harmonised/26192919-GCST003043-EFO_0003767-Build37.f.tsv.gz -O data/sumstat/IBD.tsv.gz

# Publication: Interpreting type 1 diabetes risk with genetics and single-cell epigenomics [nature 2021]
# T1D
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359-Build38.f.tsv.gz -O data/sumstat/T1D.tsv.gz

