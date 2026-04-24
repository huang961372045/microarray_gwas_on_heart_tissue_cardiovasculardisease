# step1: transform .map and .ped files to .bed file


plink --file M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward --make-bed --out ./plink
# Logging to plink.log.
# Options in effect:
#   --file M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# .ped scan complete (for binary autoconversion).
# Performing single-pass .bed write (696873 variants, 370 people).
# --file: plink.bed + plink.bim + plink.fam written.


# step2: QC


plink --bfile plink --geno 0.05 --mind 0.1 --make-bed --out qcdata # 去除缺失率超过5%的SNP和缺失率超过10%的个体
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to qcdata.log.
# Options in effect:
#   --bfile plink
#   --geno 0.05
#   --make-bed
#   --mind 0.1
#   --out qcdata

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# 696873 variants loaded from .bim file.
# 370 people (266 males, 88 females, 16 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to qcdata.nosex .
# 7 people removed due to missing genotype data (--mind).
# IDs written to qcdata.irem .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 62211 het. haploid genotypes present (see qcdata.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate in remaining samples is 0.976251.
# 45396 variants removed due to missing genotype data (--geno).
# 651477 variants and 363 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to qcdata.bed + qcdata.bim + qcdata.fam ... done.
plink --bfile qcdata --hwe 0.0001 --make-bed --out hwe_filtered # 去除不符合Hardy-Weinberg平衡的SNP [https://gwaslab.org/2021/04/04/%E5%93%88%E8%BF%AA%E6%B8%A9%E4%BC%AF%E6%A0%BC%E5%B9%B3%E8%A1%A1-hardy-weinberg-equilibrium/]
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to hwe_filtered.log.
# Options in effect:
#   --bfile qcdata
#   --hwe 0.0001
#   --make-bed
#   --out hwe_filtered

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# 651477 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to hwe_filtered.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 7120 het. haploid genotypes present (see hwe_filtered.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.995234.
# Warning: --hwe observation counts vary by more than 10%, due to the X
# chromosome.  You may want to use a less stringent --hwe p-value threshold for X
# chromosome variants.
# --hwe: 587 variants removed due to Hardy-Weinberg exact test.
# 650890 variants and 363 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to hwe_filtered.bed + hwe_filtered.bim + hwe_filtered.fam ... done.
plink --bfile hwe_filtered --maf 0.05 --make-bed --out maf_filtered # 去除次要等位基因频率（MAF）小于5%的SNP
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to maf_filtered.log.
# Options in effect:
#   --bfile hwe_filtered
#   --maf 0.05
#   --make-bed
#   --out maf_filtered

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# 650890 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to maf_filtered.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 6793 het. haploid genotypes present (see maf_filtered.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.995244.
# 376474 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 274416 variants and 363 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to maf_filtered.bed + maf_filtered.bim + maf_filtered.fam ... done.


plink --bfile maf_filtered --freq --out freq_stat # MAF 计算
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to freq_stat.log.
# Options in effect:
#   --bfile maf_filtered
#   --freq
#   --out freq_stat

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to freq_stat.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see freq_stat.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# --freq: Allele frequencies (founders only) written to freq_stat.frq .
$ head freq_stat.frq 
 CHR                                                SNP   A1   A2          MAF  NCHROBS
   0                                     GSA-rs75263511    A    G       0.1451      696
   0                                      GSA-rs9301562    C    T       0.2251      724
   0                                          rs4317451    G    A      0.08908      696
   0                                         rs61974570    A    C       0.1563      710
   1                                    GSA-rs116587930    A    G      0.06944      720
   1                                          rs3131972    T    C       0.1751      708
   1                                         rs12127425    A    G      0.05387      724
   1                                          rs4970383    T    G       0.2583      724
   1                                          rs4970382    G    A       0.4022      726
# step3: correlation analysis


# GWAS analysis
plink --bfile maf_filtered --assoc --out assoc_results
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to assoc_results.log.
# Options in effect:
#   --assoc
#   --bfile maf_filtered
#   --out assoc_results

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to assoc_results.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see assoc_results.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# 274416 variants and 363 people pass filters and QC.
# Note: No phenotypes present.
# Warning: Skipping --assoc/--model since less than two phenotypes are present.

# ============================ #
# (Optional) Linear regression with phenotype and covariates
# plink --bfile maf_filtered --linear --pheno pheno.txt --covar covar.txt --out linear_results

# (Optional) Logistic regression with phenotype and covariates
# plink --bfile maf_filtered --logistic --pheno pheno.txt --covar covar.txt --out logistic_results



# ============================ #
# R sctript
# ============================ #


# step4: visulization


# Manhattan plot

# QQ plot



# step4: SNP annotation(ANNOVAR)

wget http://www.openbioinformatics.org/annovar/download/annovar.latest.tar.gz
tar -zxvf annovar.latest.tar.gz
cd annovar

plink --bfile maf_filtered --recode vcf --out wgas1
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to wgas1.log.
# Options in effect:
#   --bfile maf_filtered
#   --out wgas1
#   --recode vcf

# 64295 MB RAM detected; reserving 32147 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to wgas1.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see wgas1.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# 274416 variants and 363 people pass filters and QC.
# Note: No phenotypes present.
# --recode vcf to wgas1.vcf ... done.
# Warning: At least one VCF allele code violates the official specification;
# other tools may not accept the file.  (Valid codes must either start with a
# '<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or
# represent a breakend.)

# 下载RefGene注释数据库
./annotate_variation.pl -buildver hg38 -downdb refGene humandb/

# 下载其他需要的数据库，如cytoBand, dbnsfp, clinvar, gnomad等
./annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
./annotate_variation.pl -buildver hg38 -downdb dbnsfp35a humandb/
./annotate_variation.pl -buildver hg38 -downdb clinvar_20210501 humandb/
./annotate_variation.pl -buildver hg38 -downdb gnomad211_genome humandb/


# NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz ... OK
# NOTICE: Uncompressing downloaded files
# NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the 'humandb' directory
# ---------------------------ADDITIONAL PROCEDURE---------------------------
# --------------------------------------------------------------------------
# NOTICE: the FASTA file for the genome is not available to download but can be generated by the ANNOVAR software.
# PLEASE RUN THE FOLLOWING TWO COMMANDS CONSECUTIVELY TO GENERATE THE FASTA FILES (you may need to change -seqdir to -seqfile for some genomes):

#         annotate_variation.pl --buildver hg38 --downdb seq humandb/hg38_seq
#         retrieve_seq_from_fasta.pl humandb/hg38_refGene.txt -seqdir humandb/hg38_seq -format refGene -outfile humandb/hg38_refGeneMrna.fa
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz ... OK
# NOTICE: Uncompressing downloaded files
# NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the 'humandb' directory
# NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/dbnsfp35a.txt.gz ... Failed
# WARNING: Some files cannot be downloaded, including http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/dbnsfp35a.txt.gz
# NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/clinvar_20210501.txt.gz ... Failed
# WARNING: Some files cannot be downloaded, including http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/clinvar_20210501.txt.gz
# NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gnomad211_genome.txt.gz ... Failed
# WARNING: Some files cannot be downloaded, including http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gnomad211_genome.txt.gz


./annotate_variation.pl --buildver hg38 --downdb seq humandb/hg38_seq
# NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.zip ... Failed
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz ... Failed
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz ... 
# OK
# NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ... OK
# NOTICE: Uncompressing downloaded files
# NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the 'humandb/hg38_seq' directory
./retrieve_seq_from_fasta.pl humandb/hg38_refGene.txt -seqdir humandb/hg38_seq -format refGene -outfile humandb/hg38_refGeneMrna.fa