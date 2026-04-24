# step1: transform .map and .ped files to .bed file ====================================================================================================

# 使用新的表型文件进行分析
plink --file M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward --pheno pheno.txt --make-bed --out phenotype_data
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to phenotype_data.log.
# Options in effect:
#   --file M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward
#   --make-bed
#   --out phenotype_data
#   --pheno pheno.txt

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# .ped scan complete (for binary autoconversion).
# Performing single-pass .bed write (696873 variants, 370 people).
# --file: phenotype_data-temporary.bed + phenotype_data-temporary.bim +
# phenotype_data-temporary.fam written.
# 696873 variants loaded from .bim file.
# 370 people (266 males, 88 females, 16 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to phenotype_data.nosex .
# 364 phenotype values present after --pheno.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 370 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 62211 het. haploid genotypes present (see phenotype_data.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.973352.
# 696873 variants and 370 people pass filters and QC.
# Phenotype data is quantitative.
# --make-bed to phenotype_data.bed + phenotype_data.bim + phenotype_data.fam ...
# done.

# step2: QC ====================================================================================================


plink --bfile phenotype_data --geno 0.05 --mind 0.1 --make-bed --out qcdata # 去除缺失率超过5%的SNP和缺失率超过10%的个体
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to qcdata.log.
# Options in effect:
#   --bfile phenotype_data
#   --geno 0.05
#   --make-bed
#   --mind 0.1
#   --out qcdata

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 696873 variants loaded from .bim file.
# 370 people (266 males, 88 females, 16 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to qcdata.nosex .
# 364 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
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
# Phenotype data is quantitative.
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

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 651477 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to hwe_filtered.nosex .
# 359 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
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
# Phenotype data is quantitative.
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

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 650890 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to maf_filtered.nosex .
# 359 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
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
# Phenotype data is quantitative.
# --make-bed to maf_filtered.bed + maf_filtered.bim + maf_filtered.fam ... done.

# step3: MAF calculation ====================================================================================================

plink --bfile maf_filtered --freq --out freq_stat # MAF calculation
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to freq_stat.log.
# Options in effect:
#   --bfile maf_filtered
#   --freq
#   --out freq_stat

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to freq_stat.nosex .
# 359 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
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



# step4: Population clustring analysis（PCA）  ====================================================================================================

plink --bfile maf_filtered --indep-pairwise 200 50 0.25 --out indep_snp #  LD pruning [https://gwaslab.org/2021/05/18/pruning-and-clumping/]
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to indep_snp.log.
# Options in effect:
#   --bfile maf_filtered
#   --indep-pairwise 200 50 0.25
#   --out indep_snp

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to indep_snp.nosex .
# 359 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see indep_snp.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# 274416 variants and 363 people pass filters and QC.
# Phenotype data is quantitative.
# --indep-pairwise: Ignoring 4 chromosome 0 variants.
# Pruned 12188 variants from chromosome 1, leaving 8349.
# Pruned 13146 variants from chromosome 2, leaving 8180.
# Pruned 10705 variants from chromosome 3, leaving 7021.
# Pruned 9284 variants from chromosome 4, leaving 6610.
# Pruned 8933 variants from chromosome 5, leaving 6318.
# Pruned 13818 variants from chromosome 6, leaving 6263.
# Pruned 8445 variants from chromosome 7, leaving 5805.
# Pruned 8205 variants from chromosome 8, leaving 5273.
# Pruned 6516 variants from chromosome 9, leaving 4784.
# Pruned 7964 variants from chromosome 10, leaving 5449.
# Pruned 8219 variants from chromosome 11, leaving 5152.
# Pruned 7656 variants from chromosome 12, leaving 5222.
# Pruned 5725 variants from chromosome 13, leaving 3948.
# Pruned 4868 variants from chromosome 14, leaving 3611.
# Pruned 4861 variants from chromosome 15, leaving 3465.
# Pruned 4896 variants from chromosome 16, leaving 3583.
# Pruned 3878 variants from chromosome 17, leaving 3223.
# Pruned 4415 variants from chromosome 18, leaving 3501.
# Pruned 2373 variants from chromosome 19, leaving 2537.
# Pruned 3614 variants from chromosome 20, leaving 2812.
# Pruned 1995 variants from chromosome 21, leaving 1647.
# Pruned 1731 variants from chromosome 22, leaving 1708.
# Pruned 12143 variants from chromosome 23, leaving 3541.
# Pruned 355 variants from chromosome 24, leaving 7.
# Pruned 193 variants from chromosome 25, leaving 234.
# Pruned 33 variants from chromosome 26, leaving 10.
# Pruning complete.  166159 of 274412 variants removed.
# Marker lists written to indep_snp.prune.in and indep_snp.prune.out .

plink --bfile maf_filtered --extract indep_snp.prune.in --pca 10 --out pca_results # use pruned SNPs to do PCA analysis
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to pca_results.log.
# Options in effect:
#   --bfile maf_filtered
#   --extract indep_snp.prune.in
#   --out pca_results
#   --pca 10

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to pca_results.nosex .
# 359 phenotype values loaded from .fam.
# --extract: 108253 variants remaining.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using up to 8 threads (change this with --threads).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 211 het. haploid genotypes present (see pca_results.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993926.
# 108253 variants and 363 people pass filters and QC.
# Phenotype data is quantitative.
# Excluding 3558 variants on non-autosomes from relationship matrix calc.
# Relationship matrix calculation complete.
# --pca: Results saved to pca_results.eigenval and pca_results.eigenvec .


# step4: correlation analysis  ====================================================================================================
# plink --bfile maf_filtered --pheno pheno.txt --covar pca_results.eigenvec --covar-number 1-3 --assoc --out dcm_control_assoc_adjusted
# plink --bfile maf_filtered --pheno pheno.txt --covar pca_results.eigenvec --covar-number 1-3 --linear --allow-no-sex  --out all_control_assoc_adjusted
# plink --bfile maf_filtered --pheno pheno_DCM.txt --covar pca_results.eigenvec --linear --allow-no-sex  --out dcm_control_assoc_adjusted

plink --bfile maf_filtered --pheno pheno.txt --covar pca_results.eigenvec --linear --allow-no-sex  --out dcm_control_assoc_adjusted
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to dcm_control_assoc_adjusted.log.
# Options in effect:
#   --allow-no-sex
#   --bfile maf_filtered
#   --covar pca_results.eigenvec
#   --linear
#   --out dcm_control_assoc_adjusted
#   --pheno pheno.txt

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to dcm_control_assoc_adjusted.nosex .
# 359 phenotype values present after --pheno.
# Using 1 thread (no multithreaded calculations invoked).
# --covar: 10 covariates loaded.
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see dcm_control_assoc_adjusted.hh
# ); many commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# 274416 variants and 363 people pass filters and QC.
# Phenotype data is quantitative.
# Writing linear model association results to
# dcm_control_assoc_adjusted.assoc.linear ... done.
plink --bfile maf_filtered --pheno pheno.txt --linear --allow-no-sex  --out control_assoc_adjusted_npPCA
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to control_assoc_adjusted_npPCA.log.
# Options in effect:
#   --allow-no-sex
#   --bfile maf_filtered
#   --linear
#   --out control_assoc_adjusted_npPCA
#   --pheno pheno.txt

# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to control_assoc_adjusted_npPCA.nosex .
# 359 phenotype values present after --pheno.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see
# control_assoc_adjusted_npPCA.hh ); many commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# 274416 variants and 363 people pass filters and QC.
# Phenotype data is quantitative.
# Writing linear model association results to
# control_assoc_adjusted_npPCA.assoc.linear ... done.



# step5: python generate VCF files  ====================================================================================================

''' Python Code'''
# 读取关联分析结果
assoc_results = pd.read_csv(dir_in+'dcm_control_assoc_adjusted.assoc.linear', sep='\s+')

# 按P值排序
assoc_results = assoc_results.sort_values('P',ascending=True)
# 设置显著性阈值（比如Bonferroni校正或其他标准）
bonferroni_threshold = 0.05 / len(assoc_results)  # Bonferroni校正
suggestive_threshold = 1e-4  # 建议性阈值

# 设置BETA和STAT的阈值
beta_threshold = 0.2  # 例如，BETA大于0.2
stat_threshold = 2.0   # 例如，STAT大于2.0
significant_snps = assoc_results[
    (assoc_results['P'] < suggestive_threshold) &
    (assoc_results['BETA'].abs() > beta_threshold) &
    (assoc_results['STAT'] > stat_threshold)
]

# 输出显著 SNP 的数量
print(f"Number of significant SNPs: {len(significant_snps)}")
significant_snps
significant_snps.to_csv(dir_in+'significant_snps', index=False, sep='\t')


assoc_results = pd.read_csv(dir_in+'dcm_control_assoc_adjusted.assoc', sep='\s+')
assoc_results
significant_snps['SNP']
df_assoc = assoc_results[assoc_results['SNP'].isin(significant_snps['SNP'])]
df_assoc

# 创建一个输出VCF文件
vcf_file = f'{dir_out}significant_SNP.vcf'

# VCF文件头部信息
vcf_header = """##fileformat=VCFv4.2
##source=PLINK
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
"""

# 写入VCF文件头
with open(vcf_file, 'w') as f:
    f.write(vcf_header)

    # 遍历每一行，提取关键字段写入VCF文件
    for index, row in df_assoc.iterrows():
        chrom = row['CHR']
        pos = row['BP']
        snp_id = row['SNP']
        ref = row['A1']  # A1通常作为参考等位基因
        alt = row['A2']  # A2作为变异等位基因
        qual = '.'       # 质量字段暂时用'.'占位
        filter = '.'     # 过滤字段暂时用'.'占位
        info = '.'       # INFO字段可以根据需求进一步填写
        
        # 写入VCF文件的一行
        f.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n")

print(f"VCF file '{vcf_file}' has been created.")
''' End '''
# step6: ANNOVAR  ====================================================================================================
grep -v "^#" output.vcf > output_noheader.vcf
convert2annovar.pl -format vcf4 output_noheader.vcf > input.avinput
# NOTICE: Finished reading 17 lines from VCF file
# NOTICE: A total of 17 locus in VCF file passed QC threshold, representing 17 SNPs (14 transitions and 3 transversions) and 0 indels/substitutions
# NOTICE: Finished writing 17 SNP genotypes (14 transitions and 3 transversions) and 0 indels/substitutions for 1 sample

table_annovar.pl input.avinput ../../../data/annovar_dataset/humandb -buildver hg19 -out output --nopolish -remove -protocol refGene,cytoBand -operation g,r -nastring . -vcfinput
# NOTICE: Running with system command <convert2annovar.pl  -includeinfo -allsample -withfreq -format vcf4 input.avinput > output.avinput>
# NOTICE: Finished reading 27 lines from VCF file
# NOTICE: A total of 27 locus in VCF file passed QC threshold, representing 27 SNPs (21 transitions and 6 transversions) and 0 indels/substitutions
# NOTICE: Finished writing allele frequencies based on 27 SNP genotypes (21 transitions and 6 transversions) and 0 indels/substitutions for 0 samples

# NOTICE: Running with system command </home/zihuang/data2/tools/annovar/table_annovar.pl output.avinput ../../../data/annovar_dataset/humandb -buildver hg19 -outfile output --nopolish -remove -protocol refGene,cytoBand -operation g,r -nastring . -otherinfo>
# -----------------------------------------------------------------
# NOTICE: Processing operation=g protocol=refGene

# NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile output.refGene -exonsort -nofirstcodondel output.avinput ../../../data/annovar_dataset/humandb>
# NOTICE: Output files are written to output.refGene.variant_function, output.refGene.exonic_variant_function
# NOTICE: Reading gene annotation from ../../../data/annovar_dataset/humandb/hg19_refGene.txt ... Done with 78239 transcripts (including 18578 without coding sequence annotation) for 28293 unique genes
# NOTICE: Processing next batch with 27 unique variants in 27 input lines
# -----------------------------------------------------------------
# NOTICE: Processing operation=r protocol=cytoBand

# NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 -outfile output output.avinput ../../../data/annovar_dataset/humandb>
# NOTICE: Output file is written to output.hg19_cytoBand
# NOTICE: Reading annotation database ../../../data/annovar_dataset/humandb/hg19_cytoBand.txt ... Done with 862 regions
# NOTICE: Finished region-based annotation on 27 genetic variants
# -----------------------------------------------------------------
# NOTICE: Multianno output file is written to output.hg19_multianno.txt
# NOTICE: Reading from output.hg19_multianno.txt
# -----------------------------------------------------------------
# NOTICE: VCF output is written to output.hg19_multianno.vcf


output.hg19_multianno.vcf is the final resultr
# step7: convert output.hg19_multianno.vcf to csv  ====================================================================================================

# step4: correlation analysis  ====================================================================================================
# plink --bfile maf_filtered --assoc --out assoc_results # GWAS analysis
# PLINK v1.90b7.2 64-bit (11 Dec 2023)
# Options in effect:
#   --bfile maf_filtered
#   --freq
#   --out freq_stat

# Hostname: srvcorem2
# Working directory: /data2/core-med1/zhuang/code/20240516_microarray_fromhomefile/M01144_01-04_PLINK_copy
# Start time: Thu Oct 24 14:18:01 2024

# Random number seed: 1729772281
# 15953 MB RAM detected; reserving 7976 MB for main workspace.
# 274416 variants loaded from .bim file.
# 363 people (266 males, 88 females, 9 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to freq_stat.nosex .
# 359 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 363 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 782 het. haploid genotypes present (see freq_stat.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993941.
# --freq: Allele frequencies (founders only) written to freq_stat.frq .

# End time: Thu Oct 24 14:18:01 2024

# ============================ #
# (Optional) Linear regression with phenotype and covariates
# plink --bfile maf_filtered --linear --pheno pheno.txt --covar covar.txt --out linear_results

# (Optional) Logistic regression with phenotype and covariates
# plink --bfile maf_filtered --logistic --pheno pheno.txt --covar covar.txt --out logistic_results



# ============================ #
# R sctript
# ============================ #
# step4: visulization ====================================================================================================


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






# 分层分析
# 对每个主要疾病类型单独进行分析：

# bashCopy# DCM组与对照组比较
# plink --bfile phenotype_data --pheno pheno.txt --pheno-name DCM --make-bed --out dcm_analysis

# # ICM组与对照组比较
# plink --bfile phenotype_data --pheno pheno.txt --pheno-name ICM --make-bed --out icm_analysis

# 多表型文件
# 创建一个包含多个表型的文件：

# bashCopy# 创建多表型文件(multiple_pheno.txt):
# # FID IID DCM ICM HCM ARVC
# # 每个表型用0/1编码
# 使用命令：
# bashCopyplink --bfile phenotype_data --mpheno 1 --make-bed --out analysis_pheno1
# plink --bfile phenotype_data --mpheno 2 --make-bed --out analysis_pheno2

# 协变量分析
# 如果需要考虑其他因素（如年龄、性别等）：

# bashCopy# 创建协变量文件(covar.txt)
# plink --bfile phenotype_data --covar covar.txt --covar-name AGE,SEX --assoc --out results_with_covariates
# 建议的分析步骤：

# 首先将表型分为主要几类：


# DCM及其亚型（包括postmyokarditisch等）
# ICM
# HCM/HOCM
# ARVC/ARVD
# 其他心肌病
# 对照组(N)


# 创建相应的表型文件
# 对每个主要疾病组进行：


# 质量控制
# 关联分析
# 差异分析

# 你想重点分析哪些表型组？或者你更倾向于用哪种方式来处理这些多表型数据？我可以帮你设计更详细的分析流程。