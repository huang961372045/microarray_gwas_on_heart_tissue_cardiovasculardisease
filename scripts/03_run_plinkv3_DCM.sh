# step0: 准备新的phenotype数据，只保留DCM进行分析 ====================================================================================================
# Extract DCM and Control samples from the pheno.txt 直接从 pheno.txt 提取 DCM 和 Control 样本
awk 'NR==1 {print} 
     NR>1 && ($3==1 || $3==2) {print $1, $2, $3}' pheno.txt > dcm_control_pheno.txt
# Slice the ra w data 原始数据切片
plink --file M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward \
      --keep dcm_control_samples.txt \
      --make-bed \
      --out DCM_vs_Control_FINAL

# step1: transform .map and .ped files to .bed file ====================================================================================================
# 使用新的表型文件进行分析
# 创建最终数据集
plink --bfile DCM_vs_Control_FINAL \
      --pheno dcm_control_pheno.txt \
      --make-bed \
      --out DCM_vs_Control_FINAL_with_pheno
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to DCM_vs_Control_FINAL_with_pheno.log.
# Options in effect:
#   --bfile DCM_vs_Control_FINAL
#   --make-bed
#   --out DCM_vs_Control_FINAL_with_pheno
#   --pheno dcm_control_pheno.txt

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 696873 variants loaded from .bim file.
# 214 people (158 males, 49 females, 7 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to DCM_vs_Control_FINAL_with_pheno.nosex .
# 214 phenotype values present after --pheno.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 214 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 36872 het. haploid genotypes present (see
# DCM_vs_Control_FINAL_with_pheno.hh ); many commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.975067.
# 696873 variants and 214 people pass filters and QC.
# Among remaining phenotypes, 177 are cases and 37 are controls.
# --make-bed to DCM_vs_Control_FINAL_with_pheno.bed +
# DCM_vs_Control_FINAL_with_pheno.bim + DCM_vs_Control_FINAL_with_pheno.fam ...
# done.

# step2: QC ====================================================================================================

plink --bfile DCM_vs_Control_FINAL_with_pheno --geno 0.05 --mind 0.1 --make-bed --out qcdata # 去除缺失率超过5%的SNP和缺失率超过10%的个体
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to qcdata.log.
# Options in effect:
#   --bfile DCM_vs_Control_FINAL_with_pheno
#   --geno 0.05
#   --make-bed
#   --mind 0.1
#   --out qcdata

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 696873 variants loaded from .bim file.
# 214 people (158 males, 49 females, 7 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to qcdata.nosex .
# 214 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# 2 people removed due to missing genotype data (--mind).
# IDs written to qcdata.irem .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 212 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 36872 het. haploid genotypes present (see qcdata.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate in remaining samples is 0.976228.
# 49007 variants removed due to missing genotype data (--geno).
# 647866 variants and 212 people pass filters and QC.
# Among remaining phenotypes, 175 are cases and 37 are controls.
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

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 647866 variants loaded from .bim file.
# 212 people (158 males, 49 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to hwe_filtered.nosex .
# 212 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 212 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 4229 het. haploid genotypes present (see hwe_filtered.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.995764.
# Warning: --hwe observation counts vary by more than 10%.  Consider using
# --geno, and/or applying different p-value thresholds to distinct subsets of
# your data.
# --hwe: 22 variants removed due to Hardy-Weinberg exact test.
# 647844 variants and 212 people pass filters and QC.
# Among remaining phenotypes, 175 are cases and 37 are controls.
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

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 647844 variants loaded from .bim file.
# 212 people (158 males, 49 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to maf_filtered.nosex .
# 212 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 212 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 4229 het. haploid genotypes present (see maf_filtered.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.995764.
# 374231 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 273613 variants and 212 people pass filters and QC.
# Among remaining phenotypes, 175 are cases and 37 are controls.
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

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 273613 variants loaded from .bim file.
# 212 people (158 males, 49 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to freq_stat.nosex .
# 212 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 212 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 790 het. haploid genotypes present (see freq_stat.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.994891.
# --freq: Allele frequencies (founders only) written to freq_stat.frq .
$ head freq_stat.frq 
 CHR                                                SNP   A1   A2          MAF  NCHROBS
   0                                     GSA-rs75263511    A    G       0.1584      404
   0                                      GSA-rs9301562    C    T       0.2146      424
   0                                          rs4317451    G    A       0.0936      406
   0                                         rs61974570    A    C       0.1602      412
   1                                    GSA-rs116587930    A    G      0.07346      422
   1                                          rs3131972    T    C       0.1875      416
   1                                          rs4970383    T    G       0.2642      424
   1                                          rs4970382    G    A       0.4009      424
   1                                          rs3935066    C    T       0.1114      422

 
# step4: 相关性分析 & 个体去重（IBD）  ====================================================================================================

# 性别一致性检查
plink --bfile maf_filtered --check-sex --out sexcheck
# 相关性分析 & 个体去重（IBD）
plink --bfile maf_filtered --genome --min 0.185 --out genome_check

# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to genome_check.log.
# Options in effect:
#   --bfile maf_filtered
#   --genome
#   --min 0.185
#   --out genome_check

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 273613 variants loaded from .bim file.
# 212 people (158 males, 49 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to genome_check.nosex .
# 212 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using up to 8 threads (change this with --threads).
# Before main variant filters, 212 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 790 het. haploid genotypes present (see genome_check.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.994891.
# 273613 variants and 212 people pass filters and QC.
# Among remaining phenotypes, 171 are cases and 36 are controls.  (5 phenotypes
# are missing.)
# Excluding 15907 variants on non-autosomes from IBD calculation.
# IBD calculations complete.
# Finished writing genome_check.genome .
$ cat genome_check.genome 
 FID1       IID1 FID2       IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO
  20        D27  23        D30 UN    NA  0.0000  0.0008  0.9992  0.9996   1  0.999870  1.0000      NA
  77       D100  78       D101 UN    NA  0.0000  0.0006  0.9994  0.9997   1  0.999898  1.0000      NA
 125       D150 229       D262 UN    NA  0.2250  0.5012  0.2738  0.5244   1  0.856255  1.0000 12.5405
 273     D308-2 334       D372 UN    NA  0.0000  0.0003  0.9997  0.9999   1  0.999955  1.0000      NA

# 然后查看 genome_check.genome，选取高相关个体
awk '$10 > 0.185 {print $1, $2; print $3, $4}' genome_check.genome | sort | uniq > high_IBD.txt
$ cat high_IBD.txt
125 D150
20 D27
229 D262
23 D30
273 D308-2
334 D372
77 D100
78 D101
FID1 IID1
FID2 IID2

# 若有高相关个体，用如下命令剔除高相关个体
plink --bfile maf_filtered --remove high_IBD.txt --make-bed --out related_removed
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to related_removed.log.
# Options in effect:
#   --bfile maf_filtered
#   --make-bed
#   --out related_removed
#   --remove high_IBD.txt

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 273613 variants loaded from .bim file.
# 212 people (158 males, 49 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to related_removed.nosex .
# 212 phenotype values loaded from .fam.
# --remove: 204 people remaining.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 204 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 760 het. haploid genotypes present (see related_removed.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate in remaining samples is 0.994882.
# 273613 variants and 204 people pass filters and QC.
# Among remaining phenotypes, 167 are cases and 37 are controls.
# --make-bed to related_removed.bed + related_removed.bim + related_removed.fam
# ... done.

# 主成分分析（PCA）
plink --bfile related_removed --pca 20 --out pca_result
# plink --bfile maf_filtered --extract indep_snp.prune.in --pca 10 --out pca_results # use pruned SNPs to do PCA analysis
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to pca_result.log.
# Options in effect:
#   --bfile related_removed
#   --out pca_result
#   --pca 20

# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 273613 variants loaded from .bim file.
# 204 people (153 males, 46 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to pca_result.nosex .
# 204 phenotype values loaded from .fam.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using up to 8 threads (change this with --threads).
# Before main variant filters, 204 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 760 het. haploid genotypes present (see pca_result.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.994882.
# 273613 variants and 204 people pass filters and QC.
# Among remaining phenotypes, 163 are cases and 36 are controls.  (5 phenotypes
# are missing.)
# Excluding 15907 variants on non-autosomes from relationship matrix calc.
# Relationship matrix calculation complete.
# --pca: Results saved to pca_result.eigenval and pca_result.eigenvec .


# step4: correlation analysis  ====================================================================================================
# plink --bfile maf_filtered --pheno pheno.txt --covar pca_results.eigenvec --covar-number 1-3 --assoc --out dcm_control_assoc_adjusted
# plink --bfile maf_filtered --pheno pheno.txt --covar pca_results.eigenvec --covar-number 1-3 --linear --allow-no-sex  --out all_control_assoc_adjusted
# plink --bfile maf_filtered --pheno pheno_DCM.txt --covar pca_results.eigenvec --linear --allow-no-sex  --out dcm_control_assoc_adjusted
# 关联分析（GWAS）
plink --bfile related_removed \
      --pheno dcm_control_pheno.txt \
      --covar pca_result.eigenvec \
      --covar-number 1-10 \
      --logistic \
      --hide-covar \
      --ci 0.95 \
      --out gwas_result
# PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to gwas_result.log.
# Options in effect:
#   --bfile related_removed
#   --ci 0.95
#   --covar pca_result.eigenvec
#   --covar-number 1-10
#   --hide-covar
#   --logistic
#   --out gwas_result
#   --pheno dcm_control_pheno.txt

# Note: --hide-covar flag deprecated.  Use e.g. "--linear hide-covar".
# 15946 MB RAM detected; reserving 7973 MB for main workspace.
# 273613 variants loaded from .bim file.
# 204 people (153 males, 46 females, 5 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to gwas_result.nosex .
# 204 phenotype values present after --pheno.
# Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
# phenotypes to be ignored, use the --allow-no-sex flag.
# Using 1 thread (no multithreaded calculations invoked).
# --covar: 10 out of 20 covariates loaded.
# Before main variant filters, 204 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 760 het. haploid genotypes present (see gwas_result.hh ); many
# commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.994882.
# 273613 variants and 204 people pass filters and QC.
# Among remaining phenotypes, 163 are cases and 36 are controls.  (5 phenotypes
# are missing.)
# Writing logistic model association results to gwas_result.assoc.logistic ...
# done.

# step5: Results Screening  ====================================================================================================
# 1. 第一步：提取每个样本的显著 SNP（从 GWAS 结果中）
awk '$P < 1e-5' gwas_result.assoc.logistic > significant_snps.txt
## 保留 SNP 一列：
awk '$P < 1e-5 {print $2}' gwas_result.assoc.logistic > sig_snp_ids.txt
$ head sig_snp_ids.txt
SNP
GSA-rs75263511
GSA-rs9301562
rs4317451
rs61974570
GSA-rs116587930
rs3131972
rs4970383
rs4970382
rs3935066

# 2. 第二步：获取 SNP 到基因的注释（Gene Symbol）
# 使用 Ensembl VEP Web 或命令行批量注释 SNP
# 网址：https://www.ensembl.org/Tools/VEP
# 上传 sig_snp_ids.txt
# 下载带 Gene Symbol 的注释表格（如包含 SYMBOL 字段）
#=============#
#Uploaded_variation	Location	Allele	Consequence	IMPACT	SYMBOL	Gene	Feature_type	Feature	BIOTYPE	EXON	INTRON	HGVSc	HGVSp	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	REF_ALLELE	UPLOADED_ALLELE	DISTANCE	STRAND	FLAGS	SYMBOL_SOURCE	HGNC_ID	MANE	MANE_SELECT	MANE_PLUS_CLINICAL	TSL	APPRIS	SIFT	PolyPhen	AF	CLIN_SIG	SOMATIC	PHENO	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	TRANSCRIPTION_FACTORS
rs3131972	1:817341-817341	C	upstream_gene_variant	MODIFIER	FAM87B	ENSG00000177757	Transcript	ENST00000326734.3	lncRNA	-	-	-	-	-	-	-	-	-	rs3131972	A	A/C/G/T	22	1	-	HGNC	HGNC:32236	-	-	-	2	-	-	-	-	-	-	-	-	-	-	-	-	-
#=============#

# 3 用.raw 或 .ped 文件能查每个 Sample 是否有对应等位基因
plink --bfile related_removed --recode H --extract sig_snp_ids.txt --out snp_per_sample

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
grep -v "^#" significant_SNP.vcf > output_noheader.vcf
convert2annovar.pl -format vcf4 output_noheader.vcf > input.avinput
# NOTICE: Finished reading 27 lines from VCF file
# NOTICE: A total of 27 locus in VCF file passed QC threshold, representing 27 SNPs (22 transitions and 5 transversions) and 0 indels/substitutions
# NOTICE: Finished writing 27 SNP genotypes (22 transitions and 5 transversions) and 0 indels/substitutions for 1 sample

table_annovar.pl input.avinput ../../../data/annovar_dataset/humandb -buildver hg19 -out output --nopolish -remove -protocol refGene,cytoBand -operation g,r -nastring . -vcfinput
# NOTICE: Running with system command <convert2annovar.pl  -includeinfo -allsample -withfreq -format vcf4 input.avinput > output.avinput>
# NOTICE: Finished reading 27 lines from VCF file
# NOTICE: A total of 27 locus in VCF file passed QC threshold, representing 27 SNPs (22 transitions and 5 transversions) and 0 indels/substitutions
# NOTICE: Finished writing allele frequencies based on 27 SNP genotypes (22 transitions and 5 transversions) and 0 indels/substitutions for 0 samples

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

# step7: get SNP-patient info  ====================================================================================================
./extract_snps.sh
# ============================ #
#!/bin/bash

# 读取显著 SNP ID 文件
snp_file="significant_snps_forexact.txt"

# 检查文件是否存在
if [[ ! -f "$snp_file" ]]; then
    echo "Error: File $snp_file does not exist."
    exit 1
fi

# 输出文件
output_file="patient_iids.txt"
echo "SNP_ID    Patient_IID" > "$output_file"  # 添加表头

# 遍历每个 SNP ID
while IFS= read -r snp_id; do
    # 生成提取文件的名称
    extract_file="${snp_id}_genotype"

    # 使用 PLINK 提取 SNP 基因型
    plink --bfile ../maf_filtered --extract <(echo "$snp_id") --recode --out "$extract_file"

    # 读取生成的PED文件（假设文件名为 extract_file.ped）
    ped_file="${extract_file}.ped"

    # 检查PED文件是否存在
    if [[ -f "$ped_file" ]]; then
        # 提取病人 IID，假设 IID 在 PED 文件的第二列（第二列索引为 2）
        iids=$(awk '{print $2}' "$ped_file" | tr '\n' ' ')  # 将 IID 以空格分隔

        # 将结果写入输出文件
        echo "$snp_id    $iids" >> "$output_file"
    else
        echo "File $ped_file does not exist."
    fi
done < "$snp_file"

echo "Patient IIDs have been saved to $output_file."
# ============================ #



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