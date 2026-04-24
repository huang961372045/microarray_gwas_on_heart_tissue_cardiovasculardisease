# Microarray GWAS Analysis on HTX Heart Tissue Samples
# 心脏组织 HTX 样本微阵列全基因组关联分析

[![PLINK](https://img.shields.io/badge/PLINK-v1.90b7.2-blue.svg)](https://www.cog-genomics.org/plink/1.9/)
[![Genome](https://img.shields.io/badge/Genome-GRCh38/hg38-green.svg)](https://www.ncbi.nlm.nih.gov/grc/human)
[![Array](https://img.shields.io/badge/Array-Illumina%20GSAv4-orange.svg)](https://www.illumina.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> A PLINK-based GWAS pipeline applied to 370 human heart tissue samples genotyped on the Illumina Infinium Global Screening Array v4 (GSAv4, hg38), followed by Ensembl VEP functional annotation.
>
> 基于 PLINK 的 GWAS 分析流程，对 370 例人心脏组织样本（Illumina GSAv4 芯片，hg38）进行质控、关联分析与 Ensembl VEP 功能注释。

---

## 📑 Table of Contents / 目录

- [Background / 背景](#background--背景)
- [Dataset / 数据集](#dataset--数据集)
- [Pipeline Overview / 流程概览](#pipeline-overview--流程概览)
- [Repository Structure / 仓库结构](#repository-structure--仓库结构)
- [Requirements / 依赖环境](#requirements--依赖环境)
- [Reproducibility / 复现步骤](#reproducibility--复现步骤)
- [Key Results / 主要结果](#key-results--主要结果)
- [Data Availability / 数据可用性](#data-availability--数据可用性)
- [Citation / 引用](#citation--引用)
- [Contact / 联系方式](#contact--联系方式)

---

## Background / 背景

**EN** — Heart transplantation (HTX) samples offer a rare opportunity to study the genetic architecture of end-stage cardiac disease in tissue-of-origin DNA. This project performs a genome-wide association study (GWAS) on microarray-genotyped HTX samples to identify SNPs associated with cardiomyopathy phenotypes, with a focus on dilated cardiomyopathy (DCM) versus control comparisons.

**中文** — 心脏移植（HTX）样本为研究终末期心脏疾病的遗传结构提供了宝贵的组织源 DNA。本项目对微阵列基因分型的 HTX 样本进行全基因组关联分析，识别与心肌病相关的 SNP 位点，重点关注扩张型心肌病（DCM）与对照组之间的比较。

### Phenotypes Included / 涉及表型

| Code | Full name (EN) | 中文名 |
|------|---------------|--------|
| DCM  | Dilated Cardiomyopathy | 扩张型心肌病 |
| ICM  | Ischemic Cardiomyopathy | 缺血性心肌病 |
| HCM  | Hypertrophic Cardiomyopathy | 肥厚型心肌病 |
| HLTX | Heart-Lung Transplant | 心肺联合移植 |
| Control | Non-diseased donor | 健康供体对照 |

---

## Dataset / 数据集

| Item | Value |
|------|-------|
| Platform | Illumina Infinium Global Screening Array v4 (GSAv4) |
| Reference genome | GRCh38 / hg38 (forward strand) |
| Raw variants | 696,873 |
| Raw samples | 370 (266 M / 88 F / 16 ambiguous sex) |
| Post-QC variants | 274,416 |
| Post-QC samples | 363 |
| Input files | `M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward.{ped,map}` |

> ⚠️ **Raw genotype and phenotype data contain personal health information and are NOT included in this repository.** See [Data Availability](#data-availability--数据可用性) below.
>
> ⚠️ **原始基因型与表型数据含个人健康信息，未包含在本仓库内。** 请见下方[数据可用性](#data-availability--数据可用性)。

---

## Pipeline Overview / 流程概览

```
┌──────────────────────┐    ┌──────────────────────┐    ┌──────────────────────┐
│  Raw .ped / .map     │───▶│  PLINK binary .bed   │───▶│      QC filters      │
│  696,873 variants    │    │      .bim .fam       │    │ --geno --mind --hwe  │
│     370 samples      │    │                      │    │       --maf          │
└──────────────────────┘    └──────────────────────┘    └──────────┬───────────┘
                                                                   │
┌──────────────────────┐    ┌──────────────────────┐    ┌──────────▼───────────┐
│   Ensembl VEP        │◀───│   Significant SNPs   │◀───│  Association test    │
│     annotation       │    │    (p < 1e-4)        │    │  --assoc / --linear  │
└──────────────────────┘    └──────────────────────┘    └──────────────────────┘
```

### Pipeline steps / 流程步骤

1. **Format conversion / 格式转换** — `.ped/.map` → PLINK binary `.bed/.bim/.fam`
2. **Quality control / 质控**
   - `--geno 0.05` — remove SNPs with > 5% missingness / 去除缺失率 > 5% 的 SNP
   - `--mind 0.1` — remove samples with > 10% missingness / 去除缺失率 > 10% 的样本
   - `--hwe 1e-4` — remove SNPs failing Hardy–Weinberg equilibrium / 去除不符合哈迪-温伯格平衡的 SNP
   - `--maf 0.05` — remove SNPs with MAF < 5% / 去除次要等位基因频率 < 5% 的 SNP
3. **Association analysis / 关联分析** — `--assoc` or `--linear` with phenotype and covariates
4. **Visualization / 可视化** — Manhattan plot, QQ plot (R scripts, see `scripts/`)
5. **Annotation / 注释** — VCF export → Ensembl VEP web tool → functional annotation

---

## Repository Structure / 仓库结构

```
microarray-gwas-htx/
├── README.md                       # This file / 本文件
├── LICENSE
├── .gitignore                      # Ignore PLINK intermediate files
├── .gitattributes                  # Git LFS tracking rules
│
├── scripts/                        # Analysis scripts / 分析脚本
│   ├── 01_run_plink.sh             # Base QC pipeline, no phenotype
│   ├── 02_run_plink_v2.sh          # With phenotype file
│   └── 03_run_plink_v3_DCM.sh      # DCM vs Control comparison
│
├── data/
│   ├── phenotype/
│   │   └── Phenotype_Information_HTX_samples.csv   # ⚠️ anonymized version only
│   └── README.md                   # Where to request raw genotype data
│
├── results/
│   ├── gwas/
│   │   ├── significant_snps.csv    # 27 top SNPs (p < 1e-4)
│   │   └── sig_snp_ids.txt         # Full SNP ID list (LFS)
│   ├── per_sample/
│   │   └── snp_per_sample.raw      # Per-sample genotypes (LFS)
│   └── annotation/
│       ├── vep_output.txt          # VEP tab output (LFS)
│       └── vep_output.vcf          # VEP annotated VCF (LFS)
│
├── docs/
│   ├── 20240516_microarray.pptx    # Progress presentation
│   ├── conclusion_microarray_DNA.xlsx
│   └── plink_help.txt              # PLINK 1.9 full help reference
│
└── figures/                        # Manhattan / QQ plots (to add)
```

---

## Requirements / 依赖环境

### Software / 软件

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | 1.90b7.2 (Dec 2023) | Genotype QC and association analysis |
| R | ≥ 4.2 | Manhattan / QQ plot visualization |
| Ensembl VEP | Web tool (current release) | Variant effect prediction |
| Git LFS | ≥ 3.0 | Large file version control |

### R packages / R 包

```r
install.packages(c("qqman", "ggplot2", "data.table", "dplyr"))
```

### Git LFS setup / Git LFS 配置

Because some result files exceed GitHub's 100 MB limit, this repository uses Git LFS.
由于部分结果文件超过 GitHub 100 MB 限制，本仓库使用 Git LFS。

```bash
# Install Git LFS (first time only)
git lfs install

# After cloning:
git lfs pull
```

Tracked file patterns (see `.gitattributes`):
```
*.raw    filter=lfs diff=lfs merge=lfs -text
*.vcf    filter=lfs diff=lfs merge=lfs -text
results/annotation/*.txt  filter=lfs diff=lfs merge=lfs -text
```

---

## Reproducibility / 复现步骤

### 1. Clone the repository / 克隆仓库

```bash
git clone https://github.com/huang961372045/microarray-gwas-htx.git
cd microarray-gwas-htx
git lfs pull
```

### 2. Place raw genotype data / 放置原始基因型数据

Place `M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward.{ped,map}` in `data/raw/` (this directory is gitignored).
将 `M01144_01-04_PLINK_genotypes_GSAv4_hg38_forward.{ped,map}` 放入 `data/raw/` 目录（此目录已被 gitignore）。

### 3. Run the pipeline / 运行流程

```bash
# Base QC pipeline
bash scripts/01_run_plink.sh

# With phenotype
bash scripts/02_run_plink_v2.sh

# DCM vs Control analysis
bash scripts/03_run_plink_v3_DCM.sh
```

### 4. Annotate significant variants / 注释显著变异

Upload `results/gwas/significant_snps.vcf` to [Ensembl VEP Web](https://www.ensembl.org/Tools/VEP) and save outputs to `results/annotation/`.
上传 VCF 到 Ensembl VEP 网页工具，结果保存到 `results/annotation/`。

---

## Key Results / 主要结果

### Post-QC summary / 质控后概要

| Stage | Variants | Samples |
|-------|----------|---------|
| Raw input | 696,873 | 370 |
| After `--geno` / `--mind` | 651,477 | 363 |
| After `--hwe` | 650,890 | 363 |
| After `--maf` | **274,416** | **363** |

### Top significant SNPs / 显著 SNP 位点

27 SNPs passed the suggestive significance threshold (p < 1e-4). The top hits are:
27 个 SNP 位点通过提示性显著性阈值（p < 1e-4），其中顶端结果为：

| CHR | SNP | BP | A1 | BETA | STAT | P-value |
|-----|-----|----|----|------|------|---------|
| X (23) | rs4986503  | 17,369,377  | C | 0.5142 | 5.034 | 7.65e-07 |
| 12     | rs1211610  | 131,122,509 | G | 0.7970 | 5.006 | 8.77e-07 |
| 9      | GSA-rs943170 | 7,468,664 | C | 0.7873 | 4.542 | 7.65e-06 |
| 8      | rs2721266  | 2,895,605   | T | 0.3925 | 4.512 | 8.73e-06 |
| 9      | GSA-rs17256298 | 80,489,443 | A | 0.6787 | 4.324 | 1.99e-05 |

> Full table: [`results/gwas/significant_snps.csv`](results/gwas/significant_snps.csv)

### Notes on chromosome X / 关于 X 染色体说明

PLINK reported warnings about heterozygous haploid genotypes on chromosome X and non-missing non-male Y genotypes. These are expected for GSA arrays and are treated as missing by downstream commands. For formal X-chromosome analysis, consider a less stringent HWE threshold on chromosome X as recommended by PLINK.
PLINK 在 X 染色体上报告了杂合单倍型基因型警告。这在 GSA 芯片中是常见现象，下游分析会将其视为缺失值。若正式做 X 染色体分析，建议对 X 染色体单独设置较宽松的 HWE 阈值。

---

## Data Availability / 数据可用性

| Data | Access |
|------|--------|
| Analysis scripts | This repository, MIT license / 本仓库，MIT 许可证 |
| Summary statistics (top SNPs) | This repository / 本仓库 |
| VEP annotations | This repository (Git LFS) / 本仓库（Git LFS） |
| Raw genotype data (`.ped/.map`) | Available upon reasonable request to the corresponding author, subject to data use agreement. / 根据合理请求向通讯作者申请，需签署数据使用协议。 |
| Sample phenotype / clinical info | Restricted — contact corresponding author / 受限 — 请联系通讯作者 |

---

## Citation / 引用

If you use the scripts or results from this repository, please cite:
如使用本仓库的脚本或结果，请引用：

```bibtex
@misc{htx_microarray_gwas_2024,
  author       = {Ziming Huang},
  title        = {Microarray GWAS Analysis on HTX Heart Tissue Samples},
  year         = {2024},
  howpublished = {\url{https://github.com/huang961372045/microarray-gwas-htx}}
}
```

And the underlying tools:
以及所使用的工具：

- **PLINK 1.9** — Chang CC, et al. *GigaScience* 2015;4:7. doi:10.1186/s13742-015-0047-8
- **Ensembl VEP** — McLaren W, et al. *Genome Biology* 2016;17:122. doi:10.1186/s13059-016-0974-4

---

## Contact / 联系方式

| Role | Name | Contact |
|------|------|---------|
| Maintainer / 维护者 | Ziming | <huangziming222@gmail.com> |
| Principal Investigator | *to fill* | — |
| Institution | *to fill* | — |

---

## License / 许可证

Code in this repository is released under the MIT License (see [`LICENSE`](LICENSE)). Derived genomic summary statistics are released under CC-BY-4.0. Raw patient-level data are **not** released and remain subject to the original ethics approval and data protection regulations (GDPR).
本仓库代码采用 MIT 许可证发布（详见 [`LICENSE`](LICENSE)）。衍生基因组汇总统计数据采用 CC-BY-4.0 发布。原始病人级数据**不公开**，受原伦理审批及数据保护法规（GDPR）约束。
