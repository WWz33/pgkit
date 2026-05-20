# pgkit - 泛基因组分析工具箱

<!-- README-I18N:START -->

[English](./README.md) | **简体中文**

<!-- README-I18N:END -->

一个面向植物泛基因组流程的轻量 Python 工具箱，基于 OrthoFinder 输出。

如果你已经有 OrthoFinder 结果目录，默认先用 `pgkit run` 跑 PAV + 曲线 + 可视化。若还要看分组差异，再加 `-m metadata.tsv`。

## 你可以做什么

| 你的目标 | 从哪里开始 | 结果 |
| --- | --- | --- |
| 把 OrthoFinder 输出转成 PAV 和基因家族分类 | `pgkit pav Orthogroups/ -o results` | `pav_matrix.tsv`、频率表、分类列表、图形 |
| 跑默认端到端流程 | `pgkit run Orthogroups/ -o results` | PAV、饱和曲线、可视化、可选分组汇总 |
| 比较野生/栽培或其他样本组 | `pgkit group results/pav_matrix.tsv metadata.tsv -f results/frequency_table.tsv -o results` | 分组频率、组特异 orthogroup、配对 Fisher 检验 |
| 汇总分类数量和分布 | `pgkit stats results/frequency_table.tsv -g results/gene_count_matrix.tsv -o results` | 分类统计和每个物种的汇总 |
| 看 PAV 结构 | `pgkit heatmap results/pav_matrix.tsv -f results/frequency_table.tsv -o results` | 带行注释的热图 |
| 计算 orthogroup 或 AXT 的 Ka/Ks | `pgkit kaks Orthogroups/ all.cds.fa -t 8 -m MA -k` | Ka/Ks 表格和汇总文件 |

## 安装

```bash
git clone https://github.com/WWz33/pgkit.git
cd pgkit
mamba env create -f environment.yml -n pgkit
mamba activate pgkit
```

## 快速开始

**输入: OrthoFinder 输出目录**

```
OrthoFinder_Results/
鈹溾攢鈹€ Orthogroups/
鈹?  鈹溾攢鈹€ Orthogroups.tsv                    # Main orthogroup assignments
鈹?  鈹溾攢鈹€ Orthogroups.txt                    # Orthogroup gene lists
鈹?  鈹溾攢鈹€ Orthogroups.GeneCount.tsv          # Gene count matrix
鈹?  鈹溾攢鈹€ Orthogroups_SingleCopyOrthologues.txt  # Single-copy orthogroups
鈹?  鈹斺攢鈹€ Orthogroups_UnassignedGenes.tsv    # Unassigned genes
鈹溾攢鈹€ Orthogroup_Sequences/                  # Protein FASTA files
鈹?  鈹溾攢鈹€ OG0000000.fa
鈹?  鈹溾攢鈹€ OG0000001.fa
鈹斺攢鈹€ ...
鈹溾攢鈹€ Comparative_Genomics_Statistics/
鈹溾攢鈹€ Gene_Duplication_Events/
鈹溾攢鈹€ Gene_Trees/
鈹溾攢鈹€ MultipleSequenceAlignments/
鈹斺攢鈹€ Orthologues/
```

**命令:**

**注意**: `all.cds.fa` 是所有物种全部基因 CDS 的超集。kaks.py 只会提取 Orthogroups 中存在的基因对应 CDS。

```bash
# Full pipeline (PAV + curve + visualization)
pgkit run Orthogroups/ -o results

# Or step by step:
pgkit pav Orthogroups/ -o results
pgkit curve results/pav_matrix.tsv -o results -s 100
pgkit stats results/frequency_table.tsv -g results/gene_count_matrix.tsv -o results
pgkit group results/pav_matrix.tsv metadata.tsv -f results/frequency_table.tsv -o results

# Ka/Ks (separate command)
pgkit kaks Orthogroups/ all.cds.fa -t 8 -m MA -k

# Heatmap with population annotation
pgkit heatmap results/pav_matrix.tsv -f results/frequency_table.tsv -P pop.tsv
```

## 命令

### run - Full Pipeline

Run complete analysis: PAV + saturation curve + visualization.

```bash
pgkit run <input> [options]

Positional:
  input                 OrthoFinder output directory or Orthogroups.tsv file

Options:
  -o, --output          Output directory (default: pgkit_output)
  -t, --threshold       Soft-core threshold (default: 0.9)
  -f, --format          Image format: png, pdf, svg (default: png)
  -s, --simulations     Simulations for curve (default: 100)
  -r, --save-r          Save R scripts to output directory
  -n, --no-plot         Skip visualization step
  -m, --metadata        Optional sample metadata TSV for group summaries
  --sample-col          Sample column name in metadata
  --group-col           Group column name in metadata
```

**Example:**
```bash
pgkit run Orthogroups/ -o results -f pdf -s 200
```

### pav - PAV Matrix Construction

Build PAV matrix, classify gene families, and auto-generate visualizations.

```bash
pgkit pav <input> [options]

Positional:
  input                 OrthoFinder output directory or Orthogroups.tsv file

Options:
  -o, --output          Output directory (default: pgkit_output)
  -t, --threshold       Soft-core threshold (default: 0.9)
  -f, --format          Image format: png, pdf, svg (default: png)
  -r, --save-r          Save R scripts to output directory
  -n, --no-plot         Skip visualization step
  -s, --stacked         Use stacked bar chart
```

**Example:**
```bash
pgkit pav Orthogroups/ -o results -f pdf -r -s
```

**Output:**
```
results/
鈹溾攢鈹€ pav_matrix.tsv           # PAV matrix (1/0)
鈹溾攢鈹€ frequency_table.tsv      # Frequency table
鈹溾攢鈹€ gene_category.tsv        # Gene-category-species table
鈹溾攢鈹€ gene_count_matrix.tsv    # Gene count per species per category
鈹溾攢鈹€ core_orthogroups.txt     # Core orthogroup list
鈹溾攢鈹€ soft_core_orthogroups.txt
鈹溾攢鈹€ dispensable_orthogroups.txt
鈹溾攢鈹€ private_orthogroups.txt
鈹溾攢鈹€ pgkit.pie.png            # Pie chart
鈹溾攢鈹€ pgkit.bar.png            # Bar chart
鈹溾攢鈹€ pgkit.heatmap.png        # Heatmap
鈹斺攢鈹€ r_scripts/               # R scripts (if -r)
```

### curve - Saturation Curve

Generate Core/Pan gene family saturation curve with curve fitting.

**Method**:
1. For each sample count k (1 to n), randomly sample k species
2. Calculate core (present in all k) and pan (present in at least 1) gene families
3. Repeat N times (default: 100) to get mean 卤 SD
4. Fit Heaps' law for pan-genome: `Pan = P1 * n^gamma + P2`
5. Fit exponential decay for core-genome: `Core = C1 * exp(-C2 * n) + C3`

```bash
pgkit curve <pav_matrix> [options]

Positional:
  pav_matrix            PAV matrix file (pav_matrix.tsv)

Options:
  -o, --output          Output directory (default: pgkit_output)
  -s, --simulations     Simulations per sample count (default: 100)
```

### pie - Pie Chart

Generate pie chart showing gene family category proportions.

```bash
pgkit pie <frequency_table> [options]

Positional:
  frequency_table       Frequency table file (frequency_table.tsv)

Options:
  -o, --output          Output directory (default: pgkit_output)
```

### bar - Bar Chart

Generate bar chart showing gene family composition by sample.

```bash
pgkit bar <gene_count> [options]

Positional:
  gene_count            Gene count matrix (gene_count_matrix.tsv)

Options:
  -o, --output          Output directory (default: pgkit_output)
  -s, --stacked         Use stacked bar chart
```

### heatmap - Heatmap

Generate heatmap visualization of PAV matrix.

```bash
pgkit heatmap <pav_matrix> [options]

Positional:
  pav_matrix            PAV matrix file (pav_matrix.tsv)

Options:
  -f, --frequency       Frequency table for row annotation
  -P, --pop             Population file (2-col TSV: species, group)
  -o, --output          Output directory (default: pgkit_output)
```

**Pop file format** (no header):
```
SpeciesA    Group1
SpeciesB    Group1
SpeciesC    Group2
```

**Example:**
```bash
pgkit heatmap results/pav_matrix.tsv -f results/frequency_table.tsv -o results
pgkit heatmap results/pav_matrix.tsv -f results/frequency_table.tsv -P pop.tsv -o results
```

### group - 分组 PAV 比较

比较不同样本组之间的基因家族存在频率，例如野生/栽培、亚群、地理来源或材料类型。这是轻量 PAV 频率比较，不是表型关联或 GWAS。

```bash
pgkit group <pav_matrix> <metadata> [options]

Positional:
  pav_matrix            PAV matrix file (pav_matrix.tsv)
  metadata              Sample metadata TSV with sample and group columns

Options:
  -o, --output          Output directory (default: pgkit_output)
  -f, --frequency       Frequency table for category annotation
  --sample-col          Sample column name (default: auto/first)
  --group-col           Group column name (default: auto/second)
  --specific-min        Minimum within-group frequency for group-specific calls (default: 1.0)
  --outside-max         Maximum outside-group frequency for group-specific calls (default: 0.0)
```

**Metadata format**:
```text
Sample    Group
SpeciesA  Wild
SpeciesB  Wild
SpeciesC  Cultivar
```

**Example:**
```bash
pgkit group results/pav_matrix.tsv metadata.tsv -f results/frequency_table.tsv -o results
pgkit run Orthogroups/ -o results -m metadata.tsv
```

**Output:**
```
results/
├── group_frequency.tsv              # Per-group PAV frequency per orthogroup
├── group_specific_orthogroups.tsv   # Group-specific or threshold-specific orthogroups
└── group_pairwise.tsv               # Pairwise Fisher exact tests with FDR
```

### stats - Statistics Report

Generate statistics report.

```bash
pgkit stats <frequency_table> [options]

Positional:
  frequency_table       Frequency table file (frequency_table.tsv)

Options:
  -g, --gene-count      Gene count matrix (gene_count_matrix.tsv)
  -o, --output          Output directory (default: pgkit_output)
```

### kaks - Ka/Ks Calculation (Separate Command)

Calculate Ka/Ks (non-synonymous/synonymous substitution rates). This is a separate command, not part of the `pav` workflow.

**Modes**:
1. **Standalone mode** (`-i input.axt`): Direct AXT input, equivalent to KaKs_Calculator 3.0
2. **Pan-genome mode** (orthogroups_dir + cds_file): Random sampling by gene family category

**Pan-genome mode method** (Sun et al., 2022 Nature):
```
1. Classify orthogroups into core/soft-core/dispensable/private
2. Random sample N orthogroups from each category (default: 50)
3. For each orthogroup, randomly select P species pairs (default: 50)
4. For each pair:
   - Load protein from Orthogroup_Sequences/ (consistent with clustering)
   - Extract matching CDS from all.cds.fa (by gene ID)
   - Protein alignment (MUSCLE/MAFFT)
   - Back-translate to CDS alignment
   - Calculate Ka/Ks using selected method
5. Compare Ka/Ks distributions across categories (Kruskal-Wallis test)
```

**Note**: `all.cds.fa` is a superset containing all CDS sequences for all genes across all species. kaks.py extracts only the CDS corresponding to genes present in Orthogroups.

**Note**: `kaks.py` is a Python reimplementation of [KaKs_Calculator 3.0](https://ngdc.cncb.ac.cn/biocode/tools/BT000001) (Zhang et al., 2021). It supports all methods from KaKs_Calculator 3.0 including Model Averaging (MA) and Model Selection (MS), with built-in Python fallback (Nei-Gojobori) when KaKs_Calculator is not available. Use `-k` flag to call external KaKs_Calculator 3.0 executable for maximum accuracy.

**AXT Format:**
```
seq1_name seq2_name
ATGCGTACGTAGCTAGC...
ATGCGTACGTAGCTAGC...
<blank line>
```

```bash
# Standalone mode
pgkit kaks -i pairs.axt -o kaks_output [options]

# Pan-genome mode
pgkit kaks <orthogroups_dir> <cds_file> [options]

Options:
  -i, --input           Standalone: input AXT file
  -o, --output          Output directory (default: kaks_results)
  -n, --n-genes         Pan-genome: orthogroups to sample per category (default: 50)
  -p, --n-pairs         Pan-genome: species pairs per orthogroup (default: 50)
  -t, --threads         Number of threads (default: 1)
  -s, --seed            Random seed (default: 42)
  -T, --threshold       Soft-core threshold (default: 0.9)
  -c, --genetic-code    Genetic code table 1-33 (default: 1=universal)
  -m, --method          Ka/Ks method: NG, LWL, LPB, GY, YN, MYN, MS, MA (default: MA)
  -k, --use-kaks-calculator   Call external KaKs_Calculator 3.0 (more accurate)
  -C, --calculator-path       Path to KaKs_Calculator 3.0 executable
  --check-ids           Only check CDS/protein ID matching, then exit
```

**Example:**
```bash
# Standalone mode (like KaKs_Calculator 3.0)
pgkit kaks -i pairs.axt -o kaks_output -m MA
pgkit kaks -i pairs.axt -o kaks_output -m YN -t 8

# Pan-genome mode (Python fallback)
pgkit kaks Orthogroups/ all.cds.fa -n 50 -p 50

# Pan-genome mode with KaKs_Calculator 3.0
pgkit kaks Orthogroups/ all.cds.fa -t 8 -m MA -k
```

**Output:**
```
kaks_results/
鈹溾攢鈹€ kaks_values.tsv      # All Ka/Ks values
鈹溾攢鈹€ kaks_summary.tsv     # Summary statistics by category
鈹溾攢鈹€ kaks_invalid.tsv     # Skipped sequences (invalid CDS)
鈹溾攢鈹€ kaks_boxplot.R       # R visualization script
鈹斺攢鈹€ tmp/                 # Temporary files
```

## 分类标准

| Category | Definition | Example (46 samples) |
|----------|------------|----------------------|
| **Core** | Present in 100% samples | 42 orthogroups |
| **Soft-core** | Present in >=90% samples | 20 orthogroups |
| **Dispensable** | Present in some samples | 62 orthogroups |
| **Private** | Present in single sample | 25 orthogroups |

## KaKs_Calculator 3.0 兼容方法

`kaks.py` 支持 KaKs_Calculator 3.0 (Zhang et al., 2021) 的全部方法：

| Method | Reference | Description | v3.0 |
|--------|-----------|-------------|------|
| **NG** | Nei & Gojobori (1986) | Simple, fast | 鉁?|
| **LWL** | Li, Wu & Luo (1985) | Weighted sites | 鉁?|
| **LPB** | Li (1993), Pamilo & Bianchi (1993) | Improved weighting | 鉁?|
| **GY** | Goldman & Yang (1994) | ML, codon model | 鉁?|
| **YN** | Yang & Nielsen (2000) | ML, HKY model | 鉁?|
| **MYN** | Zhang, Li & Yu (2006) | Modified YN | 鉁?|
| **MS** | Zhang et al. (2006) | Model Selection | 鉁?(v3.0) |
| **MA** | Zhang et al. (2006) | Model Averaging | 鉁?(v3.0) [DEFAULT] |

## R 可视化脚本

位于 `pgkit/src/scripts/`：

| Script | Description |
|--------|-------------|
| `plot_pie.R` | 饼图 |
| `plot_bar.R` | 柱状图 |
| `plot_heatmap.R` | 热图 (pheatmap) |
| `plot_heatmap_enhanced.R` | 热图 (ComplexHeatmap, enhanced) |
| `plot_curve.R` | 饱和曲线 |
| `plot_curve_enhanced.R` | 带拟合的饱和曲线 |
| `plot_stackbar_enhanced.R` | 带树状图的堆叠柱状图 |
| `plot_hist_ring.R` | 直方图 + 环形图 |
| `plot_halfviolin.R` | 半小提琴 + 抖动点 |

**手动运行 R 脚本:**
```bash
Rscript pgkit/src/scripts/plot_pie.R frequency_table.tsv out_prefix png
Rscript pgkit/src/scripts/plot_heatmap.R pav_matrix.tsv out_prefix frequency_table.tsv png
```

## 项目结构

```
pgkit/
鈹溾攢鈹€ pgkit.py                 # CLI entry point
鈹溾攢鈹€ lib/
鈹?  鈹溾攢鈹€ __init__.py
鈹?  鈹溾攢鈹€ parser.py            # OrthoFinder output parser
鈹?  鈹溾攢鈹€ classify.py          # Gene family classification
鈹斺攢鈹€ utils.py             # Utility functions
鈹斺攢鈹€ src/
    鈹溾攢鈹€ __init__.py
    鈹溾攢鈹€ pav.py               # PAV construction + classification
    鈹溾攢鈹€ curve.py             # Saturation curve
    鈹溾攢鈹€ pie.py               # Pie chart
    鈹溾攢鈹€ bar.py               # Bar chart
    鈹溾攢鈹€ heatmap.py           # Heatmap
    鈹溾攢鈹€ stats.py             # Statistics report
    鈹溾攢鈹€ kaks.py              # Ka/Ks calculation
    鈹斺攢鈹€ scripts/
        鈹溾攢鈹€ plot_pie.R
        鈹溾攢鈹€ plot_bar.R
        鈹溾攢鈹€ plot_heatmap.R
        鈹溾攢鈹€ plot_heatmap_enhanced.R
        鈹溾攢鈹€ plot_curve.R
        鈹溾攢鈹€ plot_curve_enhanced.R
        鈹溾攢鈹€ plot_stackbar_enhanced.R
        鈹溾攢鈹€ plot_hist_ring.R
        鈹斺攢鈹€ plot_halfviolin.R
```

## 参考文献

1. **OrthoFinder**: Emms & Kelly (2019) Genome Biology
2. **KaKs_Calculator 3.0**: Zhang et al. (2021) Genomics Proteomics Bioinformatics
3. **ParaAT**: Zhang et al. (2012) PLoS ONE
4. **Potato pan-genome**: Sun et al. (2022) Nature
5. **Grapevine pan-genome**: (2023) Bioinformatics
6. **APAVplot**: Visualization R package for PAV analysis

## 许可证

MIT

## 联系方式

如有问题或建议，请在 GitHub 提交 issue。
