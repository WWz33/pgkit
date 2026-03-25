# pgkit - Pan-gene Family Analysis Toolkit

A comprehensive Python toolkit for pan-genome analysis based on OrthoFinder output.

## Features

- **PAV Matrix Construction**: Auto-detect OrthoFinder output, merge UnassignedGenes
- **Gene Family Classification**: Core, Soft-core, Dispensable, Private
- **Saturation Curve**: Core/Pan gene family growth curve
- **Visualization**: Pie, Bar, Heatmap, Saturation Curve (R scripts included)
- **Ka/Ks Calculation**: Selection pressure analysis (standalone + pan-genome mode)
- **Statistics Report**: Comprehensive summary

## Installation

```bash
git clone https://github.com/WWz33/pgkit.git
cd pgkit
mamba env create -f environment.yml -n pgkit
mamba activate pgkit
```

## Quick Start

**Input: OrthoFinder output directory**

```
Orthogroups/
├── Orthogroups.tsv                    # Main orthogroup assignments
├── Orthogroups.txt                    # Orthogroup gene lists
├── Orthogroups.GeneCount.tsv          # Gene count matrix
├── Orthogroups_SingleCopyOrthologues.txt  # Single-copy orthogroups
├── Orthogroups_UnassignedGenes.tsv    # Unassigned genes
└── Orthogroup_Sequences/              # Protein FASTA files (for Ka/Ks)
    ├── OG0000000.fa
    ├── OG0000001.fa
    └── ...
```

**Commands:**

**Note**: `all.cds.fa` is a superset containing all CDS sequences for all genes across all species. kaks.py extracts only the CDS corresponding to genes present in Orthogroups.

```bash
# 1. Build PAV matrix + classification + auto visualization
python pgkit/pgkit.py pav Orthogroups/ -o results

# 2. Generate saturation curve
python pgkit/pgkit.py curve results/pav_matrix.tsv -o results -s 100

# 3. Generate statistics report
python pgkit/pgkit.py stats results/frequency_table.tsv -g results/gene_count_matrix.tsv -o results

# 4. Calculate Ka/Ks
python pgkit/src/kaks.py Orthogroups/ all.cds.fa -t 8 -m MA -k
```

## Commands

### pav - PAV Matrix Construction

Build PAV matrix, classify gene families, and auto-generate visualizations.

```bash
python pgkit/pgkit.py pav <input> [options]

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
python pgkit/pgkit.py pav Orthogroups/ -o results -f pdf -r -s
```

**Output:**
```
results/
├── pav_matrix.tsv           # PAV matrix (1/0)
├── frequency_table.tsv      # Frequency table
├── gene_category.tsv        # Gene-category-species table
├── gene_count_matrix.tsv    # Gene count per species per category
├── core_orthogroups.txt     # Core orthogroup list
├── soft_core_orthogroups.txt
├── dispensable_orthogroups.txt
├── private_orthogroups.txt
├── pgkit.pie.png            # Pie chart
├── pgkit.bar.png            # Bar chart
├── pgkit.heatmap.png        # Heatmap
└── r_scripts/               # R scripts (if -r)
```

### curve - Saturation Curve

Generate Core/Pan gene family saturation curve.

```bash
python pgkit/pgkit.py curve <pav_matrix> [options]

Positional:
  pav_matrix            PAV matrix file (pav_matrix.tsv)

Options:
  -o, --output          Output directory (default: pgkit_output)
  -s, --simulations     Simulations per sample count (default: 100)
```

**Example:**
```bash
python pgkit/pgkit.py curve results/pav_matrix.tsv -o results -s 200
```

### pie - Pie Chart

Generate pie chart showing gene family category proportions.

```bash
python pgkit/pgkit.py pie <frequency_table> [options]

Positional:
  frequency_table       Frequency table file (frequency_table.tsv)

Options:
  -o, --output          Output directory (default: pgkit_output)
```

### bar - Bar Chart

Generate bar chart showing gene family composition by sample.

```bash
python pgkit/pgkit.py bar <gene_count> [options]

Positional:
  gene_count            Gene count matrix (gene_count_matrix.tsv)

Options:
  -o, --output          Output directory (default: pgkit_output)
  -s, --stacked         Use stacked bar chart
```

### heatmap - Heatmap

Generate heatmap visualization of PAV matrix.

```bash
python pgkit/pgkit.py heatmap <pav_matrix> [options]

Positional:
  pav_matrix            PAV matrix file (pav_matrix.tsv)

Options:
  -f, --frequency       Frequency table for row annotation (optional)
  -o, --output          Output directory (default: pgkit_output)
```

### stats - Statistics Report

Generate statistics report.

```bash
python pgkit/pgkit.py stats <frequency_table> [options]

Positional:
  frequency_table       Frequency table file (frequency_table.tsv)

Options:
  -g, --gene-count      Gene count matrix (gene_count_matrix.tsv)
  -o, --output          Output directory (default: pgkit_output)
```

### kaks - Ka/Ks Calculation

Calculate Ka/Ks (non-synonymous/synonymous substitution rates). Supports two modes:

1. **Standalone mode** (`-i input.axt`): Direct AXT input, equivalent to KaKs_Calculator 3.0
2. **Pan-genome mode** (orthogroups_dir + cds_file): Random sampling by gene family category

**Pan-genome mode workflow:**
```
1. Load protein from Orthogroup_Sequences/ (consistent with OrthoFinder clustering)
2. Extract matching CDS from all.cds.fa (by gene ID)
3. Protein alignment (MUSCLE/MAFFT)
4. Back-translate to CDS alignment
5. Calculate Ka/Ks
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
python pgkit/src/kaks.py -i pairs.axt -o kaks_output [options]

# Pan-genome mode
python pgkit/src/kaks.py <orthogroups_dir> <cds_file> [options]

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
python pgkit/src/kaks.py -i pairs.axt -o kaks_output -m MA
python pgkit/src/kaks.py -i pairs.axt -o kaks_output -m YN -t 8

# Pan-genome mode (Python fallback)
python pgkit/src/kaks.py Orthogroups/ all.cds.fa -n 50 -p 50

# Pan-genome mode with KaKs_Calculator 3.0
python pgkit/src/kaks.py Orthogroups/ all.cds.fa -t 8 -m MA -k
```

**Output:**
```
kaks_results/
├── kaks_values.tsv      # All Ka/Ks values
├── kaks_summary.tsv     # Summary statistics by category
├── kaks_invalid.tsv     # Skipped sequences (invalid CDS)
├── kaks_boxplot.R       # R visualization script
└── tmp/                 # Temporary files
```

## Classification Criteria

| Category | Definition | Example (46 samples) |
|----------|------------|----------------------|
| **Core** | Present in 100% samples | 42 orthogroups |
| **Soft-core** | Present in >=90% samples | 20 orthogroups |
| **Dispensable** | Present in some samples | 62 orthogroups |
| **Private** | Present in single sample | 25 orthogroups |

## KaKs_Calculator 3.0 Compatible Methods

`kaks.py` supports all methods from KaKs_Calculator 3.0 (Zhang et al., 2021):

| Method | Reference | Description | v3.0 |
|--------|-----------|-------------|------|
| **NG** | Nei & Gojobori (1986) | Simple, fast | ✓ |
| **LWL** | Li, Wu & Luo (1985) | Weighted sites | ✓ |
| **LPB** | Li (1993), Pamilo & Bianchi (1993) | Improved weighting | ✓ |
| **GY** | Goldman & Yang (1994) | ML, codon model | ✓ |
| **YN** | Yang & Nielsen (2000) | ML, HKY model | ✓ |
| **MYN** | Zhang, Li & Yu (2006) | Modified YN | ✓ |
| **MS** | Zhang et al. (2006) | Model Selection | ✓ (v3.0) |
| **MA** | Zhang et al. (2006) | Model Averaging | ✓ (v3.0) [DEFAULT] |

## R Visualization Scripts

Located in `pgkit/src/scripts/`:

| Script | Description |
|--------|-------------|
| `plot_pie.R` | Pie chart |
| `plot_bar.R` | Bar chart |
| `plot_heatmap.R` | Heatmap (pheatmap) |
| `plot_heatmap_enhanced.R` | Heatmap (ComplexHeatmap, enhanced) |
| `plot_curve.R` | Saturation curve |
| `plot_curve_enhanced.R` | Saturation curve with fitting |
| `plot_stackbar_enhanced.R` | Stacked bar with dendrogram |
| `plot_hist_ring.R` | Histogram + ring chart |
| `plot_halfviolin.R` | Half-violin + jitter |

**Run R scripts manually:**
```bash
Rscript pgkit/src/scripts/plot_pie.R frequency_table.tsv out_prefix png
Rscript pgkit/src/scripts/plot_heatmap.R pav_matrix.tsv out_prefix frequency_table.tsv png
```

## Project Structure

```
pgkit/
├── pgkit.py                 # CLI entry point
├── lib/
│   ├── __init__.py
│   ├── parser.py            # OrthoFinder output parser
│   ├── classify.py          # Gene family classification
│   └── utils.py             # Utility functions
└── src/
    ├── __init__.py
    ├── pav.py               # PAV construction + classification
    ├── curve.py             # Saturation curve
    ├── pie.py               # Pie chart
    ├── bar.py               # Bar chart
    ├── heatmap.py           # Heatmap
    ├── stats.py             # Statistics report
    ├── kaks.py              # Ka/Ks calculation
    └── scripts/
        ├── plot_pie.R
        ├── plot_bar.R
        ├── plot_heatmap.R
        ├── plot_heatmap_enhanced.R
        ├── plot_curve.R
        ├── plot_curve_enhanced.R
        ├── plot_stackbar_enhanced.R
        ├── plot_hist_ring.R
        └── plot_halfviolin.R
```

## References

1. **OrthoFinder**: Emms & Kelly (2019) Genome Biology
2. **KaKs_Calculator 3.0**: Zhang et al. (2021) Genomics Proteomics Bioinformatics
3. **ParaAT**: Zhang et al. (2012) PLoS ONE
4. **Potato pan-genome**: Sun et al. (2022) Nature
5. **Grapevine pan-genome**: (2023) Bioinformatics
6. **APAVplot**: Visualization R package for PAV analysis

## License

MIT

## Contact

For questions or issues, please open an issue on GitHub.
