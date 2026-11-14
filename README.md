# NGS Bisulfite-Seq DMR Analyzer & Genomic Annotation Suite

A command-line Python tool for comparative epigenomics. It compares a Control and a Treatment bisulfite-seq methylation track, finds Differentially Methylated Regions (DMRs), labels each one with its genomic context (CpG island, promoter or open sea), and writes a CSV report and a plot.

![Methylation change per interval](docs/images/methylation_graph.png)

## What it does

1. **Reads two BED tracks** (chromosome, start, end, methylation %) for Control and Treatment.
2. **Matches intervals by exact coordinates.** The Control track is loaded into a dictionary keyed by `(chrom, start, end)` for constant-time lookup.
3. **Quantifies the change:** `Δβ = Treatment % − Control %`. A positive value is hypermethylated, a negative value is hypomethylated.
4. **Filters by a cutoff** you set with `-d`: a region is kept if `|Δβ| >= min-diff` (default 25.0).
5. **Annotates each DMR** against a feature track using interval overlap (`start <= feature_end and end >= feature_start`). The first overlapping feature is reported; if none overlaps, the label is `Intergenic_Open_Sea`.
6. **Exports** a CSV report and a PNG plot.

![Pipeline overview](docs/images/work_explanation.png)

## Installation

Requires Python 3.8 or newer.

```bash
git clone https://github.com/<your-username>/dmr-analyzer.git
cd dmr-analyzer
pip install -r requirements.txt
```

## Usage

Run with the built-in test data (the script creates the test tracks and a small annotation file if they are not found):

```bash
python cpg_pipeline.py
```

On Windows you can also use `py cpg_pipeline.py`.

Run on the example files in this repository:

```bash
python cpg_pipeline.py -c examples/data/control_methylation.bed -t examples/data/treatment_methylation.bed -a examples/data/genomic_features.bed
```

Change the cutoff and output names:

```bash
python cpg_pipeline.py -d 40.0 -o strict_report.csv -p strict_plot.png
```

### Command-line arguments

| Flag | Default | Description |
|---|---|---|
| `-c`, `--control` | `control_methylation.bed` | Control BED track |
| `-t`, `--treatment` | `treatment_methylation.bed` | Treatment BED track |
| `-a`, `--annotation` | `genomic_features.bed` | Reference feature BED track (4th column is the label) |
| `-d`, `--min-diff` | `25.0` | Minimum absolute methylation difference (%) to call a DMR |
| `-o`, `--output` | `dmr_report.csv` | CSV report path |
| `-p`, `--plot` | `volcano_plot.png` | Plot image path |

### Input format

Tab- or space-separated BED-style lines:

```
chr1    1000    1200    30.0
```

Control and Treatment files must contain the same intervals (exact coordinate match).

### Output format

`dmr_report.csv` has one row per retained DMR:

`Chromosome, Start, End, Control_Methylation_%, Treatment_Methylation_%, Delta_Beta_%, Status, Genomic_Context`

## Results on the example data

The example tracks contain three intervals on chr1.

| Interval | Control % | Treatment % | Δβ | Status | Context |
|---|---|---|---|---|---|
| chr1:1000-1200 | 30.0 | 85.0 | +55.0 | Hypermethylated | CpG_Island |
| chr1:1200-1400 | 85.0 | 88.0 | +3.0 | filtered out (noise) | n/a |
| chr1:1400-1600 | 40.0 | 10.0 | −30.0 | Hypomethylated | Promoter_TSS |

| Cutoff (`-d`) | DMRs retained |
|---|---|
| 25.0 (default) | 2 |
| 30.0 | 2 |
| 40.0 | 1 (+55.0 only) |

Example outputs are in [`examples/output/`](examples/output/).

![Plot output](examples/output/volcano_plot.png)

## Scope and known limitations

- Validated only on small synthetic tracks (3 intervals), not on real sequencing data.
- Intervals are matched by exact coordinates only; overlap-based matching is not implemented.
- The plot's y-axis is |Δβ|, and its dash-dot boundary curves are a fixed visual guide. They are not a statistical significance test and do not change with `-d`.
- No statistical testing or read-coverage weighting yet.
- If an input file is missing, the script currently exits without printing an error message.

## Roadmap

- Overlap-based interval matching
- Statistical testing (for example Fisher's exact test with multiple-testing correction) and coverage weighting
- A clear error message for missing input files
- Testing on real bisulfite-seq data

## Author

Garv Gupta (Roll No: 25BIB013)
