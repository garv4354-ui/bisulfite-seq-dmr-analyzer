NGS Bisulfite-Seq DMR Analyzer & Genomic Annotation Suite is a command-line bioinformatics tool, written from scratch in Python, that identifies Differentially Methylated Regions (DMRs) between two experimental conditions.

The problem: DNA methylation changes can switch genes on or off without altering the sequence. Finding where methylation differs between a Control and a Treatment sample, and deciding which differences are large enough to matter, is a core step in epigenomic analysis.

What the tool does: it reads Control and Treatment bisulfite-seq tracks in BED format, matches regions by genomic coordinates, and calculates the methylation shift for each one (Δβ = Treatment % − Control %). Regions whose shift passes a cutoff set from the command line (-d) are called DMRs and labelled hypermethylated or hypomethylated. Each DMR is then checked against a reference feature track by interval overlap and tagged as CpG island, promoter, or open sea. The tool writes a structured CSV report and a plot.

Results on test data: on a three-interval validation set, it identified a +55.0% hypermethylated region in a CpG island and a −30.0% hypomethylated region at a promoter, and filtered out a +3.0% change as background noise. Raising the cutoff from 25 to 40 reduced the output from two DMRs to one, confirming the threshold works as intended.

How it was built: it grew in stages, from a single-track filter to paired comparison, a six-flag argparse interface, input handling, CSV reporting, annotation, and plotting. It uses Python, argparse, csv, numpy and matplotlib.

Scope: so far it has only been validated on small synthetic tracks. Planned next steps are overlap-based interval matching, statistical testing with multiple-testing correction, coverage weighting, and testing on real bisulfite-seq data.

