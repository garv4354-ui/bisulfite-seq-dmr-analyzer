import os
import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def generate_test_bed_files(control_file, treatment_file):
    if not os.path.exists(control_file):
        with open(control_file, 'w') as f:
            f.write("chr1\t1000\t1200\t30.0\n")
            f.write("chr1\t1200\t1400\t85.0\n")
            f.write("chr1\t1400\t1600\t40.0\n")
    if not os.path.exists(treatment_file):
        with open(treatment_file, 'w') as f:
            f.write("chr1\t1000\t1200\t85.0\n")
            f.write("chr1\t1200\t1400\t88.0\n")
            f.write("chr1\t1400\t1600\t10.0\n")

def generate_annotation_file(annotation_path='genomic_features.bed'):
    if not os.path.exists(annotation_path):
        with open(annotation_path, 'w') as f:
            f.write("chr1\t950\t1250\tCpG_Island\n")
            f.write("chr1\t1350\t1450\tPromoter_TSS\n")
            f.write("chr1\t1550\t1650\tOpen_Sea\n")

def load_annotations(annotation_path='genomic_features.bed'):
    features = []
    if os.path.exists(annotation_path):
        with open(annotation_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    features.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return features

def annotate_dmr(chrom, start, end, features):
    for feat_chrom, feat_start, feat_end, label in features:
        if chrom == feat_chrom:
            if start <= feat_end and end >= feat_start:
                return label
    return 'Intergenic_Open_Sea'

def generate_volcano_plot(dmrs, plot_path='volcano_plot.png'):
    if not dmrs:
        print('-> No DMRs found to plot.')
        return
    deltas = [d[3] for d in dmrs]
    magnitudes = [abs(d[3]) for d in dmrs]
    statuses = [d[4] for d in dmrs]
    contexts = [d[5] for d in dmrs]
    coords = [f"{d[0][0]}:{d[0][1]}-{d[0][2]}" for d in dmrs]

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = ['#e74c3c' if s == 'Hypermethylated' else '#3498db' for s in statuses]
    
    ax.scatter(deltas, magnitudes, c=colors, alpha=0.9, edgecolors='#2c3e50', s=160, zorder=3)
    ax.axvline(0, color='#7f8c8d', linestyle='--', linewidth=1.2, zorder=2)

    x_pos = np.linspace(25, 65, 100)
    y_pos = (200 / (x_pos - 20)) + 24
    ax.plot(x_pos, y_pos, color='#e74c3c', linestyle='-.', linewidth=1.5, alpha=0.6, zorder=2, label='Significance Threshold Curve')
    
    x_neg = np.linspace(-65, -25, 100)
    y_neg = (200 / (np.abs(x_neg) - 20)) + 24
    ax.plot(x_neg, y_neg, color='#3498db', linestyle='-.', linewidth=1.5, alpha=0.6, zorder=2)

    for i, txt in enumerate(coords):
        ax.annotate(f"{txt}\n({contexts[i]})",
                    (deltas[i], magnitudes[i]),
                    xytext=(12, 12 if deltas[i] > 0 else -18), textcoords='offset points',
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#bdc3c7", alpha=0.95),
                    fontsize=9, fontweight='bold', color='#2c3e50')

    ax.set_title('NGS Bisulfite-Seq DMR Analyzer: Methylation Delta vs Magnitude', fontsize=13, fontweight='bold', pad=15)
    ax.set_xlabel('Methylation Difference (Delta Beta %)', fontsize=11, fontweight='semibold')
    ax.set_ylabel('Magnitude (|Delta %|)', fontsize=11, fontweight='semibold')
    ax.grid(True, linestyle=':', alpha=0.6, zorder=1)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Hypermethylated (CpG_Island)', markerfacecolor='#e74c3c', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Hypomethylated (Promoter_TSS)', markerfacecolor='#3498db', markersize=10),
        Line2D([0], [0], color='#e74c3c', linestyle='-.', lw=1.5, label='Significance Threshold Curve')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, facecolor='#f8f9fa', edgecolor='#bdc3c7')
    
    fig.text(0.99, 0.01, 'Prepared by: Garv Gupta | Roll No: 25BIB013 | BTech Bioinformatics', ha='right', va='bottom', fontsize=8, color='#7f8c8d', style='italic')

    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'-> Exported curved-threshold volcano plot successfully to: {plot_path}')

def compare_dmrs(control_path, treatment_path, annotation_path='genomic_features.bed', min_diff=25.0, output_csv=None, output_plot='volcano_plot.png'):
    if not os.path.exists(control_path) or not os.path.exists(treatment_path):
        return
    features = load_annotations(annotation_path)
    print(f'\nComparing Control ({control_path}) vs. Treatment ({treatment_path})...')
    control_data = {}
    with open(control_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                control_data[(parts[0], int(parts[1]), int(parts[2]))] = float(parts[3])
    dmrs = []
    with open(treatment_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            key = (parts[0], int(parts[1]), int(parts[2]))
            treat_meth = float(parts[3])
            if key in control_data:
                ctrl_meth = control_data[key]
                delta = treat_meth - ctrl_meth
                if abs(delta) >= min_diff:
                    status = 'Hypermethylated' if delta > 0 else 'Hypomethylated'
                    context = annotate_dmr(key[0], key[1], key[2], features)
                    dmrs.append((key, ctrl_meth, treat_meth, delta, status, context))
    print(f'-> Found {len(dmrs)} significant DMRs (|Delta| >= {min_diff}%):')
    for idx, (coord, ctrl, treat, delta, status, context) in enumerate(dmrs, 1):
        print(f'   [{idx}] {coord[0]}:{coord[1]}-{coord[2]} | Ctrl: {ctrl}% -> Treat: {treat}% | Delta: {delta:+.1f}% ({status}) | Context: {context}')
    if output_csv:
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Chromosome', 'Start', 'End', 'Control_Methylation_%', 'Treatment_Methylation_%', 'Delta_Beta_%', 'Status', 'Genomic_Context'])
            for coord, ctrl, treat, delta, status, context in dmrs:
                writer.writerow([coord[0], coord[1], coord[2], ctrl, treat, round(delta, 2), status, context])
        print(f'-> Exported annotated DMR report successfully to: {output_csv}')
    generate_volcano_plot(dmrs, output_plot)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--control', default='control_methylation.bed')
    parser.add_argument('-t', '--treatment', default='treatment_methylation.bed')
    parser.add_argument('-a', '--annotation', default='genomic_features.bed')
    parser.add_argument('-d', '--min-diff', type=float, default=25.0)
    parser.add_argument('-o', '--output', default='dmr_report.csv')
    parser.add_argument('-p', '--plot', default='volcano_plot.png')
    args = parser.parse_args()
    if args.control == 'control_methylation.bed' and args.treatment == 'treatment_methylation.bed':
        generate_test_bed_files(args.control, args.treatment)
    generate_annotation_file(args.annotation)
    compare_dmrs(args.control, args.treatment, args.annotation, args.min_diff, args.output, args.plot)
