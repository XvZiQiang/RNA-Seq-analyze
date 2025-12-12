#!/opt/conda/bin/python

'''
Calculate the RNA-seq reads coverage over gene body.

Note:
1) Only input sorted and indexed BAM file(s). SAM format is not supported.
2) Genes/transcripts with mRNA length < 100 will be skipped (Number specified to "-l" cannot be < 100).
'''

import os
import sys

if sys.version_info[0] != 3:
    print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(
        sys.version_info[1]) + " This version of RSeQC needs python3!\n", file=sys.stderr)
    sys.exit()

import re
from optparse import OptionParser
import collections
from time import strftime
from os.path import basename

from numpy import std, mean
import pysam
from multiprocessing import Pool, cpu_count

__author__ = "NAN"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "NAN"
__email__ = "NAN"
__status__ = "FAKE"


def extract_sample_name(filename):
    '''Extract the sample name from the filename, removing extensions like .sorted.bam, .bam, etc.'''
    return re.sub(r'(\.sorted)?\.bam$', '', filename)


def printlog(mesg):
    '''print progress into stderr'''
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    print(mesg, file=sys.stderr)


def pearson_moment_coefficient(lst):
    '''measure skewness'''
    mid_value = lst[int(len(lst) / 2)]
    sigma = std(lst, ddof=1)
    tmp = []
    for i in lst:
        tmp.append(((i - mid_value) / sigma) ** 3)
    return mean(tmp)


def percentile_list(lst):
    '''Get 100 percentile points from a list of gene coordinates'''
    if len(lst) < 100:
        return lst
    return [lst[int(len(lst) * p / 100)] if int(len(lst) * p / 100) < len(lst) else lst[-1] for p in range(1, 101)]


def genebody_percentile(refbed, mRNA_len_cut=100):
    '''
    return percentile points of gene body
    mRNA length < mRNA_len_cut will be skipped
    '''
    if refbed is None:
        print("You must specify a bed file representing gene model\n", file=sys.stderr)
        exit(0)

    g_percentiles = {}
    transcript_count = 0
    for line in open(refbed, 'r'):
        try:
            if line.startswith(('#', 'track', 'browser')): continue
            fields = line.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            tx_end = int(fields[2])
            geneName = fields[3]
            strand = fields[5]
            geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])

            exon_starts = list(map(int, fields[11].rstrip(',\n').split(',')))
            exon_starts = list(map((lambda x: x + tx_start), exon_starts))
            exon_ends = list(map(int, fields[10].rstrip(',\n').split(',')))
            exon_ends = list(map((lambda x, y: x + y), exon_starts, exon_ends))
            transcript_count += 1
        except:
            print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
            continue
        gene_all_base = []
        for st, end in zip(exon_starts, exon_ends):
            gene_all_base.extend(list(range(st + 1, end + 1)))  # 1-based coordinates on genome
        if len(gene_all_base) < mRNA_len_cut:
            continue
        g_percentiles[geneID] = (chrom, strand,
                                 percentile_list(gene_all_base))  # get 100 points from each gene's coordinates
    printlog("Total " + str(transcript_count) + ' fragment loaded')
    return g_percentiles


def genebody_coverage(bam, position_list):
    '''
    position_list is dict returned from genebody_percentile
    position is 1-based genome coordinate
    '''
    samfile = pysam.AlignmentFile(bam, "rb")
    aggreagated_cvg = collections.defaultdict(int)

    gene_finished = 0
    for chrom, strand, positions in list(position_list.values()):
        coverage = {i: 0.0 for i in positions}
        chrom_start = positions[0] - 1
        chrom_end = positions[-1]
        for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
            ref_pos = pileupcolumn.pos + 1
            if ref_pos in coverage:
                cover_read = sum(1 for pileupread in pileupcolumn.pileups
                                 if not pileupread.is_del and not pileupread.alignment.is_qcfail
                                 and not pileupread.alignment.is_secondary
                                 and not pileupread.alignment.is_unmapped
                                 and not pileupread.alignment.is_duplicate)
                coverage[ref_pos] = cover_read
        tmp = [coverage[k] for k in sorted(coverage)]
        if strand == '-':
            tmp = tmp[::-1]
        for i, cov in enumerate(tmp):
            aggreagated_cvg[i] += cov
        gene_finished += 1

        if gene_finished % 100 == 0:
            print("\t%d fragments finished\r" % (gene_finished), end=' ', file=sys.stderr)
    return aggreagated_cvg


def visualize_coverage(coverages, output_prefix, sample_names):
    '''
    Visualize the coverage data.
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    # Create a DataFrame for coverage data
    coverage_df = pd.DataFrame(coverages, index=sample_names).T
    coverage_df.index.name = 'Percentile'

    # Create a colormap
    cmap = plt.get_cmap('tab20')
    norm = Normalize(vmin=0, vmax=len(coverage_df.columns))
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Plotting the gene body coverage
    plt.figure(figsize=(12, 8))
    plt.rcParams['font.family'] = 'Arial'
    for i, column in enumerate(coverage_df.columns):
        plt.plot(coverage_df.index, coverage_df[column], label=column, color=cmap(norm(i)))

    plt.xlabel("Gene body percentile (5'->3')", fontsize=14)
    plt.ylabel('Coverage', fontsize=14)
    plt.title('Gene Body Coverage', fontweight='bold', fontsize=16)

    # Customizing the x-axis and y-axis
    plt.xticks(ticks=range(0, 101, 20), labels=range(0, 101, 20), fontsize=12)
    plt.yticks(fontsize=12)
    plt.minorticks_on()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.grid(True, which='major', linestyle='-', linewidth=1)
    plt.legend(title='Sample', frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10,
               title_fontsize=12)
    plt.tight_layout()

    # Save as png and pdf with 300 DPI
    plt.savefig(f'{output_prefix}.geneBodyCoverage.curves.png', dpi=300, transparent=False, facecolor='white')
    plt.savefig(f'{output_prefix}.geneBodyCoverage.curves.pdf', dpi=300)
    plt.close()


def process_bam_file(args):
    '''
    Process a single BAM file to calculate coverage.
    '''
    bamfile, gene_percentiles = args
    return genebody_coverage(bamfile, gene_percentiles)


def main():
    usage = "%prog [options]" + '\n' + __doc__ + "\n"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-i", "--input", action="store", type="string", dest="input_files",
                      help='Input file(s) in BAM format. "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files. 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam file (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools.')
    parser.add_option("-r", "--refgene", action="store", type="string", dest="ref_gene_model",
                      help="Reference gene model in bed format. [required]")
    parser.add_option("-l", "--minimum_length", action="store", type="int", default=100, dest="min_mRNA_length",
                      help="Minimum mRNA length (bp). mRNA smaller than \"min_mRNA_length\" will be skipped. default=%default")
    parser.add_option("-f", "--format", action="store", type="string", dest="output_format", default='pdf',
                      help="Output file format, 'pdf', 'png' or 'jpeg'. default=%default")
    parser.add_option("-o", "--out-prefix", action="store", type="string", dest="output_prefix",
                      help="Prefix of output files(s). [required]")
    (options, args) = parser.parse_args()

    if not (options.output_prefix and options.input_files and options.ref_gene_model):
        parser.print_help()
        sys.exit(0)

    if not os.path.exists(options.ref_gene_model):
        print('\n\n' + options.ref_gene_model + " does NOT exist" + '\n', file=sys.stderr)
        sys.exit(0)
    if options.min_mRNA_length < 100:
        print('The number specified to "-l" cannot be smaller than 100.' + '\n', file=sys.stderr)
        sys.exit(0)

    OUT1 = open(options.output_prefix + ".geneBodyCoverage.txt", 'w')
    print("Percentile\t" + '\t'.join([str(i) for i in range(1, 101)]), file=OUT1)

    printlog("Read BED file (reference gene model) ...")
    gene_percentiles = genebody_percentile(refbed=options.ref_gene_model, mRNA_len_cut=options.min_mRNA_length)

    printlog("Get BAM file(s) ...")
    bamfiles = []
    if os.path.isdir(options.input_files):
        bamfiles = [os.path.join(options.input_files, f) for f in os.listdir(options.input_files) if f.endswith('.bam')]
    else:
        bamfiles = [options.input_files]

    for f in bamfiles:
        print("\t" + f, file=sys.stderr)

    file_container = []
    args = [(bamfile, gene_percentiles) for bamfile in bamfiles]

    # Use multiprocessing Pool with the number of CPUs available
    with Pool(cpu_count()) as pool:
        coverages = pool.map(process_bam_file, args)

    sample_names = []
    coverage_dict = {}

    for bamfile, cvg in zip(bamfiles, coverages):
        if len(cvg) == 0:
            print("\nCannot get coverage signal from " + basename(bamfile) + ' ! Skip', file=sys.stderr)
            continue
        sample_name = extract_sample_name(basename(bamfile))
        normalized_cvg = [(v - min(cvg.values())) / (max(cvg.values()) - min(cvg.values())) for v in cvg.values()]
        coverage_dict[sample_name] = normalized_cvg
        file_container.append(sample_name)

    # Sort sample names
    sorted_sample_names = sorted(coverage_dict.keys())

    # Write sorted data to file
    for sample_name in sorted_sample_names:
        print(sample_name + '\t' + '\t'.join(map(str, coverage_dict[sample_name])), file=OUT1)
    OUT1.close()

    dataset = [(name, coverage_dict[name], pearson_moment_coefficient(coverage_dict[name])) for name in
               sorted_sample_names]

    print("\n\n", file=sys.stderr)
    print("\tSample\tSkewness", file=sys.stderr)
    for a, b, c in dataset:
        print('\t' + a + '\t' + str(c), file=sys.stderr)

    visualize_coverage([d[1] for d in dataset], options.output_prefix, sorted_sample_names)


if __name__ == '__main__':
    main()