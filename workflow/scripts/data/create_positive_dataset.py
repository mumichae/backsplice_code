"""Create Positive Dataset

This script creates the circRNA positive dataset using:
- the human genome sequence in hg 38,
- the bed file from circBase (liftover to hg38 already done with
  the UCSC tool (http://genome.ucsc.edu/cgi-bin/hgLiftOver)), and
- the bed file from circDeep, which originates from circRNADb.

This script requires to be run within our conda environment.

Parameters
----------
    -g : str
        The file location of the genome sequence (must be in fasta format)
    -cb : str
        The file location of the circBase bed file
    -cd : str
        The file location of the circDeep bed file
    -o : str
        The output folder

Returns
-------
    positive_dataset.fa
        Sequences of all different circRNAs found in the input files
    positive_dataset.bed
        bed file of the different circRNAs found in the input files
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def parser():
    # uses argparse to parse the command line arguments
    # returns the parsed arguments
    parser = argparse.ArgumentParser(description='Create Positive Dataset')
    parser.add_argument('-g', type=str, required=True, dest='genome',
                        help='The file location of the genome sequence (must be in fasta format)')
    parser.add_argument('-cb', type=str, required=True, dest='circBase_bed',
                        help='The file location of the circBase bed file')
    parser.add_argument('-cd', type=str, required=True, dest='circDeep_bed',
                        help='The file location of the circDeep bed file')
    parser.add_argument('-o', type=str, required=True, dest='output_path',
                        help='The output folder')

    return parser.parse_args()


def process_bed(circBase_bed, circDeep_bed):
    # reads bed files, processes the input data and removes duplicates based on their coordinates
    # note: circBase_bed's format: chr  start   end name    [1000]    strand  [...] (12 columns)
    #       circDeep_bed's format: chr  start   end strand  name    (5 columns)
    # and the circRNA's names are not consistent.
    # rule: use circ_base's annotation if available. If not, use circDeep's name, but with suffix '*'.
    circRNAs = {}

    for line in circBase_bed.readlines():
        line = line.strip("\n").split('\t')
        chr = line[0]
        start = line[1]
        end = line[2]
        strand = line[5]
        name = line[3]
        # in the bed files, the chromosome is sometimes eg "chr17_GL000258v2_alt" instead of "chr17"
        # this would lead to key errors later on
        # -> take substring if required
        if "alt" in chr:
            chr = chr.partition("_")[0]     # use string before the first occurrence of "_"
        circRNAs[chr + '_' + start + '-' + end + '_' + strand] = (chr, start, end, strand, name)

    for line in circDeep_bed.readlines():
        line = line.strip("\n").split("\t")
        chr = line[0]
        start = line[1]
        end = line[2]
        strand = line[3]
        name = line[4]
        # in the bed files, the chromosome is sometimes eg "chr17_GL000258v2_alt" instead of "chr17"
        # this would lead to key errors later on
        # -> take substring if required
        if "alt" in chr:
            chr = chr.partition("_")[0]  # use string before the first occurrence of "_"
        key = chr + '_' + start + '-' + end + '_' + strand
        if key not in circRNAs:
            circRNAs[key] = (chr, start, end, strand, name+'*')

    return circRNAs


def get_circRNA_sequences (chromosomes, circRNAs, output_path):
    # reads the circRNA sequences from the genome dictionary 'chromosomes'
    # and writes them as fasta file to the output destination
    # note: bed files use 0-based, but exclusive start position
    output = open(output_path + "/positive_dataset.fa", "w")
    line_length = 80
    for key in circRNAs:
        circ = circRNAs[key]    # being (chr, start, end, strand, name)
        sequence = chromosomes[circ[0]][int(circ[1])+1:int(circ[2])]
        if circ[3] == '-':
            # reverse complement of sequence
            sequence = sequence.reverse_complement()


        # write to output file
        sequence = sequence.seq
        output.write(">" + circ[4] + "\n")
        while len(sequence) > line_length:
            output.write(str(sequence[:line_length]) + "\n")
            sequence = sequence[line_length:]
        output.write(str(sequence) + "\n")
    output.flush()
    output.close()



if __name__ == "__main__":

    # 1. parse parameters
    args = parser()

    # 2. read bed files and process input data
    circBase_bed = open(args.circBase_bed, "r")
    circDeep_bed = open(args.circDeep_bed, "r")
    circRNAs = process_bed(circBase_bed, circDeep_bed)
    print("# different circRNAs: " + str(len(circRNAs)))

    # 3. print common bed file
    output_bed = open(args.output_path + '/positive_dataset.bed', "w")
    for key in circRNAs:
        value = circRNAs[key]
        output_bed.write(value[0] + "\t" + value[1] + "\t" + value[2] + "\t" + value[3] + "\t" + value[4] + "\n")
    output_bed.flush()
    output_bed.close()

    # 4. get sequences and write them to the output file
    chromosomes = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    get_circRNA_sequences(chromosomes, circRNAs, args.output_path)






