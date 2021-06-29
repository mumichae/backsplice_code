"""Get DeepCirCode Input

This script extracts the input files for DeepCirCode, using:
- the human genome sequence in hg 38,
- the bed file containing the positive data, i.e. circRNA junctions, and
- the bed file containing the negative data, i.e. linear junctions.

This script requires to be run within our conda environment.

Parameters
----------
    -g : str
        The file location of the genome sequence (must be in fasta format)
    -pos : str
        The file location of the positive data's bed file
    -neg : str
        The file location of the negative data's bed file
    -o : str
        The output folder

Returns
-------
    deepCirCode_test.tsv
        table of the input data, containing the label, sequence, encoded sequence and junction name
    deepCirCode_x_test.txt
        only the encoded sequences
    deepCirCode_y_test.txt
        only the label: 01 for positive, 10 for negative

"""

import argparse
from Bio import SeqIO


def parser():
    # uses argparse to parse the command line arguments
    # returns the parsed arguments
    parser = argparse.ArgumentParser(description='Get DeepCirCode Input')
    parser.add_argument('-g', type=str, required=True, dest='genome',
                        help='The file location of the genome sequence (must be in fasta format)')
    parser.add_argument('-pos', type=str, required=True, dest='positive',
                        help=' The file location of the positive data´s bed file')
    parser.add_argument('-neg', type=str, required=True, dest='negative',
                        help=' The file location of the positive data´s bed file')
    parser.add_argument('-o', type=str, required=True, dest='output_path',
                        help='The output folder')

    return parser.parse_args()


def process_bed(positive, negative):
    # reads bed files, processes the input data and puts negative and positive data together
    dataset = []

    for line in positive.readlines():
        # 0-based, both start and end are inclusive
        # format: chr \t start \t end \t name \t '.' \t strand \n
        if line[0] != "#":
            line = line.strip("\n").split('\t')
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            strand = line[5]
            name = line[3]
            # in the bed files, the chromosome is sometimes eg "chr17_GL000258v2_alt" instead of "chr17"
            # this would lead to key errors later on
            # -> take substring if required
            # if "alt" in chr:
            #    chr = chr.partition("_")[0]     # use string before the first occurrence of "_"

            # '1' is the positive label
            dataset.append(['1', chr, int(start), int(end), strand, name])

    for line in negative.readlines():
        # 1-based, both start and end are inclusive
        # format: chr \t start \t end \t name \t '.' \t strand \n
        if line[0] != "#":
            line = line.strip("\n").split("\t")
            chr = line[0]
            start = int(line[1]) - 1  # -1 to make it 0-based
            end = int(line[2]) - 1    # -1 to make it 0-based
            strand = line[5]
            name = line[3]
            # in the bed files, the chromosome is sometimes eg "chr17_GL000258v2_alt" instead of "chr17"
            # this would lead to key errors later on
            # -> take substring if required
            # if "alt" in chr:
            #    chr = chr.partition("_")[0]  # use string before the first occurrence of "_"

            # '0' is the negative dataset label
            dataset.append(['0', chr, int(start), int(end), strand, name])

    return dataset


def get_sequences(chromosomes, dataset):
    # reads the sequences flanking the junctions from the genome dictionary 'chromosomes'
    # input is 0-based, and inclusive on start and end

    output = []
    for junction in dataset:
        up_i = chromosomes[junction[1]][int(junction[2])-50: int(junction[2])]
        down_e = chromosomes[junction[1]][int(junction[2]): int(junction[2])+50]
        up_e = chromosomes[junction[1]][int(junction[3])-49: int(junction[3])+1]
        down_i = chromosomes[junction[1]][int(junction[3])+1: int(junction[3]+51)]
        sequence = up_i + down_e + up_e + down_i

        if junction[4] == '-':
            # reverse complement of sequence
            sequence = sequence.reverse_complement()

        sequence = str(sequence.seq).upper()
        output.append([junction[0], sequence, junction[5]])
    return output


def encode_and_write(dataset, output, output_x, output_y):
    # creates one-hot-encoding for the extracted sequences and writes everything to the output file
    encoder = {"A": "1000",
               "T": "0100",
               "G": "0010",
               "C": "0001"}
    output.write("label\tRNA_input_seq\tencoded_seq\tname\n")
    for data_point in dataset:
        enc = ""
        for residue in data_point[1]:
            enc += encoder[residue]
        output.write(data_point[0] + "\t" + data_point[1] + "\t" + enc + "\t" + data_point[2] + "\n")
        output_x.write(enc + "\n")
        output_y.write(str(abs(int(data_point[0])-1)) + data_point[0] + "\n")


if __name__ == "__main__":

    # 1. parse parameters
    args = parser()

    # 2. read bed files and process input data
    positive_bed = open(args.positive, "r")
    negative_bed = open(args.negative, "r")
    dataset = process_bed(positive_bed, negative_bed)

    # 3. get sequences
    chromosomes = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    dataset = get_sequences(chromosomes, dataset)

    # 4. get one-hot-encoding and write output
    output = open(args.output_path + "/deepCirCode_test.tsv", "w")
    output_x = open(args.output_path + "/deepCirCode_x_test.txt", "w")
    output_y = open(args.output_path + "/deepCirCode_y_test.txt", "w")
    encode_and_write(dataset, output, output_x, output_y)
    output.flush()
    output_x.flush()
    output_y.flush()
    output.close()
    output_x.close()
    output_y.close()
