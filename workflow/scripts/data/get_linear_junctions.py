""" Get Linear Junctions

This script creates the coordinates of the circRNA negative dataset using:
- the high confidence circRNA coordinates from Andre's dataset.
- the gtf annotation of the human genome.

This script requires to be run within our conda environment.

Parameters
----------
    -gtf : str
        The file location of the gtf annotation
    -circ : str
        The file location of the high confidence circRNA positions
    -o : str
        The output file

Returns
-------
    negative_dataset.bed
        positions of all introns in the transcripts overlapping with the circRNA positions
"""

import argparse


def parser():
    # uses argparse to parse the command line arguments
    # returns the parsed arguments
    parser = argparse.ArgumentParser(description='Create Positive Dataset')
    parser.add_argument('-gtf', type=str, required=True, dest='gtf_file',
                        help='The file location of the gtf annotation')
    parser.add_argument('-circ', type=str, required=True, dest='circ_file',
                        help='The file location of the high confidence circRNA positions')
    parser.add_argument('-o', type=str, required=True, dest='output_path',
                        help='The output BED file')

    return parser.parse_args()


def process_circ(file):
    # reads the high confidence circRNA locations and puts them into a dictionary
    # format: {'chr:start-end:strand':
    #           chr, circ_start, circ_end, strand, linear1_start, linear1_end, linear2_start, linear2_end}
    circ = {}
    for line in file.readlines()[1:]:
        tabs = line.split("\t")
        circ_entry = tabs[0].split(":")
        val = [circ_entry[0],  # chr
               int(circ_entry[1].split("-")[0]),  # start, 0-based, inclusive
               int(circ_entry[1].split("-")[1]),  # end, 0-based, inclusive
               circ_entry[2],  # strand
               int(tabs[1].split(":")[1].split("-")[0]),  # linear1_start, 1-based, inclusive
               int(tabs[1].split(":")[1].split("-")[1]),  # linear1_end, 1-based, inclusive
               int(tabs[2].split(":")[1].split("-")[0]),  # linear2_start, 1-based, inclusive
               int(tabs[2].split(":")[1].split("-")[1])]  # linear2_end, 1-based, inclusive
        circ[tabs[0]] = val
    return circ


def get_introns(circ, gtf_file):
    # iterates through gtf file, spots transcripts that overlap with our circRNAs, writes down these introns,
    # as long as they're not already in the set of introns found.
    # introns: [chr, start, end, strand, transcript_id, intron_number]; 1-based, start and end inclusive
    introns = {}
    circ_flag = False
    end_last_exon = -1
    intron_number = 0

    for line in gtf_file.readlines():
        tabs = line.split("\t")
        chr = tabs[0]
        start = int(tabs[3])  # start and end are 1-based in GTF, and both are inclusive
        end = int(tabs[4])
        strand = tabs[6]

        if tabs[2] == "transcript":  # new transcript, check if it overlaps with a circRNA
            circ_flag = False
            end_last_exon = -1
            for circRNA in circ.values():
                if circRNA[0] == chr and circRNA[1] >= start - 1 and circRNA[2] <= end - 1 and circRNA[3] == strand:
                    # means, we found a circRNA in this transcript
                    circ_flag = True
                    intron_number = 0
                    # print(
                    #    "in transcript " + chr + " " + str(start) + "-" + str(end) + ":" + strand +
                    #    " there is the circRNA " + str(circRNA))
                    break

        elif tabs[2] == "exon" and circ_flag:
            if end_last_exon != -1:
                # means, this is not the first exon of the transcript, so write down the intron
                transcript_id = tabs[8].split("\"")[3]
                intron_number += 1
                introns[chr + ":" + str(end_last_exon + 1) + "-" + str(start - 1) + ":" + strand] =\
                    [chr, end_last_exon + 1, start - 1, strand, transcript_id, intron_number]
                # due to the dictionary structure duplicate exons are overwritten and don't appear twice in the output
            end_last_exon = end

    return introns


if __name__ == "__main__":
    # 1. parse parameters
    args = parser()

    # 2. read files and process input data
    circRNA_file = open(args.circ_file, "r")
    gtf_file = open(args.gtf_file, "r")
    output = open(args.output_path, "w")

    circRNAs = process_circ(circRNA_file)

    # 3. get the introns from overlapping transcripts
    introns = get_introns(circRNAs, gtf_file)

    # 4. write the introns to the output file
    output.write("# chr\tstart\tend\tstrand\ttranscript_id\tintron_number\n")
    for intron in introns.values():
        # print(intron)
        output.write(str(intron[0]) + "\t" + str(intron[1]) + "\t" + str(intron[2]) + "\t" + str(intron[3]) + "\t" +
                     str(intron[4]) + "\t" + str(intron[5]) + "\n")
    output.flush()
    output.close()
