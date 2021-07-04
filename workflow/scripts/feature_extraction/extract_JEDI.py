import argparse
import pandas as pd
import gtfparse
from pybedtools import BedTool
from pybedtools.helpers import set_tempdir
from Bio import SeqIO
import json

BED_COLNAMES = ['chrom', 'start', 'end', 'name', 'score', 'strand']


def parser():
    # uses argparse to parse the command line arguments
    # returns the parsed arguments
    parser = argparse.ArgumentParser(description='Extract sequences and junction locations for JEDI')
    parser.add_argument('-pos', type=str, required=True, dest='out_pos',
                        help='The output file for circRNA junctions (positive data)')
    parser.add_argument('-neg', type=str, required=True, dest='out_neg',
                        help='The output file for transcripts overlapping with circRNA (negative data)')
    parser.add_argument('-gtf', type=str, required=True, dest='gtf',
                        help='The file location of the genome sequence (must be in fasta format)')
    parser.add_argument('-circ', type=str, required=True, dest='positive',
                        help='The file location of the positive data\'s bed file')
    parser.add_argument('-fasta', type=str, required=True, dest='fasta',
                        help='FASTA file')
    parser.add_argument('-tx', type=str, default='transcript_id', dest='tx_col',
                        help='Attribute key in GTF for transcript ID')
    parser.add_argument('-tmp', type=str, default='/tmp', dest='tmpdir',
                        help='Temporary directory path')

    return parser.parse_args()


def extract_from_bed(transcripts, exons_df, fasta_file, output_path):
    """
    Writes one dictionary per line with
        sequence: sequence of transcript (including introns)
        strand: + or -
        junctions:
            head: list of splice acceptor sites relative to transcript
            tail: list of splice donor sites relative to transcript

    :param transcripts: pybedtools.BedTool or pandas dataframe (0-based) of transcripts/circRNA BSJ
    :param exons_df: pandas dataframe (0-based) of all exons in annotation
    :param output_path: path of file to write parsed information to
    """
    if isinstance(transcripts, BedTool):
        tx_df = transcripts.to_dataframe()
        tx_bed = transcripts
    elif isinstance(transcripts, pd.DataFrame):
        tx_df = transcripts
        tx_bed = BedTool.from_dataframe(transcripts)
    else:
        raise TypeError('invalid type for transcripts')

    # extract fasta sequences of transcripts
    tx_bed = tx_bed.sequence(fi=fasta_file)
    fasta_records = SeqIO.parse(tx_bed.seqfn, "fasta")

    with open(output_path, 'w') as output_file:
        for ((_, transcript), record) in zip(tx_df.iterrows(), fasta_records):
            # 2. get sequences of transcripts + strand
            sequence = record.seq.__str__().upper()
            strand = transcript['strand'][0]

            # 3. get exons regions per transcript
            exons = exons_df[
                (exons_df['chrom'] >= transcript['chrom']) &
                (exons_df['start'] >= transcript['start']) &
                (exons_df['end'] <= transcript['end']) &
                (exons_df['name'] == transcript['name'])
            ]
            if exons.shape[0] == 0:  # skip if no exons overlap
                continue

            # 4. translate exon regions into relative positions & determine head/tail
            head = []
            tail = []
            for _, exon in exons.iterrows():
                head.append(exon['start'] - transcript['start'])
                tail.append(exon['end'] - transcript['start'])
            head.sort()
            tail.sort()

            # 5. write JSON entry to file
            entry = dict(
                name=transcript['name'],  # optional
                chr=transcript['chrom'],  # optional
                start=transcript['start'],  # optional
                end=transcript['end'],  # optional
                strand=strand,
                junctions=dict(
                    head=head,
                    tail=tail
                ),
                seq=sequence,
            )
            output_file.write(json.dumps(entry))
            output_file.write('\n')


if __name__ == "__main__":
    args = parser()

    tx_col = args.tx_col
    fasta_file = args.fasta
    set_tempdir(args.tmpdir)  # set tmp dir for pybedtools

    # read circRNA BED file
    circ = pd.read_table(args.positive, header=None, sep='\t')
    circ.columns = BED_COLNAMES
    circ_bed = BedTool.from_dataframe(circ)

    # parse GTF
    gtf = gtfparse.read_gtf(args.gtf)
    gtf['start'] = gtf['start'] - 1  # 0-based for later processing

    # separate exons and transcripts
    gtf['chrom'] = gtf['seqname']
    gtf['name'] = gtf[tx_col]
    exons = gtf[gtf.feature == 'exon'][BED_COLNAMES].drop_duplicates()
    exons_bed = BedTool.from_dataframe(exons)
    transcripts = gtf[gtf.feature == 'transcript'][BED_COLNAMES].drop_duplicates()
    transcripts_bed = BedTool.from_dataframe(transcripts)

    # subset transcripts to only overlapping with circRNA junctions
    transcripts_bed = transcripts_bed.intersect(circ_bed, s=True, u=True)

    # map linear transcript IDs to circRNA junctions
    anno_circ_df = circ_bed.intersect(transcripts_bed, s=True, wa=True, wb=True, loj=True).to_dataframe()
    anno_circ_df['name'] = anno_circ_df[['blockCount']]  # blockCount contains name from second BED
    circ_bed = anno_circ_df[BED_COLNAMES].drop_duplicates([x for x in BED_COLNAMES if x != 'name'])

    print('Extract circRNA features')
    extract_from_bed(
        transcripts=circ_bed,
        exons_df=exons,
        fasta_file=fasta_file,
        output_path=args.out_pos
    )

    print('Extract linear transcript features')
    extract_from_bed(
        transcripts=transcripts_bed,
        exons_df=exons,
        fasta_file=fasta_file,
        output_path=args.out_neg
    )
