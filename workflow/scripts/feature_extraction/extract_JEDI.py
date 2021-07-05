import argparse
import pandas as pd
import gtfparse
from pybedtools import BedTool
from pybedtools.helpers import set_tempdir
from Bio import SeqIO
import json
from joblib import Parallel, delayed

BED_COLNAMES = ['chrom', 'start', 'end', 'name', 'score', 'strand']


def parser():
    # uses argparse to parse the command line arguments
    # returns the parsed arguments
    parser = argparse.ArgumentParser(description='Extract sequences and junction locations for JEDI')
    parser.add_argument('-pos', type=str, required=True, dest='out_pos',
                        help='The output file for circRNA junctions (positive data)')
    parser.add_argument('-neg', type=str, required=True, dest='out_neg',
                        help='The output file for transcripts overlapping with circRNA (negative data)')
    parser.add_argument('-exons', type=str, required=True, dest='exons',
                        help='The file location of the genome sequence (must be in fasta format)')
    parser.add_argument('-transcripts', type=str, required=True, dest='transcripts',
                        help='The file location of the genome sequence (must be in fasta format)')
    parser.add_argument('-circ', type=str, required=True, dest='positive',
                        help='The file location of the positive data\'s bed file')
    parser.add_argument('-fasta', type=str, required=True, dest='fasta',
                        help='FASTA file')
    parser.add_argument('-id_key', type=str, default='transcript_id', dest='id_key',
                        help='Attribute key in GTF for transcript ID')
    parser.add_argument('-tmp', type=str, default='/tmp', dest='tmpdir',
                        help='Temporary directory path')
    parser.add_argument('-threads', type=int, default=2, dest='threads',
                        help='Number of threads for parallel processing')

    return parser.parse_args()


def parse_gtf(filename, id_key):
    gtf = gtfparse.read_gtf(filename)
    gtf['start'] = gtf['start'] - 1  # 0-based for later processing
    gtf['chrom'] = gtf['seqname']
    gtf['name'] = gtf[id_key]
    return gtf[BED_COLNAMES].drop_duplicates()


def extract_from_bed(transcripts, exons_df, fasta_file, output_path, n_jobs):
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

    def extract_features(input_tuple, exons_df):
        # unpack input
        tx_row, record = input_tuple
        _, transcript = tx_row

        # get exons regions per transcript
        exons = exons_df[
            (exons_df['chrom'] == transcript['chrom']) &
            (exons_df['strand'] == transcript['strand']) &
            (exons_df['start'] >= transcript['start']) &
            (exons_df['end'] <= transcript['end']) &
            (exons_df['name'] == transcript['name'])
            ]
        if exons.shape[0] == 0:
            exons = exons.append(transcript, ignore_index=True)

        # translate exon regions into relative positions & determine head/tail
        head = exons['start'] - transcript['start']
        tail = exons['end'] - transcript['start']

        del exons

        return dict(
            name=transcript['name'],  # optional
            chr=transcript['chrom'],  # optional
            start=transcript['start'],  # optional
            end=transcript['end'],  # optional
            strand=transcript['strand'][0],
            junctions=dict(
                head=head.sort_values().tolist(),
                tail=tail.sort_values().tolist()
            ),
            seq=record.seq.__str__().upper(),
        )

    features = Parallel(n_jobs=n_jobs, verbose=1, backend="loky", max_nbytes=10)(
        delayed(extract_features)(input_tuple, exons_df)
        for input_tuple in zip(tx_df.iterrows(), fasta_records)
    )
    print(f'{len(features)} entries computed')

    with open(output_path, 'w') as output_file:
        for feature in features:
            output_file.write(json.dumps(feature))
            output_file.write('\n')


if __name__ == "__main__":
    args = parser()

    id_key = args.id_key
    fasta_file = args.fasta
    set_tempdir(args.tmpdir)  # set tmp dir for pybedtools

    # read circRNA BED file
    circ = pd.read_table(args.positive, header=None, sep='\t')
    circ.columns = BED_COLNAMES
    circ_bed = BedTool.from_dataframe(circ)

    exons = parse_gtf(args.exons, id_key)
    transcripts = parse_gtf(args.transcripts, id_key)

    exons_bed = BedTool.from_dataframe(exons)
    transcripts_bed = BedTool.from_dataframe(transcripts)

    # subset transcripts to only overlapping with circRNA junctions
    transcripts_bed = transcripts_bed.intersect(circ_bed, s=True, u=True)

    # map linear transcript IDs to circRNA junctions
    anno_circ_df = circ_bed.intersect(transcripts_bed, s=True, wa=True, wb=True, loj=True).to_dataframe()
    anno_circ_df['name'] = anno_circ_df[['blockCount']]  # blockCount contains name from second BED
    circ_df = anno_circ_df[BED_COLNAMES].drop_duplicates([x for x in BED_COLNAMES if x != 'name'])

    print('Extract circRNA features')
    extract_from_bed(
        transcripts=circ_df,
        exons_df=exons,
        fasta_file=fasta_file,
        output_path=args.out_pos,
        n_jobs=args.threads
    )

    print('Extract linear transcript features')
    extract_from_bed(
        transcripts=transcripts_bed,
        exons_df=exons,
        fasta_file=fasta_file,
        output_path=args.out_neg,
        n_jobs=args.threads
    )
