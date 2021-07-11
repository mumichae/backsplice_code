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
                        help='The output file for gene sequences (negative data)')
    parser.add_argument('-exons', type=str, required=True, dest='exons',
                        help='GTF of exons for overlapping gene/circRNA regions')
    parser.add_argument('-genes', type=str, required=True, dest='genes',
                        help='GTF of genes considered for negative dataset')
    parser.add_argument('-circ', type=str, required=True, dest='positive',
                        help='The file location of the positive data\'s bed file')
    parser.add_argument('-fasta', type=str, required=True, dest='fasta',
                        help='FASTA file')
    parser.add_argument('-id_key', type=str, default='gene_id', dest='id_key',
                        help='Attribute key in GTF for gene ID')
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


def extract_from_bed(genes, exons_df, fasta_file, output_path, n_jobs):
    """
    Writes one dictionary per line with
        sequence: sequence of gene (including introns)
        strand: + or -
        junctions:
            head: list of splice acceptor sites relative to gene
            tail: list of splice donor sites relative to gene

    :param genes: pybedtools.BedTool or pandas dataframe (0-based) of genes/circRNA BSJ
    :param exons_df: pandas dataframe (0-based) of all exons in annotation
    :param output_path: path of file to write parsed information to
    """
    if isinstance(genes, BedTool):
        gene_df = genes.to_dataframe()
        gene_bed = genes
    elif isinstance(genes, pd.DataFrame):
        gene_df = genes
        gene_bed = BedTool.from_dataframe(genes)
    else:
        raise TypeError('invalid type for genes')

    # extract fasta sequences of genes
    gene_bed = gene_bed.sequence(fi=fasta_file)
    fasta_records = SeqIO.parse(gene_bed.seqfn, "fasta")

    def extract_features(input_tuple, exons_df):
        # unpack input
        gene, record = input_tuple
        _, gene = gene

        # get exons regions per gene
        exons = exons_df[
            (exons_df['chrom'] == gene['chrom']) &
            (exons_df['strand'] == gene['strand']) &
            (exons_df['start'] >= gene['start']) &
            (exons_df['end'] <= gene['end']) &
            (exons_df['name'] == gene['name'])
            ]
        if exons.shape[0] == 0:
            exons = exons.append(gene, ignore_index=True)

        # translate exon regions into relative positions & determine head/tail
        head = exons['start'] - gene['start']
        tail = exons['end'] - gene['start']

        del exons

        return dict(
            name=gene['name'],  # optional
            chr=gene['chrom'],  # optional
            start=gene['start'],  # optional
            end=gene['end'],  # optional
            strand=gene['strand'][0],
            junctions=dict(
                head=head.sort_values().tolist(),
                tail=tail.sort_values().tolist()
            ),
            seq=record.seq.__str__().upper(),
        )

    features = Parallel(n_jobs=n_jobs, verbose=1, backend="loky", max_nbytes=10, batch_size=300)(
        delayed(extract_features)(input_tuple, exons_df)
        for input_tuple in zip(gene_df.iterrows(), fasta_records)
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
    genes = parse_gtf(args.genes, id_key)

    exons_bed = BedTool.from_dataframe(exons)
    genes_bed = BedTool.from_dataframe(genes)

    # map linear gene IDs to circRNA junctions
    anno_circ_df = circ_bed.intersect(genes_bed, s=True, wa=True, wb=True, loj=True).to_dataframe()
    anno_circ_df['name'] = anno_circ_df[['blockCount']]  # blockCount contains name from second BED
    circ_df = anno_circ_df[BED_COLNAMES].drop_duplicates([x for x in BED_COLNAMES if x != 'name'])

    print('Extract circRNA features')
    extract_from_bed(
        genes=circ_df,
        exons_df=exons,
        fasta_file=fasta_file,
        output_path=args.out_pos,
        n_jobs=args.threads
    )

    print('Extract linear gene features')
    extract_from_bed(
        genes=genes_bed,
        exons_df=exons,
        fasta_file=fasta_file,
        output_path=args.out_neg,
        n_jobs=args.threads
    )
