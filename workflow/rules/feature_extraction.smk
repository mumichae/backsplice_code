include: "data.smk"

feature_pattern = config['processed_data'] + '/features/{method}/{source}'


rule extract_JEDI_raw:
    input:
        positive_bed=get_positive_data,
        gtf=rules.canonical_gtf.output[0],
        fasta=get_fasta
    output:
        positive=expand(feature_pattern + '_human_gene.pos', method='JEDI', allow_missing=True),
        negative=expand(feature_pattern + '_human_gene.neg', method='JEDI', allow_missing=True)
    script: '../scripts/feature_extraction/extract_JEDI.py'

