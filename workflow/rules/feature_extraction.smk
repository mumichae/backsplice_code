include: "data.smk"

feature_pattern = config['processed_data'] + '/features/{method}/{source}'


rule extract_data_JEDI:
    input:
        script='workflow/scripts/feature_extraction/extract_JEDI.py',
        circ=get_positive_data,
        gtf=rules.canonical_gtf.output[0],
        fasta=get_fasta
    output:
        positive=expand(feature_pattern + '/human_gene.pos',method='JEDI',allow_missing=True),
        negative=expand(feature_pattern + '/human_gene.neg',method='JEDI',allow_missing=True)
    params:
        tx=config['gene_annotations'][assembly]['transcript_column']
    shell:
        """
        python {input.script} \
            -gtf {input.gtf} \
            -circ {input.circ} \
            -fasta {input.fasta} \
            -pos {output.positive} \
            -neg {output.negative} \
            -tx {tx}
        """


rule all_extract_JEDI_data:
    input: expand(rules.extract_data_JEDI.output,source=['DiLiddo2019', 'Wang2019'])
