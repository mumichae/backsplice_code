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


# rule get_Wang2019_training_set:
#     """
#     Subtract the test data from DiLiddo2019 from Wang2019's positive training set
#     Create the negative trining dataset
#     Create the DeepCirCode Input for training
#     """
#     input:
#         wang_positive=config['processed_data'] + '/datasets/Wang2019/circRNA.bed',
#         diLiddo_positive=config['processed_data'] + '/datasets/DiLiddo2019/circRNA.bed',
#         gtf=config['processed_data'] + '/reference/hg38/hg38.ensembl.canonical.gtf',
#         fasta=config['processed_data'] + '/reference/hg38/hg38.fa',
#         dcc_dir=config['processed_data'] + '/DeepCirCode_input'
#     output:
#         wang_without_diLiddo=config['processed_data'] + '/datasets/Wang2019/Wang_without_DiLiddo/circRNA.bed',
#         negative_training=config['processed_data'] + '/negative_dataset/Wang_training_negative.bed',
#         dcc_input_tsv=config['processed_data'] + '/DeepCirCode_input/deepCirCode_train.tsv',
#         dcc_input_x=config['processed_data'] + '/DeepCirCode_input/deepCirCode_x_train.txt',
#         dcc_input_y=config['processed_data'] + '/DeepCirCode_input/deepCirCode_y_train.txt'
#     shell:
#         """
#         bedtools subtract -s -A -a {input.wang_positive} -b {input.diLiddo_positive} > {output.wang_without_diLiddo}
#         python workflow/scripts/data/get_linear_junctions.py -gtf {input.gtf} -circ {output.wang_without_diLiddo} -o {output.negative_training}
#         python workflow/scripts/data/get_DeepCirCode_input.py -g {input.fasta} -pos {output.wang_without_diLiddo} -neg {output.negative_training} -mode train -o {input.dcc_dir}
#         """


rule extract_DeepCirCode_data:
    """
    Create the DeepCirCode Input for training
    """
    input:
        script='workflow/scripts/feature_extraction/extract_DeepCirCode_input.py',
        positive=get_positive_data,
        negative=get_negative_data,
        fasta=get_fasta
    output:
        tsv=expand(feature_pattern + '/all_data.tsv',method='DeepCirCode',allow_missing=True),
        x=expand(feature_pattern + '/x_matrix.txt',method='DeepCirCode',allow_missing=True),
        y=expand(feature_pattern + '/y_matrix.txt',method='DeepCirCode',allow_missing=True)
    shell:
        """
        python {input.script} \
            -g {input.fasta} \
            -pos {input.positive} \
            -neg {input.negative} \
            -tsv {output.tsv} \
            -x {output.x} \
            -y {output.y}
        """


rule all_extract_DeepCirCode_data:
    input: expand(rules.extract_DeepCirCode_data.output,source=['DiLiddo2019', 'Wang2019'])
