include: "data.smk"

feature_pattern = config['processed_data'] + '/features/{method}/{source}'
pred_pattern = config['evaluation'] + '/{method}/{source}'


def get_train_test(wildcards, pattern, train_test=None):
    """
    Translate train/test to dataset (source) identifier
    wildcards must contain 'source'
    """
    if train_test is None:
        try:
            train_test = wildcards.train_test
        except:
            raise LookupError(f"Invalid train_test as wildcard or parameter: {train_test}")
    if train_test == 'test':
        source = 'DiLiddo2019'
    elif train_test == 'train':
        source = wildcards.source
    else:
        raise ValueError(f'Invalid wildcard for test_train: {wildcards.train_test}')
    return expand(pattern,source=source,allow_missing=True)


rule extract_data_JEDI:
    input:
        script='workflow/scripts/feature_extraction/extract_JEDI.py',
        circ=get_positive_data,
        exons=rules.canonical_gtf.output.exons,
        transcripts=rules.canonical_gtf.output.transcripts,
        fasta=get_fasta
    output:
        positive=expand(feature_pattern + '/human_isoform.pos',method='JEDI',allow_missing=True),
        negative=expand(feature_pattern + '/human_isoform.neg',method='JEDI',allow_missing=True),
        config=expand(feature_pattern + '/config.yaml',method='JEDI',allow_missing=True),
    params:
        path_data=expand(feature_pattern,method='JEDI',allow_missing=True),
        path_pred=expand(pred_pattern,method='JEDI',allow_missing=True),
        id_key=config['gene_annotations'][assembly]['transcript_column']
    threads: 60
    run:
        shell('mkdir -p {params.path_data}')
        shell('mkdir -p {params.path_pred}')
        path_data = shell('realpath {params.path_data}',read=True)
        path_pred = shell('realpath {params.path_pred}',read=True)

        with open(output.config[0].__str__(),'w') as fn:
            fn.write(f'path_data: {path_data}')
            fn.write(f'path_pred: {path_pred}')

        shell(
            'python {input.script} '
            '-exons {input.exons} '
            '-transcripts {input.transcripts} '
            '-circ {input.circ} '
            '-fasta {input.fasta} '
            '-pos {output.positive} '
            '-neg {output.negative} '
            '-id_key {params.id_key} '
            '-threads {threads} '
        )


rule all_extract_data_JEDI:
    input: expand(rules.extract_data_JEDI.output,source=['DiLiddo2019', 'Wang2019'])


rule extract_features_JEDI:
    """
    TODO: create multiple CV folds
    """
    input:
        script='methods/JEDI/src/generate_input.py',
        positive=lambda w: get_train_test(w,rules.extract_data_JEDI.output.positive),
        negative=lambda w: get_train_test(w,rules.extract_data_JEDI.output.negative),
    output:
        features=expand(feature_pattern + '/data.0.K{K}.L{L}.{train_test}',method='JEDI',allow_missing=True),
    shell:
        """
        python {input.script} {input.positive} {input.negative} {wildcards.K} {wildcards.L} {output}
        # head -100 {input.positive} > /tmp/pos.trunc
        # head -100 {input.negative} > /tmp/neg.trunc
        # python {input.script} /tmp/pos.trunc /tmp/neg.trunc {wildcards.K} {wildcards.L} {output}
        """


rule collect_features_JEDI:
    input:
        expand(
            rules.extract_features_JEDI.output,
            method='JEDI',
            K=config['methods']['JEDI']['kmer_len'],
            L=config['methods']['JEDI']['flank_len'],
            train_test=['train', 'test'],
            allow_missing=True
        )


rule all_extract_features_JEDI:
    input:
        expand(rules.collect_features_JEDI.input,source=['Wang2019'])


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
        features=expand(feature_pattern + '/x_matrix.txt',method='DeepCirCode',allow_missing=True),
        labels=expand(feature_pattern + '/y_matrix.txt',method='DeepCirCode',allow_missing=True)
    shell:
        """
        python {input.script} \
            -g {input.fasta} \
            -pos {input.positive} \
            -neg {input.negative} \
            -tsv {output.tsv} \
            -x {output.features} \
            -y {output.labels}
        """


rule all_extract_DeepCirCode_data:
    input: expand(rules.extract_DeepCirCode_data.output,source=['DiLiddo2019', 'Wang2019'])


rule SVM_RF_features:
    """
    run the R script extracting the features for SVM and RF from the DeepCirCode input data
    """
    input:
        test_data=lambda w: get_train_test(w,rules.extract_DeepCirCode_data.output.tsv,train_test='test')[0],
        train_data=lambda w: get_train_test(w,rules.extract_DeepCirCode_data.output.tsv,train_test='train')[0],
    output:
        test_features=expand(feature_pattern + 'test.rds',method='SVM_RF',allow_missing=True),
        train_features=expand(feature_pattern + 'train.rds',method='SVM_RF',allow_missing=True),
    script:
        '../scripts/data/wang_2019.R'
