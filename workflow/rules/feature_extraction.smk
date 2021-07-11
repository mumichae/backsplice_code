include: "data.smk"


def get_train_test(wildcards, pattern, train_test=None, source=None):
    """
    Translate train/test to dataset (source) identifier
    wildcards must contain 'source'
    """
    if train_test is None:
        try:
            train_test = wildcards.train_test
        except:
            raise LookupError(f"Invalid train_test as wildcard or parameter: {train_test}")
    if source is None:
        try:
            source = wildcards.source
        except:
            raise LookupError(f"Invalid source as wildcard or parameter: {source}")

    if train_test == 'test':
        # NoChr has it's own test set
        # Default test set is DiLiddo2019
        source = 'DiLiddo2019' if source != 'NoChr' else 'NoChr_test'
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
        genes=lambda w: get_negative_data(w,method='JEDI'),
        fasta=get_fasta
    output:
        positive=expand(feature_pattern + '/human_gene.pos',method='JEDI',allow_missing=True),
        negative=expand(feature_pattern + '/human_gene.neg',method='JEDI',allow_missing=True),
        config=expand(feature_pattern + '/config.yaml',method='JEDI',allow_missing=True),
    params:
        path_data=expand(feature_pattern,method='JEDI',allow_missing=True),
        path_pred=expand(prediction_pattern,method='JEDI',allow_missing=True),
        id_key=config['gene_annotations'][assembly]['gene_column']
    threads: 40
    resources:
        mem_mb=100000
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
            '-genes {input.genes} '
            '-circ {input.circ} '
            '-fasta {input.fasta} '
            '-pos {output.positive} '
            '-neg {output.negative} '
            '-id_key {params.id_key} '
            '-threads {threads} '
        )


rule all_extract_data_JEDI:
    input: expand(rules.extract_data_JEDI.output,source=all_sources)


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
    resources:
        mem_mb=100000
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
        expand(rules.collect_features_JEDI.input,source=all_sources)


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
    input: expand(rules.extract_DeepCirCode_data.output,source=all_sources)


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
