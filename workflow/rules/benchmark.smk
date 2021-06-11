rule circDeep:
    input:
        script = 'methods/circDeep/circDeep.py',
        bigwig = expand(rules.get_phastCons.output[1], assembly=config['assembly']),
        fasta = expand(rules.genomepy.output[0],assembly=config['assembly']),
        gtf = expand(rules.canonical_gtf.output, assembly=config['assembly']),
        positive = rules.Chaabane2020.output.positive_bed,
        negative = rules.Chaabane2020.output.negative_bed,
    output:
        model = directory(config['models'] + '/circDeep'),
        intermediate_data = directory(config['processed_data'] + '/intermediate/circDeep')
    params:
        data_dir = 'methods/circDeep/data/'
    conda:
        "../envs/circDeep.yaml"
    threads: 12
    shell:
        """
        python {input.script} \
            --data_dir {output.intermediate_data} \
            --train True \
            --model_dir '{output.model}' \
            --seq True --rcm True --cons True \
            --genome '{input.fasta}' \
            --gtf '{input.gtf}' \
            --bigwig '{input.bigwig}' \
            --positive_bed '{input.positive}' \
            --negative_bed '{input.negative}'
        """


rule train:
    """
    Train and optimise model on training data
    """
    input:
        train=rules.split_train_test.output.train
    output:
        model=config['models'] + '/{method}/model_{language}.{suffix}'
    script: '../scripts/models/train.{wildcards.language}'


rule predict:
    # TODO: generalise for different models, using wildcards
    """
    Prediction and evaluation on test data
    Output metrics and processed_data
    """
    input:
        model=rules.train.output.model,
        test=rules.split_train_test.output.test
    output:
        prediction=config['models'] + '/{method}/prediction_{suffix}_{language}.tsv'
    script: '../scripts/models/predict.{wildcards.language}'


rule evaluation:
    """
    Compute evaluation metrics on all model predictions
    Collect all predictions in this rule
    """
    input:
        predictions=expand(
            rules.predict.output.prediction, zip, **get_wildcards(params)
        )
    output:
        metrics=config['evaluation'] + '/metrics.tsv'
    script: '../scripts/evaluation/metrics.py'


rule all_benchmark:
    """
    Collect all output of the benchmark
    """
    input:
        metrics=rules.evaluation.output
        # processed_data=...
