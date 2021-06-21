include: "data.smk"

intermediate_dir = config['processed_data'] + '/intermediate'

rule train_circDeep:
    input:
        script = 'methods/circDeep/circDeep.py',
        bigwig = expand(rules.get_phastCons.output[1], assembly=config['assembly']),
        fasta = expand(rules.genomepy.output[0],assembly=config['assembly']),
        gtf = expand(rules.canonical_gtf.output, assembly=config['assembly']),
        positive = rules.data_Chaabane2020.output.positive_bed,
        negative = rules.data_Chaabane2020.output.negative_bed,
    output:
        model = directory(config['models'] + '/circDeep_retrained'),
        intermediate_data = directory(intermediate_dir + '/circDeep_retrained')
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


rule predict_circDeep:
    input:
        script = 'methods/circDeep/circDeep.py',
        bigwig = expand(rules.get_phastCons.output[1], assembly=config['assembly']),
        fasta = expand(rules.genomepy.output[0],assembly=config['assembly']),
        gtf = expand(rules.canonical_gtf.output, assembly=config['assembly']),
        seq_features = 'methods/circDeep/data/seq_features.txt',
        rcm_features = 'methods/circDeep/data/rcm_features.txt',
        conservation_features = 'methods/circDeep/data/conservation_features.txt',
        class_txt = 'methods/circDeep/data/class.txt',
        test = rules.subset_data.output.positive,
        model = 'methods/circDeep/models',  # rules.train_circDeep.output.model
    output:
        seq_features = intermediate_dir + '/circDeep/seq_features.txt',
        rcm_features= intermediate_dir + '/circDeep/rcm_features.txt',
        conservation_features= intermediate_dir + '/circDeep/conservation_features.txt',
        class_txt= intermediate_dir + '/circDeep/class.txt',
        prediction = config['evaluation'] + '/circDeep/prediction.txt'
    conda:
        "../envs/circDeep.yaml"
    threads: 12
    shell:
        """
        data_dir=$(dirname {output.seq_features})
        mkdir -p $data_dir
        cp {input.seq_features} {output.seq_features}
        cp {input.rcm_features} {output.rcm_features}
        cp {input.conservation_features} {output.conservation_features}
        cp {input.class_txt} {output.class_txt}
        python {input.script} \
            --data_dir $data_dir/ \
            --model_dir {input.model}/ \
            --predict True \
            --out_file {output.prediction} \
            --seq True --rcm True --cons True \
            --genome {input.fasta} \
            --gtf {input.gtf} \
            --bigwig {input.bigwig} \
            --testing_bed {input.test}
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
