rule train:
    """
    Train and optimise model on training data
    """
    input:
        train=rules.split_train_test.output.train
    output:
        model=config['models'] + '/{method}/model_{language}.{suffix}'
    script: '../scripts/models/train.{language}'


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
    script: '../scripts/models/predict.{language}'


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
