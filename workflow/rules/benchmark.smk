rule train_DeepCirCode:
    # TODO: generalise for different methods, using wildcards
    """
    Train and optimise model on training data
    """
    input:
        train=rules.split_train_test.output.train
    output:
        model=config['models'] + '/DeepCirCode/model.RDS'
    script: '../scripts/methods/DeepCirCode_train.R'


rule predict_DeepCirCode:
    # TODO: generalise for different methods, using wildcards
    """
    Prediction and evaluation on test data
    Output metrics and processed_data
    """
    input:
        model=rules.train_DeepCirCode.output.model,
        test=rules.split_train_test.output.test
    output:
        prediction=config['models'] + '/DeepCirCode/prediction.tsv'
    script: '../scripts/methods/DeepCirCode_predict.R'


rule evaluation:
    """
    Compute evaluation metrics on all model predictions
    Collect all predictions in this rule
    """
    input:
        predictions=expand(rules.predict_DeepCirCode.output.prediction)
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
