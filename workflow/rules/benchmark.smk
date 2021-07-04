include: "data.smk"
include: "feature_extraction.smk"

tmpdir = config['processed_data'] + '/tmp'
prediction_pattern = config['evaluation'] + '/{method}/{source}_prediction.txt'


rule train_circDeep:
    input:
        script = 'methods/circDeep/circDeep.py',
        bigwig = expand(rules.get_phastCons.output[1], assembly=config['assembly']),
        fasta = expand(rules.genomepy.output[0],assembly=config['assembly']),
        gtf = expand(rules.canonical_gtf.output, assembly=config['assembly']),
        positive = rules.data_Chaabane2020.output.positive_bed,
        negative = rules.data_Chaabane2020.output.negative_bed,
    output:
        model = config['models'] + '/circDeep/{source}/bestmodel_ACNN_BLSTM_3 1 408000.hdf5',
        intermediate_data = directory(tmpdir + '/train/circDeep/{source}')
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
        test = rules.data_dev.output.positive,
        model = rules.train_circDeep.output.model # 'methods/circDeep/models'
    output:
        seq_features = tmpdir + '/predict/circDeep/{source}/seq_features.txt',
        rcm_features=tmpdir + '/predict/circDeep/{source}/predict/rcm_features.txt',
        conservation_features=tmpdir + '/predict/circDeep/{source}/predict/conservation_features.txt',
        class_txt=tmpdir + '/predict/circDeep/{source}/predict/class.txt',
        prediction = expand(prediction_pattern, method="circDeep", allow_missing=True)
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

rule train_RF:
    """
    trains the Random Forest on given data
    """
    input: 
        train_features=config['processed_data']+'/features/SVM_RF/train.rds',
        train_labels=config['processed_data']+'/features/DeepCirCode/Wang2019/y_matrix.txt'
    output:
        RF_model=config['processed_data']+'/../trained_models/RandomForest.rds'
    script: '../scripts/models/RandomForest.R'

rule test_RF:
    """
    tests the Random Forest model using test data
    """
    input:
        RF_model=config['processed_data']+'/../trained_models/RandomForest.rds',
        test_features=config['processed_data']+'/features/SVM_RF/test.rds',
        test_labels=config['processed_data']+'/features/DeepCirCode/DiLiddo2019/y_matrix.txt'
    output:
        prediction=config['processed_data']+'/../evaluation/RandomForest/prediction.txt',
        plot=config['processed_data']+'/../evaluation/RandomForest/roc.jpg'
    script: '../scripts/models/RandomForest_predict.R'


rule train_SVM:
    """
    trains the Support Vector Machine on given data
    """
    input: 
        train_features=config['processed_data']+'/features/SVM_RF/train.rds',
        train_labels=config['processed_data']+'/features/DeepCirCode/Wang2019/y_matrix.txt'
    output:
        SVM_model=config['processed_data']+'/../trained_models/SVM.rds'
    script: '../scripts/models/SVM.R'


rule test_SVM:
    """
    tests the Support Vector Machine model using test data
    """
    input:
        SVM_model=config['processed_data']+'/../trained_models/SVM.rds',
        test_features=config['processed_data']+'/features/SVM_RF/test.rds',
        test_labels=config['processed_data']+'/features/DeepCirCode/DiLiddo2019/y_matrix.txt'
    output:
        prediction=config['processed_data']+'/../evaluation/SVM/prediction.txt',
        plot=config['processed_data']+'/../evaluation/SVM/roc.jpg'
    script: '../scripts/models/SVM_predict.R'


rule performance_assessment:
    """
    perform performance assessment on all predictions
    """
    input:
        RF_prediction=config['processed_data']+'/../evaluation/RandomForest/prediction.txt',
        SVM_prediction=config['processed_data']+'/../evaluation/SVM/prediction.txt'
#       DCC_prediction=config['processed_data']+'/../evaluation/DeepCirCode/prediction.txt',
#       JEDI_prediction=config['processed_data']+'/../evaluation/JEDI/prediction.txt'
    output:
        barplot=config['processed_data']+'/../evaluation/performance.jpg',
        roc=config['processed_data']+'/../evaluation/roc.jpg'
#        prc=config['processed_data']+'/../evaluation/prc.jpg'
    script: '../scripts/evaluation/performance_assessment.R'




rule evaluation:
    """
    Compute evaluation metrics on all model predictions
    Collect all predictions in this rule
    """
    input:
        predictions=expand(
            prediction_pattern, zip, **get_wildcards(params)
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
