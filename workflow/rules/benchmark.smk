include: "feature_extraction.smk"

tmpdir = config['processed_data'] + '/tmp'
model_pattern = config['models'] + '/{method}/{source}'
evaluation_pattern = config['evaluation'] + '/{method}/{source}'


rule train_circDeep:
    input:
        script='methods/circDeep/circDeep.py',
        bigwig=expand(rules.get_phastCons.output[1],assembly=config['assembly']),
        fasta=expand(rules.genomepy.output[0],assembly=config['assembly']),
        gtf=expand(rules.canonical_gtf.output,assembly=config['assembly']),
        positive=rules.data_Chaabane2020.output.positive_bed,
        negative=rules.data_Chaabane2020.output.negative_bed,
    output:
        model=config['models'] + '/circDeep/{source}/bestmodel_ACNN_BLSTM_3 1 408000.hdf5',
        intermediate_data=directory(tmpdir + '/train/circDeep/{source}')
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
        script='methods/circDeep/circDeep.py',
        bigwig=expand(rules.get_phastCons.output[1],assembly=config['assembly']),
        fasta=expand(rules.genomepy.output[0],assembly=config['assembly']),
        gtf=expand(rules.canonical_gtf.output,assembly=config['assembly']),
        seq_features='methods/circDeep/data/seq_features.txt',
        rcm_features='methods/circDeep/data/rcm_features.txt',
        conservation_features='methods/circDeep/data/conservation_features.txt',
        class_txt='methods/circDeep/data/class.txt',
        test=rules.data_dev.output.positive,
        model=rules.train_circDeep.output.model  # 'methods/circDeep/models'
    output:
        seq_features=tmpdir + '/predict/circDeep/{source}/seq_features.txt',
        rcm_features=tmpdir + '/predict/circDeep/{source}/predict/rcm_features.txt',
        conservation_features=tmpdir + '/predict/circDeep/{source}/predict/conservation_features.txt',
        class_txt=tmpdir + '/predict/circDeep/{source}/predict/class.txt',
        prediction=expand(evaluation_pattern,method="circDeep",allow_missing=True)
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
        train_features=rules.SVM_RF_features.output.train_features,
        train_labels=lambda w: get_train_test(
            w, rules.extract_DeepCirCode_data.output.labels, train_test="train"
        )
    output:
        model=expand(model_pattern + '/model.rds',method='RandomForest',allow_missing=True)[0]
    script: '../scripts/models/RandomForest.R'


rule test_RF:
    """
    tests the Random Forest model using test data
    """
    input:
        model=rules.train_RF.output.model,
        test_features=rules.SVM_RF_features.output.test_features,
        test_labels=lambda w: get_train_test(
            w, rules.extract_DeepCirCode_data.output.labels, train_test="test"
        )
    output:
        prediction=expand(evaluation_pattern + '/prediction.tsv',method='RandomForest',allow_missing=True)[0],
        plot=expand(evaluation_pattern + '/roc.jpg',method='RandomForest',allow_missing=True)[0],
    script: '../scripts/models/RandomForest_predict.R'


rule train_SVM:
    """
    trains the Support Vector Machine on given data
    """
    input:
        train_features=rules.SVM_RF_features.output.train_features,
        train_labels=lambda w: get_train_test(
            w, rules.extract_DeepCirCode_data.output.labels, train_test="train"
        )
    output:
        model=expand(model_pattern + '/model.rds',method='SVM',allow_missing=True)[0]
    script: '../scripts/models/SVM.R'


rule test_SVM:
    """
    tests the Support Vector Machine model using test data
    """
    input:
        model=rules.train_SVM.output.model,
        test_features=rules.SVM_RF_features.output.test_features,
        test_labels=lambda w: get_train_test(
            w, rules.extract_DeepCirCode_data.output.labels, train_test="test"
        )
    output:
        prediction=expand(evaluation_pattern + '/prediction.tsv',method='SVM',allow_missing=True)[0],
        plot=expand(evaluation_pattern + '/roc.jpg',method='SVM',allow_missing=True)[0],
    script: '../scripts/models/SVM_predict.R'


rule train_JEDI:
    """
    Train and predict JEDI on given data
    """
    input:
        script='methods/JEDI/src/run.py',
        data=rules.collect_features_JEDI.input,
        config=rules.extract_data_JEDI.output.config
    output:
        # model=expand(model_pattern + '/model.tf',method='JEDI',allow_missing=True),
        train_eval=expand(evaluation_pattern + '/train_eval.tsv',method='JEDI',allow_missing=True),
        test_eval=expand(evaluation_pattern + '/test_eval.tsv',method='JEDI',allow_missing=True),
        prediction=expand(evaluation_pattern + '/prediction.json',method='JEDI',allow_missing=True)[0]
    params:
        K=config['methods']['JEDI']['kmer_len'],
        L=config['methods']['JEDI']['flank_len'],
        epochs=config['methods']['JEDI']['epochs'],
    conda: '../envs/JEDI.yaml'
    resources:
        gpu=1,
        threads=10
    shell:
        """
        python {input.script} --cv=0 --K={params.K} --L={params.L} \
            --emb_dim=128 --rnn_dim=128 --att_dim=16 --hidden_dim=128 \
            --num_epochs={params.epochs} --learning_rate=1e-3 --l2_reg=1e-3 \
            --config {input.config}
        out_path=$(grep path_pred {input.config} | cut -f2 -d' ')
        mv $out_path/pred.0.K{params.K}.L{params.L} {output.prediction}
        mv $out_path/train_loss.0.K{params.K}.L{params.L} {output.train_eval}
        mv $out_path/test_loss.0.K{params.K}.L{params.L} {output.test_eval}
        """


rule predict_JEDI:
    """
    Translate predictions from JSON to TSV
    """
    input:
        prediction=rules.train_JEDI.output.prediction
    output:
        prediction=expand(evaluation_pattern + '/prediction.tsv',method='JEDI',allow_missing=True)[0]
    run:
        import ujson

        with open(input.prediction,'r') as jpred:
            preds = ujson.load(jpred)

        with open(output.prediction,'w') as out:
            out.write('label\tscore\n')
            for score, label in preds:
                out.write(f'{label}\t{score}\n')


rule performance_assessment:
    """
    perform performance assessment on all predictions
    """
    input:
        RF_prediction=config['processed_data']+'/../evaluation/RandomForest/Wang2019/prediction.tsv',
        SVM_prediction=config['processed_data']+'/../evaluation/SVM/Wang2019/prediction.tsv'
#       DCC_prediction=config['processed_data']+'/../evaluation/DeepCirCode/Wang2019/prediction.tsv',
#       JEDI_prediction=config['processed_data']+'/../evaluation/JEDI/Wang2019/prediction.tsv'
    output:
        barplot=config['processed_data']+'/../evaluation/performance.jpg',
        roc=config['processed_data']+'/../evaluation/roc.jpg',
        pr=config['processed_data']+'/../evaluation/pr.jpg'
    script: '../scripts/evaluation/performance_assessment.R'


rule evaluation:
    """
    Compute evaluation metrics on all model predictions
    Collect all predictions in this rule
    """
    input:
        predictions=expand(evaluation_pattern + '/prediction.tsv',zip,**get_wildcards(params_df))
    params:
        methods=params_df[['method']],
        sources=params_df[['source']]
    output:
        # metrics=config['evaluation'] + '/metrics.tsv',
        barplot = config['evaluation'] + '/performance.jpg',
        roc = config['evaluation'] + '/roc.jpg',
        pr = config['evaluation'] + '/pr.jpg'
    script: '../scripts/evaluation/performance_assessment.R'


rule all_benchmark:
    """
    Collect all output of the benchmark
    """
    input:
        metrics=rules.evaluation.output
