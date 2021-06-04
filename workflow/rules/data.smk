# circRNA databases (positive dataset)

rule extract_circBase:
    """
    Retrieve circBase data and save it in a more usable format for training/test data
    """
    # input: # TODO circbase file/download
    output: config['processed_data'] + '/db/circbase.tsv'
    # script: # TODO R or python script
    # shell: # (alternative to script)
    #    """
    #    <shell commands>
    #   """

# TODO add more databases

# Reference files
# TODO: download reference sequence (for features), gene annotation (for negative dataset)


# Create datasets for training and evaluation

rule split_train_test:
    """
    Assemble labels from multiple circRNA & gene annotation datasets
    """
    input: rules.extract_circBase.output[0]
    output:
        train=config['processed_data'] + '/train.tsv',
        test=config['processed_data'] + '/test.tsv'
    script: '../data/split_train_test.py'
