# Reference files
# Download reference sequence (for features), gene annotation (for negative dataset)

canonical_chroms = [f'chr{i}' for i in list(range(1,23)) + ['X', 'Y']]

rule genomepy:
    # Download reference sequence using genomepy
    output:
        multiext(
            config['processed_data'] + "/reference/{assembly}/{assembly}",
            ".fa",".fa.fai",".fa.sizes",".annotation.gtf.gz",".annotation.bed.gz"
        )
    params:
        provider="UCSC"  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
    cache: True  # mark as eligible for between workflow caching
    shell:
        """
        genome_dir=$(dirname "{output[0]}")
        genome_dir=$(dirname "$genome_dir")
        echo downloading genome to: $genome_dir
        genomepy install \
            {wildcards.assembly} \
            -p {params.provider} \
            --annotation \
            -g $genome_dir
        """


rule get_phastCons:
    # Request fasta file for genome assembly specified in config
    output:
        multiext(
            config['processed_data'] + '/reference/{assembly}/phastCons20way/{assembly}.phastCons20way',
            '.wigFix.gz','.bw'
        )
    params:
        url='rsync://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/phastCons20way/{assembly}.phastCons20way.wigFix.gz',
        chrom_sizes='http://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.chrom.sizes'
    shell:
        """
        rsync -avz --progress {params.url} {output[0]}
        wigToBigWig {output[0]} {params.chrom_sizes} {output[1]}
        """


rule canonical_gtf:
    input:
        expand(rules.genomepy.output[3],assembly=config['assembly'])
    output:
        config['processed_data'] + "/reference/{assembly}/{assembly}.annotation.canonical.gtf"
    params:
        chroms='|'.join(canonical_chroms)
    shell:
        """
        zcat {input} | grep -E '{params.chroms}' > {output}
        """


rule get_genome_references:
    # Request fasta file for genome assembly specified in config
    input:
        expand(rules.genomepy.output,assembly=config['assembly']),
        expand(rules.get_phastCons.output,assembly=config['assembly'])


# circRNA databases (positive dataset)

rule extract_circBase:
    """
    Retrieve circBase data and save it in a more usable format for training/test data
    """
    # input: # TODO circbase file/download
    output: config['processed_data'] + '/db/circbase.tsv'
    # script: # TODO R or python script


rule data_Chaabane2020:
    """
    Dataset used by Chaabane et al. 2020 (circDeep)
    processed negative and positive datasets, as well as test data
    BED file format
    """
    input:
        positive_bed='methods/circDeep/data/circRNA_dataset.bed',
        negative_bed='methods/circDeep/data/negative_dataset.bed',
        test_bed='methods/circDeep/data/negative_dataset.bed'
    output:
        positive_bed=config['processed_data'] + '/datasets/Chaabane2020/circRNA_dataset.bed',
        negative_bed=config['processed_data'] + '/datasets/Chaabane2020/negative_dataset.bed',
        test_bed=config['processed_data'] + '/datasets/Chaabane2020/test_dataset.bed'
    params:
        chroms='|'.join(canonical_chroms)
    shell:
        """
        grep -E '({params.chroms})\t' {input.positive_bed} > {output.positive_bed}
        grep -E '({params.chroms})\t' {input.negative_bed} > {output.negative_bed}
        grep -E '({params.chroms})\t' {input.test_bed} > {output.test_bed}
        """


# Getters for positive and negative datasets

def get_positive_data(wildcards, source=None):
    if source == "Chaabane2020":
        return rules.data_Chaabane2020.output.positive_bed
    return rules.extract_circBase.output[0]


def get_negative_data(wildcards, source=None):
    if source == "Chaabane2020":
        return rules.data_Chaabane2020.output.negative_bed
    return rules.canonical_gtf.output[0]


# Create datasets for training and evaluation

rule split_train_test:
    """
    Assemble labels from multiple circRNA & gene annotation datasets
    Split train and test set, add validation set annotations
    """
    input:
        positive=get_positive_data,
        negative=get_negative_data
    output:
        train=config['processed_data'] + '/train.tsv',
        test=config['processed_data'] + '/test.tsv'
    script: '../data/split_train_test.py'


rule subset_data:
    """
    Subset data purely for development
    """
    input:
        positive=lambda wildcards: get_positive_data(wildcards, source='Chaabane2020'),
        negative=lambda wildcards: get_negative_data(wildcards, source='Chaabane2020')
    output:
        positive=config['processed_data'] + '/positive_dev.tsv',
        negative=config['processed_data'] + '/negative_dev.tsv'
    shell:
        """
        head {input.positive} > {output.positive}
        head {input.negative} > {output.negative}
        """


def get_test_data(wildcards):
    if config["dev"]:
        return rules.subset_data.output.positive
    return rules.split_train_test.output.test
