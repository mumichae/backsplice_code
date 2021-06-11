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
            '.wigFix.gz', '.bw'
        )
    params:
        url = 'rsync://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/phastCons20way/{assembly}.phastCons20way.wigFix.gz',
        chrom_sizes= 'http://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.chrom.sizes'
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
        chroms = '|'.join(canonical_chroms)
    shell:
        """
        zcat {input} | grep -E '{params.chroms}' > {output}
        """


rule get_genome_references:
    # Request fasta file for genome assembly specified in config
    input:
        expand(rules.genomepy.output,assembly=config['assembly']),
        expand(rules.get_phastCons.output, assembly=config['assembly'])


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


rule Chaabane2020:
    input:
        positive_bed = 'methods/circDeep/data/circRNA_dataset.bed',
        negative_bed = 'methods/circDeep/data/negative_dataset.bed'
    output:
        positive_bed = config['processed_data'] + '/datasets/Chaabane2020/circRNA_dataset.bed',
        negative_bed = config['processed_data'] + '/datasets/Chaabane2020/negative_dataset.bed'
    params:
        chroms='|'.join(canonical_chroms)
    shell:
        """
        grep -E '({params.chroms})\t' {input.positive_bed} > {output.positive_bed}
        grep -E '({params.chroms})\t' {input.negative_bed} > {output.negative_bed}
        """


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
