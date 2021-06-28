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
        config['processed_data'] + f"/reference/{config['assembly']}/{config['assembly']}.annotation.canonical.gtf"
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


# circRNA junctions (positive dataset)

rule data_DiLiddo2019:
    input: 'resources/high_conf_circ.tab'
    output:
        circRNA = config['processed_data'] + '/datasets/DiLiddo2019/circRNA.bed',
        linear_1 = config['processed_data'] + '/datasets/DiLiddo2019/linear1.bed',
        linear_2 = config['processed_data'] + '/datasets/DiLiddo2019/linear2.bed'
    run:
        import pandas as pd
        raw_data = pd.read_table(input[0], sep='\t')
        for col in raw_data:
            bed_df = raw_data[col].str.split(':', expand=True)
            if col == 'circRNA':
                bed_df.columns = ['chrom', 'range', 'strand']
            else:
                bed_df.columns = ['chrom', 'range']
                bed_df['strand'] = '.'
            bed_df[['chromStart', 'chromEnd']] = bed_df.iloc[:, 1].str.split('-',expand=True)
            bed_df['name'] = raw_data[col]
            bed_df['score'] = '.'
            bed_df = bed_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]
            bed_df.to_csv(output[col], sep='\t', header=False, index=False)


rule data_Chaabane2020:
    """
    Dataset used by Chaabane et al. 2020 (circDeep)
    Processed negative (GENCODE lncRNAs) and positive (circRNADb) datasets, as well as test data
    Subset to canonical chromosomes
    """
    input:
        positive_bed='methods/circDeep/data/circRNA_dataset.bed',
        negative_bed='methods/circDeep/data/negative_dataset.bed',
        test_bed='methods/circDeep/data/negative_dataset.bed'
    output:
        positive_bed=config['processed_data'] + '/datasets/Chaabane2020/circRNA.bed',
        negative_bed=config['processed_data'] + '/datasets/Chaabane2020/negative.bed',
        test_bed=config['processed_data'] + '/datasets/Chaabane2020/test.bed'
    params:
        chroms='|'.join(canonical_chroms)
    shell:
        """
        grep -E '({params.chroms})\t' {input.positive_bed} > {output.positive_bed}
        grep -E '({params.chroms})\t' {input.negative_bed} > {output.negative_bed}
        grep -E '({params.chroms})\t' {input.test_bed} > {output.test_bed}
        """


rule data_Wang2019:
    """
    Dataset used by Wang et al. 2019 (DeepCirCode)
    Processed negative (GENCODE lncRNAs) and positive (circRNADb & circBase) datasets
    """
    input:
        positive='resources/Wang2019/human_positive.tsv',
        negative='resources/Wang2019/human_negative.tsv'
    output:
        positive=config['processed_data'] + '/datasets/Wang2019/circRNA.bed',
        negative=config['processed_data'] + '/datasets/Wang2019/negative.bed'
    run:
        import pandas as pd
        for i, file in enumerate(input):
            df = pd.read_table(
                file,
                skiprows=1,
                usecols=[0,1,2,3,4],
                names=['chr', 'start', 'end', 'strand', 'name'],
                sep='\t'
            )
            df['score'] = '.'
            df[['chr', 'start', 'end', 'name', 'score', 'strand']].to_csv(
                output[i],
                index=False,
                header=False,
                sep="\t"
            )



def get_positive_data(wildcards, source=None):
    if source is None:
        try:
            source = wildcards.source
        except:
            raise LookupError("Must define a valid source as wildcard or parameter")

    if source == "Chaabane2020":
        return rules.data_Chaabane2020.output.positive_bed
    elif source == 'DiLiddo2019':
        return rules.data_DiLiddo2019.output.circRNA
    else:
        raise LookupError(f'"{source}" not a valid data source')


# processed negative data

rule get_linear_junctions:
    """
    Get linear junctions subset to circular junctions
    """
    input:
        gtf='methods/circDeep/data/circRNA_dataset.bed',
        circ=get_positive_data
    output:
        junctions=config['processed_data'] + '/datasets/gene_annotation/linear_junctions-{source}.bed'
    script: '../scripts/data/get_linear_junctions.py'


def get_negative_data(wildcards, source=None):
    if source is None:
        try:
            source = wildcards.source
        except:
            raise LookupError("Must define a valid source as wildcard or parameter")
    if source == "Chaabane2020":
        return rules.data_Chaabane2020.output.negative_bed
    return expand(rules.get_linear_junctions.output.junctions, source=source)[0]


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
        positive=config['processed_data'] + '/datasets/train_data/positive_{source}.tsv',
        negative=config['processed_data'] + '/datasets/train_data/negative_{source}.tsv',
        test=config['processed_data'] + '/datasets/test_data/{source}.tsv'
    script: '../scripts/data/split_train_test.py'


rule data_dev:
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


def get_test_data(wildcards, source=None):
    if config["dev"]:
        return rules.data_dev.output.positive
    if source == 'DiLiddo2019':
        return rules.data_DiLiddo2019.output.circRNA
    return rules.split_train_test.output.test
