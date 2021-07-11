# Reference files
# Download reference sequence (for features), gene annotation (for negative dataset)

canonical_chroms = [f'chr{i}' for i in list(range(1,23)) + ['X', 'Y']]
assembly = config['assembly']

rule genomepy:
    # Download reference sequence using genomepy
    output:
        protected(multiext(
            config['processed_data'] + "/reference/{assembly}/{assembly}",
            ".fa",".fa.fai",".fa.sizes",".annotation.gtf.gz",".annotation.bed.gz"
        ))
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
    # input:
    #     expand(rules.genomepy.output[3],assembly=config['assembly'])
    output:
        gtf=config['processed_data'] + f"/reference/{assembly}/{assembly}.canonical.gtf",
        exons=config['processed_data'] + f"/reference/{assembly}/{assembly}.exons.gtf",
        transcripts=config['processed_data'] + f"/reference/{assembly}/{assembly}.transcripts.gtf"
    params:
        chroms='|'.join(canonical_chroms),
        url=config['gene_annotations'][assembly]['url'],
        chr_prefix=config['gene_annotations'][assembly]['chr_prefix']
    shell:
        """
        # zcat {input} | grep -E '{params.chroms}' > {output.gtf}
        wget -nc {params.url} -O {output.gtf}.gz
        zcat {output.gtf}.gz \
            | grep -v '#' \
            | sed 's/^/{params.chr_prefix}/' > {output.gtf}.unsorted
        bedtools sort -i {output.gtf}.unsorted > {output.gtf}
        grep -P '\texon\t' {output.gtf} > {output.exons}
        grep -P '\ttranscript\t' {output.gtf} > {output.transcripts}
        """


rule download_chainfile:
    output: config['processed_data'] + '/reference/hg19ToHg38.over.chain'
    params:
        url='ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
    shell:
        """
        wget -nc {params.url} -O {output}.gz
        gunzip {output}.gz
        """


rule get_genome_references:
    # Request fasta file for genome assembly specified in config
    input:
        rules.canonical_gtf.output,
        expand(rules.genomepy.output,assembly=assembly),
        expand(rules.get_phastCons.output,assembly=assembly)


def get_fasta(wildcards):
    return expand(rules.genomepy.output[0],assembly=assembly)


# circRNA junctions (positive dataset)

rule data_DiLiddo2019:
    input: 'resources/high_conf_circ.tab'
    output:
        circRNA=config['processed_data'] + '/datasets/DiLiddo2019/circRNA.bed',
        linear_1=config['processed_data'] + '/datasets/DiLiddo2019/linear1.bed',
        linear_2=config['processed_data'] + '/datasets/DiLiddo2019/linear2.bed'
    run:
        import pandas as pd

        raw_data = pd.read_table(input[0],sep='\t')
        bed_dfs = {}
        strand = '.'

        for col in raw_data:
            df = raw_data[col].str.split(':',expand=True)
            if col == 'circRNA':
                # strand only saved for circRNA BSJ
                df.columns = ['chrom', 'range', 'strand']
                # save strand information for other dfs
                strand = df['strand']
            else:
                df.columns = ['chrom', 'range']
            bed_dfs[col] = df

        for name, bed_df in bed_dfs.items():
            bed_df[['chromStart', 'chromEnd']] = bed_df.iloc[:, 1].str.split('-',expand=True)
            bed_df['name'] = raw_data[name]
            bed_df['score'] = '.'

            if name != 'circRNA':
                # use same strand info as circRNA BSJ
                bed_df['strand'] = strand
                # linear junctions are 1-based -> convert to 0-based
                bed_df['chromStart'] = bed_df['chromStart'].astype(int) - 1

            bed_df = bed_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]
            bed_df.to_csv(output[name],sep='\t',header=False,index=False)


rule data_Chaabane2020:
    """
    Dataset used by Chaabane et al. 2020 (circDeep)
    Processed negative (GENCODE lncRNAs) and positive (circRNADb) datasets, as well as test data
    Subset to canonical chromosomes
    Remove DiLiddo Dataset from positive data
    """
    input:
        positive_bed='methods/circDeep/data/circRNA_dataset.bed',
        negative_bed='methods/circDeep/data/negative_dataset.bed',
        test_bed='methods/circDeep/data/test.bed',
        test_DiLidddo=rules.data_DiLiddo2019.output.circRNA
    output:
        positive_bed=config['processed_data'] + '/datasets/Chaabane2020/circRNA.bed',
        negative_bed=config['processed_data'] + '/datasets/Chaabane2020/negative.bed',
        test_bed=config['processed_data'] + '/datasets/Chaabane2020/test.bed'
    params:
        chroms='|'.join(canonical_chroms)
    run:
        shell(
            """
            grep -E '({params.chroms})\t' {input.positive_bed} > {output.positive_bed}
            grep -E '({params.chroms})\t' {input.negative_bed} > {output.negative_bed}
            grep -E '({params.chroms})\t' {input.test_bed} > {output.test_bed}
            """
        )
        import pandas as pd

        for bedfile in output:
            bed_df = pd.read_table(bedfile,sep='\t',header=None)
            bed_df['score'] = '.'
            bed_df[[0, 1, 2, 4, 'score', 3]].to_csv(bedfile,sep='\t',index=False,header=False)

        shell(
            """
            bedtools subtract -a {output.positive_bed} -b {input.test_DiLidddo} -s -A > tmp_chaab.bed
            mv tmp_chaab.bed {output.positive_bed}
            """
        )


rule data_Wang2019:
    """
    Dataset used by Wang et al. 2019 (DeepCirCode)
    Processed negative (linear junctions in same transcript?) and positive (circRNADb & circBase) datasets
    Remove all high-confidence circRNA test data from DiLiddo2019
    """
    input:
        positive='resources/Wang2019/human_positive.tsv',
        negative='resources/Wang2019/human_negative.tsv',
        test_data=rules.data_DiLiddo2019.output.circRNA,
        chainfile=rules.download_chainfile.output[0]
    output:
        positive=config['processed_data'] + '/datasets/Wang2019/circRNA.bed',
        negative=config['processed_data'] + '/datasets/Wang2019/negative.bed'
    run:
        import pandas as pd
        import pybedtools

        for i, file in enumerate([input.positive, input.negative]):
            df = pd.read_table(
                file,
                skiprows=1,
                usecols=[0, 1, 2, 3, 4],
                names=['chrom', 'chromStart', 'chromEnd', 'strand', 'name'],
                sep='\t'
            )
            df['score'] = '.'
            bed = pybedtools.BedTool.from_dataframe(
                df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]
            )
            if config["assembly"] == "hg38":
                bed = bed.liftover(input.chainfile)

            if i == 0:
                # remove DiLiddo data for positive data
                bed = bed.subtract(input.test_data,s=True,A=True)
            bed.saveas(output[i])


rule data_lncRNA:
    input:
        positive=rules.data_Wang2019.output.positive,
        negative=rules.data_Chaabane2020.output.negative_bed,
    output:
        positive=config['processed_data'] + '/datasets/ncRNA/circRNA.bed',
        negative_train=config['processed_data'] + '/datasets/ncRNA/lncRNA_train.bed',
        negative_test=config['processed_data'] + '/datasets/ncRNA/lncRNA_test.bed',
    run:
        shell("cp {input.positive} {output.positive}")
        # split test and train set
        n_all = shell("wc -l {input.negative} | cut -f1 -d' '",read=True)
        n_train = int(float(n_all) * 0.9)
        shell(
            """
            shuf {input.negative} | split -l {n_train}
            mv xaa {output.negative_train}
            mv xab {output.negative_test}
            """
        )


rule data_NoChr:
    """
    Dataset used by Wang et al. 2019 (DeepCirCode)
    Processed negative (GENCODE lncRNAs) and positive (circRNADb & circBase) datasets
    Remove all high-confidence circRNA test data from DiLiddo2019
    """
    input:
        Chaabane=rules.data_Wang2019.output.positive,
        DiLiddo=rules.data_DiLiddo2019.output.circRNA
    output:
        train=config['processed_data'] + '/datasets/NoChr/train.bed',
        test=config['processed_data'] + '/datasets/NoChr/test.bed'
    shell:
        """
        cat {input.Chaabane} {input.DiLiddo} | grep -v -P 'chr1\t' > {output.train}
        cat {input.Chaabane} {input.DiLiddo} | grep -P 'chr1\t' > {output.test}
        """


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
    elif source == 'Wang2019':
        return rules.data_Wang2019.output.positive
    elif source == 'lncRNA':
        return rules.data_lncRNA.output.positive
    elif source == 'NoChr':
        return rules.data_NoChr.output.train
    elif source == 'NoChr_test':
        return rules.data_NoChr.output.test
    else:
        raise LookupError(f'"{source}" not a valid data source')


# processed negative data

rule get_overlapping:
    input:
        positive=get_positive_data,
        gtf=rules.canonical_gtf.output.gtf
    output:
        genes=config['processed_data'] + '/datasets/{source}/genes_overlapping.gtf',
        transcripts=config['processed_data'] + '/datasets/{source}/transcripts_overlapping.gtf'
    shell:
        """
        grep -P '\tgene\t' {input.gtf} > tmp.gtf
        bedtools intersect -a tmp.gtf -b {input.positive} -wa -s > {output.genes}
        grep -P '\ttranscript\t' {input.gtf} > tmp.gtf
        bedtools intersect -a tmp.gtf -b {input.positive} -wa -s > {output.transcripts}
        rm tmp.gtf
        """


rule get_linear_junctions:
    """
    Get linear junctions subset to circular junctions
    """
    input:
        script='workflow/scripts/data/get_linear_junctions.py',
        gtf=rules.canonical_gtf.output.gtf,
        circ=get_positive_data,
        overlapping_transcripts=rules.get_overlapping.output.transcripts
    output:
        junctions=config['processed_data'] + '/datasets/{source}/junctions_overlapping.bed',
        flanking_junctions=config['processed_data'] + '/datasets/{source}/junctions_flanking.bed'
    threads: 10
    run:
        import pyranges as pr

        gtf = pr.read_gtf(input['gtf'])
        print('read GTF')
        transcripts = pr.read_gtf(input['overlapping_transcripts'])
        print('read transcripts')
        circ = pr.read_bed(input['circ'])
        print('read circRNA BED')

        introns = gtf.features.introns(by='transcript').set_intersect(
            transcripts,strandedness='same',how='containment'
        )
        introns.to_bed(output.junctions)
        print('wrote all junctions')

        flanking = introns.nearest(circ,strandedness='same')
        flanking[flanking.df['Distance'] <= 1].to_bed(output.flanking_junctions)
        print('wrote all flanking junctions')
# shell:
#     """
#     python {input.script} \
#         -gtf {input.gtf} \
#         -circ {input.circ} \
#         -junc {output.junctions} \
#         -flank_junc {output.flanking_junctions}
#     """


def get_negative_data(wildcards, source=None, method=None):
    if source is None:
        try:
            source = wildcards.source
        except:
            raise LookupError("Must define a valid source as wildcard or parameter")
    if method is None:
        try:
            method = wildcards.method
        except:
            raise LookupError("Must define a valid method as wildcard or parameter")

    if source.startswith('lncRNA'):
        if method == 'JEDI':
            target = rules.data_lncRNA.output.negative_train
        else:  # method = in [DeepCirCode, SVM, RandomForest]
            target = rules.get_linear_junctions.output.junctions
    else:
        if method == 'JEDI':
            target = rules.get_overlapping.output.genes
        else:  # method = in [DeepCirCode, SVM, RandomForest]
            target = rules.get_linear_junctions.output.flanking_junctions
    return expand(target,source=source)[0]


# Split datasets for training and evaluation

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
        positive=lambda wildcards: get_positive_data(wildcards,source='Chaabane2020'),
        negative=lambda wildcards: get_negative_data(wildcards,source='Chaabane2020')
    output:
        positive=config['processed_data'] + '/positive_dev.tsv',
        negative=config['processed_data'] + '/negative_dev.tsv'
    shell:
        """
        head {input.positive} > {output.positive}
        head {input.negative} > {output.negative}
        """


def get_test_data(wildcards, source=None):
    if source is None:
        try:
            source = wildcards.source
        except:
            raise LookupError("Must define a valid source as wildcard or parameter")
    if config["dev"]:
        return rules.data_dev.output.positive
    if source == 'DiLiddo2019':
        return rules.data_DiLiddo2019.output.circRNA
    return rules.split_train_test.output.test


rule all_data:
    input:
        positive=lambda w: [get_positive_data(w,source=source) for source in all_sources],
        negative=lambda w: [get_negative_data(w,source=source) for source in all_sources]
