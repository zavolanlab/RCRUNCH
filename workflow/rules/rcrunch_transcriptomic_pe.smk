# --------------------------------------------------------------------------------
# RCRUNCH
# Author : Katsantoni Maria
# Company: Mihaela Zavolan group, Biozentrum, Basel
# --------------------------------------------------------------------------------
# Transcriptomic RCRUNCH subpipeline for paired-end CLIP data
# This approach is based on the hypothesis that the RBP is binding the mature
# mRNA and thus the sequence of the mature mRNA is important for the binding.
# ________________________________________________________________________________
import os

rule TR_initial_map_to_genome:
    """
        Initial mapping to the genome using STAR.
    """
    input:
        genome = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),
        reads1 = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "cutadapt",
                    "{sample}_mate1.cutadapt_{format}.fastq.gz"),
                sample=wildcards.sample,
                format=config[wildcards.sample]['format']
                ),
        reads2 = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "cutadapt",
                    "{sample}_mate2.cutadapt_{format}.fastq.gz"),
                sample=wildcards.sample,
                format=config[wildcards.sample]['format']
                ),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}",
                "{sample}.bam")
                ),
        logfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}",
                "{sample}Log.final.out")
                ),
        outdir = temp(
            directory(
                os.path.join(
                    config["output_dir"],
                    "TR",
                    "initial_map",
                    "{sample}")
                    )
                ),

    params:
        cluster_log_path = config["cluster_log"],
        sample_id = "{sample}",
        genome = os.path.join(
            config["output_dir"],
            "STAR_index"),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}",
            "{sample}"),
        multimappers = config["multimappers"]

    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_map_to_genome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_map_to_genome.stderr.log"),

    shell:
        "(STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.genome} \
        --readFilesIn {input.reads1} {input.reads2} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 1 \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outSAMattributes All \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outFilterScoreMinOverLread 0.9 \
        --outFilterMatchNminOverLread 0.9 \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.sample_id} \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_initial_index:
    """
        Index genome bamfile using samtools.
    """
    input:
        outdir = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}"),
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}",
            "{sample}.bam"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}.sorted.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}.sorted.bam.bai")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}_temp"),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 1

    log:
       stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_index.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_index.stderr.log"),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai}) \
         1> {log.stdout} 2> {log.stderr}"


rule TR_flag_ncRNA_reads:
    """
        Make txt file of reads mapping to ncRNAs.
        If a multimapper maps to a ncRNA, all reads are excluded.
        Custom script.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}.sorted.bam"),
        bai = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}.sorted.bam.bai"),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_filter_ncRNAs.py"),
        ncRNAs = config['ncRNAs'],
        ncRNA_biotypes = expand(config['ncRNA_biotypes']),
        paired = lambda wildcards:
            config[wildcards.sample]['paired'],
        sense = lambda wildcards:
            config[wildcards.sample]['sense'],

    output:
        outfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_ncRNAs",
                "{sample}.ncrna.txt")
                ),
        flag = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_ncRNAs",
                "{sample}_done")
                ),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_ncRNA_reads.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_ncRNA_reads.stderr.log"),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --RNA_central {params.ncRNAs} \
        --outfile {output.outfile} \
        --paired {params.paired} \
        --sense {params.sense} \
        --ncRNAs {params.ncRNA_biotypes} \
        --flag {output.flag} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_remove_ncRNA_reads:
    """
        Remove the reads mapping to rRNAs using picard.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "initial_map",
            "{sample}.sorted.bam"),
        read_names = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.ncrna.txt"),

    output:
        reads = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_ncRNAs",
                "{sample}.filtered.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/picard:2.18.9"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_remove_ncRNA_reads.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_remove_ncRNA_reads.stderr.log"),

    shell:
        "(java -jar /usr/local/bin/picard.jar FilterSamReads \
        I={input.bam} \
        O={output.reads} \
        READ_LIST_FILE={input.read_names} \
        FILTER=excludeReadList) \
        1> {log.stdout} 2> {log.stderr}"


rule TR_index_ncRNA_rm:
    """
        Index bamfile without the rRNA-mapped reads using samtools.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam"),

    output:
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_ncRNAs",
                "{sample}.filtered.bam.bai")
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_index_ncRNA_rm.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_index_ncRNA_rm.stderr.log"),

    shell:
        "(samtools index \
        {input.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_flag_duplicates:
    """
        Flag duplicate reads using the STAR flags.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam"),
        bai = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam.bai"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "flag_duplicates",
                "{sample}.duplicates.Processed.out.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates."),

    singularity:
        "docker://zavolab/star:2.6.0a"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_duplicates.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_duplicates.stderr.log"),

    shell:
        "(STAR \
        --inputBAMfile {input.bam} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --runMode inputAlignmentsFromBAM \
        --outFileNamePrefix {params.outFileNamePrefix} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_remove_duplicates:
    """
        Remove reads flagged as duplicates by STAR.
        Custom script.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.Processed.out.bam"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.duplicates.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_filter_duplicates.py"),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_remove_duplicates.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_remove_duplicates.stderr.log"),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --outfile {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_no_duplicate_removal:
    """
        No removal of duplicates
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam"),
        bai = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam.bai"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.with_duplicates.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.with_duplicates.bam.bai")
                ),
    
    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://bash:5.0.16"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_no_duplicate_removal.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_no_duplicate_removal.stderr.log"),

    shell:
        "(cp {input.bam} {output.bam}; \
        cp {input.bai} {output.bai}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_umi_collapse:
    """
        Alternative duplicate removal when UMIs are present using umitools.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam"),
        bai = os.path.join(
            config["output_dir"],
            "TR",
            "remove_ncRNAs",
            "{sample}.filtered.bam.bai"),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "TR",
            "remove_duplicates",
            "{sample}.umis.bam")
            ),

    params:
        cluster_log_path = config["cluster_log"]

    singularity:
        "docker://zavolab/umi-tools:0.5.4"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_umi_collapse.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_umi_collapse.stderr.log"),

    shell:
        "(umi_tools \
        dedup \
        -I {input.bam} \
        --paired \
        -S {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_merge:
    input:
        bam = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "TR",
                    "remove_duplicates",
                    "{sample}.{dup_type}.bam"),
                sample=config[wildcards.experiment][wildcards.name],
                dup_type=config[config[wildcards.experiment][wildcards.name][0]]['dup_type'],
                ),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_merge_samples.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_merge.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_merge.stderr.log"),

    shell:
        "(samtools merge \
        {output.bam} {input.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_sort_merged:
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_merge_samples.bam"
            ),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_merge_samples.sorted.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_merge_samples.sorted.bam.bai")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_merge_samples"),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_sort_merged.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_sort_merged.stderr.log"),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_sortname_merged:
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_merge_samples.bam"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_merge_samples.sortedname.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_merge_samples_sortname"),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_sortname_merged.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_sortname_merged.stderr.log"),

    shell:
        "(samtools sort -n \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_read_frequencies:
    """
        Calculate occurence of reads that are multimappers.
        Custom script.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam"),
        bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam.bai"),

    output:
        frequencies = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.frequencies.csv")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_bam_get_read_frequencies.py"),
        paired = lambda wildcards:
                    config[config[wildcards.experiment][wildcards.name][0]]['paired'],

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_read_frequencies.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_read_frequencies.stderr.log"),
        
    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --paired {params.paired} \
        --outfile {output.frequencies}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_fastq:
    """
        Convert initially mapped reads back to fastq files to do mappings at the
        transcriptome and genome level, using bedtools.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sortedname.bam"),
        bam_coord = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam"),

    output:
        fastq_mate1 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate1.unfiltered.fastq")
                ),
        fastq_mate2 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate2.unfiltered.fastq")
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/bedtools:2.27.0"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq.stderr.log"),

    shell:
        "(bedtools bamtofastq \
        -i {input.bam} \
        -fq {output.fastq_mate1} \
        -fq2 {output.fastq_mate2}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_fastq_filtered:
    """
        Remove duplicate reads after converting bam to fastq
        (when there are multimappers the reads will appear
        multiple times in the fastq).
    """
    input:
        mate1 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate1.unfiltered.fastq"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.unfiltered.fastq"),

    output:
        mate1 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate1.fastq.gz")
                ),
        mate2 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate2.fastq.gz")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_fastq_keep_unique_reads.py"),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq_filtered.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq_filtered.stderr.log"),

    shell:
        "(python {params.script} \
        --fastq_in {input.mate1} \
        --fastq_out {output.mate1} \
        --fastq_in2 {input.mate2} \
        --fastq_out2 {output.mate2}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_salmon_index:
    """
        Build salmon index based on transcript annotation.
        Required for transcript quantification done by salmon.
    """
    input:
        transcriptome = config["transcriptome"]

    output:
        index = directory(os.path.join(
            config["output_dir"],
            "TR_salmon_index")
            ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/salmon:0.11.0"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "salmon_index.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "salmon_index.stderr.log"),

    shell:
        "(salmon index \
        -t {input.transcriptome} \
        -i {output.index} \
        --type fmd; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_tr_quantification:
    """
        Quantification of transcript isoforms using Salmon.
    """
    input:
        transcriptome = os.path.join(
            config["output_dir"],
            "TR_salmon_index"),
        mate1 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate1.fastq.gz"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.fastq.gz"),

    output:
        quantification = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_quantification_{name}",
                "quant.sf")
                ),
        outdir = temp(
            directory(
                os.path.join(
                    config["output_dir"],
                    "TR",
                    "{experiment}",
                    "tr_quantification_{name}")
                    )
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/salmon:0.11.0"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_tr_quantification.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_tr_quantification.stderr.log"),

    shell:
        "(salmon quant \
        -i {input.transcriptome} \
        -l ISR \
        -1 {input.mate1} \
        -2 {input.mate2} \
        -o {output.outdir} \
        -p {threads} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_customised_transcriptome:
    """
        Build customised transcriptome keeping
        the most expressed isoform per gene.
        Custom script.
    """
    input:
        transcriptome = config["transcriptome"],
        ensembl_csv = os.path.join(
            config["output_dir"],
            "ensembl_csv",
            "ensembl.csv"),
        ensembl_gtf = config["gtf"],
        quantification = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_quantification_{name}",
            "quant.sf"),
        outdir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_quantification_{name}"),

    output:
        output_gtf = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.transcriptome.gtf"),
        output_fasta = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.transcriptome.fa"),
        info = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.info.csv"),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_most_expressed_transcriptome.py"),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_customised_transcriptome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_customised_transcriptome.stderr.log"),

    shell:
        "(python {params.script} \
        --quantification_file {input.quantification} \
        --transcriptome {input.transcriptome}  \
        --ensembl_csv {input.ensembl_csv} \
        --ensembl_gtf {input.ensembl_gtf} \
        --out_gtf {output.output_gtf} \
        --transcript_info {output.info} \
        --out_fasta {output.output_fasta} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_transcriptome_index:
    """
        Create index of the customised transcriptome using STAR.
    """
    input:
        tr_fasta = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.transcriptome.fa"),
        tr_gtf = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.transcriptome.gtf"),

    output:
        output = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_STAR_index",
            "chrNameLength.txt"),

    params:
        cluster_log_path = config["cluster_log"],
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_STAR_index"),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_STAR_index/STAR_"),
        sjdbOverhang = config["sjdbOverhang"]

    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_transcriptome_index.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_transcriptome_index.stderr.log"),

    shell:
        "(chmod -R 777 {params.output_dir}; \
        STAR \
        --runMode genomeGenerate \
        --genomeSAindexNbases 12 \
        --genomeChrBinNbits 12 \
        --genomeDir {params.output_dir} \
        --genomeFastaFiles {input.tr_fasta} \
        --runThreadN {threads} \
        --outFileNamePrefix {params.outFileNamePrefix} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_map_to_transcriptome:
    """
        Map to the transcriptome treating each transcript
        as a "chromosome", using STAR.
    """
    input:
        transcriptome = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),
        mate1 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate1.fastq.gz"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.fastq.gz"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map",
                "{name}",
                "transcriptome.bam")
                ),
        logfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map",
                "{name}",
                "transcriptomeLog.final.out")
                ),
        outdir = temp(
            directory(
                os.path.join(
                    config["output_dir"],
                    "TR",
                    "{experiment}",
                    "tr_map",
                    "{name}")
                    )
                ),
    params:
        cluster_log_path = config["cluster_log"],
        transcriptome = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index"),
        experiment_id = "{experiment}_{name}_transcriptome",
        fragment_size = config["fragment_size"],
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_map",
            "{name}",
            "transcriptome"),
        multimappers = config['multimappers'] + 10,

    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_transcriptome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_transcriptome.stderr.log"),

    shell:
        "(STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.transcriptome} \
        --readFilesIn {input.mate1} {input.mate2} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 1 \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outFilterScoreMinOverLread 0.9 \
        --outFilterMatchNminOverLread 0.9 \
        --outFilterMismatchNoverLmax 0.05 \
        --outSAMattributes All \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.experiment_id} \
        --alignIntronMax 1 \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_index_map_to_transcriptome:
    """
        Index the transcriptomic alignment.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_map",
            "{name}",
            "transcriptome.bam"
            ),
        outdir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_map",
            "{name}")

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map_sort",
                "{name}",
                "transcriptome.sorted.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map_sort",
                "{name}",
                "transcriptome.sorted.bam.bai")
                ),
    
    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_map",
            "{name}",
            "transcriptome_temp"),

    singularity:
        "docker://zavolab/samtools:1.8"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_index_map_to_transcriptome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_index_map_to_transcriptome.stderr.log"),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_map_to_genome:
    """
        Initial mapping to the genome using STAR.
    """
    input:
        genome = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),
        mate1 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate1.fastq.gz"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.fastq.gz"),

    output:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome.bam"),
        logfile = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genomeLog.final.out"),
        outdir = temp(
            directory(
                os.path.join(
                    config["output_dir"],
                    "TR",
                    "{experiment}",
                    "gn_map",
                    "{name}")
                    )
                ),

    params:
        cluster_log_path = config["cluster_log"],
        experiment_id = "{experiment}_{name}_genome",
        genome = os.path.join(
            config["output_dir"],
            "STAR_index"),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome"),
        multimappers = config["multimappers"]

    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_genome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_genome.stderr.log"),

    shell:
        "(STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.genome} \
        --readFilesIn {input.mate1} {input.mate2} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMinOverLread 0.9 \
        --outFilterMatchNminOverLread 0.9 \
        --outFilterMismatchNoverLmax 0.05 \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outSAMattributes All \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.experiment_id} \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_index_map_to_genome:
    """
        Index the genomic alignment.
    """
    input:
        outdir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}"),
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome.bam"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "gn_map_sort",
                "{name}",
                "genome.sorted.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "gn_map_sort",
                "{name}",
                "genome.sorted.bam.bai")
                ),
    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map_sort",
            "{name}",
            "genome_temp"),

    singularity:
        "docker://zavolab/samtools:1.8"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_index_map_to_genome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_index_map_to_genome.stderr.log"),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_assign_to_gn_tr:
    """
        Remove reads that map to '-' strand in transcriptome
    """
    input:
        ensembl_csv = os.path.join(
            config["output_dir"],
            "ensembl_csv",
            "ensembl.csv"),
        chromosome_info = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),
        gn_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map_sort",
            "{name}",
            "genome.sorted.bam"),
        gn_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map_sort",
            "{name}",
            "genome.sorted.bam.bai"),
        tr_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_map_sort",
            "{name}",
            "transcriptome.sorted.bam"),
        tr_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_map_sort",
            "{name}",
            "transcriptome.sorted.bam.bai"),

    output:
        tr_bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_transcriptome.bam")
                ),
        gn_bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_genome.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_assign_read_to_gn_tr.py"),
        out_folder = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_assign_to_gn_tr_temp"),
        paired = lambda wildcards:
                    config[config[wildcards.experiment][wildcards.name][0]]['paired'],
        sense = lambda wildcards:
                    config[config[wildcards.experiment][wildcards.name][0]]['sense'],

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_assign_to_gn_tr.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_assign_to_gn_tr.stderr.log"),

    shell:
        "(mkdir -p {params.out_folder} ;\
        python {params.script} \
        --annotation {input.ensembl_csv}  \
        --tr_bam {input.tr_bam} \
        --gn_bam {input.gn_bam} \
        --out_folder {params.out_folder} \
        --filtered_tr_bam {output.tr_bam} \
        --filtered_gn_bam {output.gn_bam} \
        --paired {params.paired} \
        --sense {params.sense} \
        --chromosome_info {input.chromosome_info} \
        --threads {threads} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_transcriptome_filter_sort:
    """
        Sort and index transcriptomic reads after
        filtering wrong strand reads and
        those that map better to genome
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_transcriptome.bam"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_transcriptome.sorted.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_transcriptome.sorted.bam.bai")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_sort_samples"),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_transcriptome_filter_sort.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_transcriptome_filter_sort.stderr.log"),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_genome_filter_sort:
    """
        Sort and index genomic reads after filtering wrong strand reads and
        those that map better to genome
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_genome.bam"),

    output:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_genome.sorted.bam"),
        bai =  os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_genome.sorted.bam.bai"),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{name}_gn_sort_samples"),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_genome_filter_sort.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_genome_filter_sort.stderr.log"),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_bed:
    """
        Convert bamfile to bedfile using bedtools
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_genome.sorted.bam"),

    output:
        bed = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_genome.bed")
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/bedtools:2.27.0"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_bed.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_bed.stderr.log"),

    shell:
        "(bedtools bamtobed \
        -i {input.bam} > {output.bed}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_bed_transcriptome:
    """
        Convert bamfile to bedfile using bedtools
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}_transcriptome.sorted.bam"),

    output:
        bed = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}_transcriptome.bed")
                ),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/bedtools:2.27.0"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_bed_transcriptome.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_bed_transcriptome.stderr.log"),

    shell:
        "(bedtools bamtobed \
        -i {input.bam} > {output.bed}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_genome_coverage:
    """
        Count number of fragments in windows 
        of specific size at the genome.
        Custom script.
    """
    input:
        chromosome_info = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),
        fg_bed = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_genome.bed"),
        bg_bed = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_genome.bed"),
        fg_frequencies = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.frequencies.csv"),
        bg_frequencies = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis.frequencies.csv"),

    output:
        coverage = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}.genome_coverage.bed")
                ),
    
    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "sliding_windows.py"),
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}"),
        window_f = lambda wildcards:
            config[wildcards.experiment]["window_f"],
        window_b = lambda wildcards:
            config[wildcards.experiment]["window_b"],
        step_size = lambda wildcards:
            config[wildcards.experiment]["step_size"],
        sense_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['sense'],
        sense_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['sense'],
        paired_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['paired'],
        paired_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['paired'],
        prefix = "{experiment}.genome_coverage",
        background_type = lambda wildcards:
            config[wildcards.experiment]["background_type"],

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "genome_coverage.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "genome_coverage.stderr.log"),

    shell:
        "(python {params.script} \
        --bed_foreground {input.fg_bed} \
        --bed_background {input.bg_bed} \
        --foreground_frequencies {input.fg_frequencies} \
        --background_frequencies {input.bg_frequencies} \
        --out {params.output_dir} \
        --prefix {params.prefix} \
        --chromosomes {input.chromosome_info} \
        --sense_f {params.sense_f} \
        --sense_b {params.sense_b} \
        --paired_f {params.paired_f} \
        --paired_b {params.paired_b} \
        --window_f {params.window_f} \
        --window_b {params.window_b} \
        --step_size {params.step_size} \
        --data_type genome \
        --background {params.background_type} \
        --cutoff 1 \
        --threads {threads}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_transcriptome_coverage:
    """
        Count number of fragments in windows 
        of specific size at
        the customised transcriptome.
        Custom script.
    """
    input:
        chromosome_info = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),
        fg_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_transcriptome.sorted.bam"),
        fg_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_transcriptome.sorted.bam.bai"),
        bg_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_transcriptome.sorted.bam"),
        bg_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_transcriptome.sorted.bam.bai"),
        fg_frequencies = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.frequencies.csv"),
        bg_frequencies = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis.frequencies.csv"),

    output:
        coverage = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}.transcriptome_coverage.bed")
                ),
            
    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_bam_count_sliding_windows.py"),
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}"),
        window_f = lambda wildcards:
            config[wildcards.experiment]["window_f"],
        window_b = lambda wildcards:
            config[wildcards.experiment]["window_b"],
        step_size = lambda wildcards:
            config[wildcards.experiment]["step_size"],
        sense_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['sense'],
        sense_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['sense'],
        paired_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['paired'],
        paired_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['paired'],
        prefix = "{experiment}.transcriptome_coverage",
        background_type = lambda wildcards:
            config[wildcards.experiment]["background_type"]

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "transcriptome_coverage.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "transcriptome_coverage.stderr.log"),

    shell:
        "(python {params.script} \
        --bam_foreground {input.fg_bam} \
        --bam_background {input.bg_bam} \
        --bam_foreground_frequencies {input.fg_frequencies} \
        --bam_background_frequencies {input.bg_frequencies} \
        --out {params.output_dir} \
        --prefix {params.prefix} \
        --chromosomes {input.chromosome_info} \
        --sense_f {params.sense_f} \
        --sense_b {params.sense_b} \
        --paired_f {params.paired_f} \
        --paired_b {params.paired_b} \
        --window_f {params.window_f} \
        --window_b {params.window_b} \
        --step_size {params.step_size} \
        --data_type transcriptome \
        --background {params.background_type} \
        --cutoff 1 \
        --threads {threads}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_merge_windows:
    """
        Merge windows from foreground/ background.
        Custom script.
    """
    input:
        tr_coverage = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}.transcriptome_coverage.bed"),
        gn_coverage = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}.genome_coverage.bed"),

    output:
        merged_coverage = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}.total_coverage.bed")
            ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_merge_fg_bg.py"),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "merge_windows.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "merge_windows.stderr.log"),

    shell:
        "(python {params.script} \
        --transcriptome {input.tr_coverage} \
        --genome {input.gn_coverage} \
        --out {output.merged_coverage} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_enriched_regions:
    """
        Obtain enriched regions for binding based on EM using
        bg for noise estimation. Modified and adapted from:
        Berger, Severin M., et al. "Crunch: Integrated processing and modeling of
        ChIP-seq data in terms of regulatory motifs."
        Genome research (2019): gr-239319
    Custom script.
   """
    input:
        coverage = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}.total_coverage.bed"),

    output:
        enriched_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_enriched_regions.csv"),
        parameters = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_parameters.csv"),
        all_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_regions.csv"),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_find_enriched_regions.py"),
        plot_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "enriched_regions_plots"),
        fdr_cutoff = config['FDR']

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "enriched_regions.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "enriched_regions.stderr.log"),

    shell:
        "(mkdir -p {params.plot_dir}; \
        python {params.script} \
        --windows_table {input.coverage} \
        --out {output.enriched_regions} \
        --all_regions {output.all_regions} \
        --out_parameters {output.parameters} \
        --out_folder {params.plot_dir} \
        --fdr_cutoff {params.fdr_cutoff} \
         ) 1> {log.stdout} 2> {log.stderr}"


rule TR_fit_peaks:
    """
        Fit individual binding events in the enriched regions.
        Modifed and adapted from:
        Berger, Severin M., et al. "Crunch: Integrated processing and
        modeling of ChIP-seq data in terms of regulatory motifs."
        Genome research (2019): gr-239319.
        Custom script.
    """
    input:
        fg_gn_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_genome.sorted.bam"),
        fg_gn_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_genome.sorted.bam.bai"),
        bg_gn_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_genome.sorted.bam"),
        bg_gn_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_genome.sorted.bam.bai"),
        fg_tr_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_transcriptome.sorted.bam"),
        fg_tr_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates_transcriptome.sorted.bam.bai"),
        bg_tr_bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_transcriptome.sorted.bam"),
        bg_tr_bai = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis_transcriptome.sorted.bam.bai"),
        fg_frequencies = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.frequencies.csv"),
        bg_frequencies = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_smis.frequencies.csv"),
        enriched_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_enriched_regions.csv"),
        parameters = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_parameters.csv"),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),
        chromosomes_t = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),

    output:
        peaks = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_peaks.csv"),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_define_peaks.py"),
        plot_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "fit_peaks_plots"),
        fragment_size = config["fragment_size"],
        window_size = lambda wildcards:
            config[wildcards.experiment]['window_f'],
        paired_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['paired'],
        paired_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['paired'],
        sense_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['sense'],
        sense_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['sense'],

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "fit_peaks.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "fit_peaks.stderr.log"),

    shell:
        "(mkdir -p {params.plot_dir}; \
        python {params.script} \
        --enriched_regions {input.enriched_regions} \
        --parameters {input.parameters} \
        --fragment_size {params.fragment_size} \
        --window_size {params.window_size} \
        --chromosomes_g {input.chromosomes_g} \
        --bam_foreground_g {input.fg_gn_bam} \
        --bam_background_g {input.bg_gn_bam} \
        --chromosomes_t {input.chromosomes_t} \
        --bam_foreground_t {input.fg_tr_bam} \
        --bam_background_t {input.bg_tr_bam} \
        --bam_fg_fq {input.fg_frequencies} \
        --bam_bg_fq {input.bg_frequencies} \
        --paired_f {params.paired_f} \
        --paired_b {params.paired_b} \
        --sense_f {params.sense_f} \
        --sense_b {params.sense_b} \
        --out {output.peaks} \
        --out_folder {params.plot_dir} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_peaks_split_sets_crosslink:
    """
        Split the significant peaks into training and test
        for subsequent motif analysis
    """
    input:
        genome_fasta = config["genome"],
        transcriptome_fasta = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.transcriptome.fa"),
        all_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_regions.csv"),
        reads = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_peaks.csv"),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),
        chromosomes_t = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),

    output:
        training = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "crosslink",
                "training.fasta")
            ),
        test = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "crosslink",
                "test.fasta")
                ),
        training_bg = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "crosslink",
                "training_bg.fasta")
                ),
        test_bg = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "crosslink",
                "test_bg.fasta")
                ),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "crosslink",
                "tmpdir")
            )
        ),


    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_peaks_split_sets.py"),
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "crosslink"),
        genome_tag = config["genome_tag"],
        peak_size = 50

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}_peaks_split_sets_crosslink.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}_peaks_split_sets_crosslink.stderr.log"),

    shell:
        "(mkdir -p {output.dirtemp}; \
        python {params.script} \
        --peaks_file {input.reads} \
        --genome_fasta {input.genome_fasta} \
        --all_regions {input.all_regions} \
        --genome_tag {params.genome_tag} \
        --tempdir {output.dirtemp} \
        --transcriptome_fasta {input.transcriptome_fasta} \
        --chromosomes_g {input.chromosomes_g} \
        --chromosomes_t {input.chromosomes_t} \
        --crosslink_type 'crosslink' \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_peaks_split_sets:
    """
        Split the significant peaks into training and test
        for subsequent motif analysis
    """
    input:
        genome_fasta = config["genome"],
        transcriptome_fasta = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.transcriptome.fa"),
        all_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_regions.csv"),
        reads = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_peaks.csv"),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),
        chromosomes_t = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),

    output:
        training = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "peak_center",
                "training.fasta")
                ),
        test = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "peak_center",
                "test.fasta")
                ),
        training_bg = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "peak_center",
                "training_bg.fasta")
                ),
        test_bg = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "peak_center",
                "test_bg.fasta")
                ),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_crossvalidation",
                "{run}",
                "peak_center",
                "tmpdir")
            )
        ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_peaks_split_sets.py"),
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "peak_center"),
        genome_tag = config["genome_tag"],
        peak_size = 50

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}_peaks_split_sets_peak_center.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}_peaks_split_sets_peak_center.stderr.log"),

    shell:
        "(mkdir -p {output.dirtemp}; \
        python {params.script} \
        --peaks_file {input.reads} \
        --genome_fasta {input.genome_fasta} \
        --all_regions {input.all_regions} \
        --tempdir {output.dirtemp} \
        --genome_tag {params.genome_tag} \
        --transcriptome_fasta {input.transcriptome_fasta} \
        --chromosomes_g {input.chromosomes_g} \
        --chromosomes_t {input.chromosomes_t} \
        --crosslink_type 'peak_center' \
        --out_folder {params.output_dir} \
        ) 1> {log.stdout} 2> {log.stderr}"

rule TR_peaks_split_sets_crosslink_all_motifs:
    """
        Get all peaks crosslink
    """
    input:
        transcriptome_fasta = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.transcriptome.fa"),
        genome_fasta = config["genome"],
        all_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),
        chromosomes_t = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),

    output:
        test = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_final",
                "crosslink",
                "test.fasta")
            ),
        test_bg = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_final",
                "crosslink",
                "test_bg.fasta")
            ),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_final",
                "crosslink",
                "tmpdir")
            )
        ),
    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "total_enrichment_peakstofasta.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "motif_analysis_final",
            "crosslink",
            ),
        genome_tag = config["genome_tag"],
        peak_size = 50,
        
    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_final",
            "peaks_split_sets__all_motifs_crosslink.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_final",
            "peaks_split_sets__all_motifs_crosslink.stderr.log"
            ),

    shell:
        "(mkdir -p {params.output_dir}; \
        mkdir -p {output.dirtemp}; \
        python {params.script} \
        --peaks_file {input.reads} \
        --genome_fasta {input.genome_fasta} \
        --all_regions {input.all_regions} \
        --genome_tag {params.genome_tag} \
        --tempdir {output.dirtemp} \
        --transcriptome_fasta {input.transcriptome_fasta} \
        --chromosomes_g {input.chromosomes_g} \
        --chromosomes_t {input.chromosomes_t} \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        --crosslink_type 'crosslink' \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_peaks_split_sets_all_motifs:
    """
        get all peaks crosslink
    """
    input:
        transcriptome_fasta = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_replicates.transcriptome.fa"),
        genome_fasta = config["genome"],
        all_regions = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),
        chromosomes_t = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "replicates_STAR_index",
            "chrNameLength.txt"),


    output:
        test = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_final",
                "peak_center",
                "test.fasta")
            ),
        test_bg = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_final",
                "peak_center",
                "test_bg.fasta")
            ),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "motif_analysis_final",
                "peak_center",
                "tmpdir")
            )
        ),
    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "total_enrichment_peakstofasta.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "motif_analysis_final",
            "peak_center",
            ),
        genome_tag = config["genome_tag"],
        peak_size = 50
        
    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_final",
            "peaks_split_sets__all_motifs_standard.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "motif_analysis_final",
            "peaks_split_sets__all_motifs_standard.stderr.log"
            ),

    shell:
        "(mkdir -p {params.output_dir}; \
        mkdir -p {output.dirtemp}; \
        python {params.script} \
        --peaks_file {input.reads} \
        --genome_fasta {input.genome_fasta} \
        --all_regions {input.all_regions} \
        --genome_tag {params.genome_tag} \
        --tempdir {output.dirtemp} \
        --transcriptome_fasta {input.transcriptome_fasta} \
        --chromosomes_g {input.chromosomes_g} \
        --chromosomes_t {input.chromosomes_t} \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        --crosslink_type 'peak_center' \
        ) 1> {log.stdout} 2> {log.stderr}"

