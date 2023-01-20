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

rule TR_initial_map_to_genome_se:
    """
        Initial mapping to the genome using STAR.
    """
    input:
        genome = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),
        reads = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "cutadapt",
                    "{sample}_mate1.se.cutadapt_{format}.fastq.gz"),
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
                "{sample}.se.bam")
                ),
        logfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}",
                "{sample}.se.Log.final.out")
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
            "{sample}.se."),
        multimappers = config["multimappers"],
        mismatches = lambda wildcards: get_mismatches(wildcards.sample)
    
    shadow: "full"
    
    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_map_to_genome.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_map_to_genome.se.stderr.log"),

    shell:
        "(STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.genome} \
        --readFilesIn {input.reads} \
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
        --outFilterMismatchNoverLmax {params.mismatches} \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.sample_id} \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"

rule TR_flag_duplicates_se:
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
                "{sample}.duplicates.se.Processed.out.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.se."),

    singularity:
        "docker://zavolab/star:2.6.0a"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_duplicates.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_duplicates.se.stderr.log"),

    shell:
        "(STAR \
        --inputBAMfile {input.bam} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --runMode inputAlignmentsFromBAM \
        --outFileNamePrefix {params.outFileNamePrefix} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_remove_duplicates_se:
    """
        Remove reads flagged as duplicates by STAR.
        Custom script.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.se.Processed.out.bam"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.duplicates.se.bam")
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
            "{sample}_remove_duplicates.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_remove_duplicates.se.stderr.log"),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --outfile {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_no_duplicate_removal_se:
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
                "{sample}.with_duplicates.se.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.with_duplicates.se.bam.bai")
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
            "{sample}_no_duplicate_removal.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_no_duplicate_removal.se.stderr.log"),

    shell:
        "(cp {input.bam} {output.bam}; \
        cp {input.bai} {output.bai}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_umi_collapse_se:
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
            "{sample}.umis.se.bam")
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
            "{sample}_umi_collapse.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_umi_collapse.se.stderr.log"),

    shell:
        "(umi_tools \
        dedup \
        -I {input.bam} \
        -S {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_fastq_se:
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
                "{experiment}_{name}.mate1.se.unfiltered.fastq")
                )
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
            "{name}_bam_to_fastq.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq.se.stderr.log"),

    shell:
        "(bedtools bamtofastq \
        -i {input.bam} \
        -fq {output.fastq_mate1};) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_fastq_filtered_se:
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
            "{experiment}_{name}.mate1.se.unfiltered.fastq"),

    output:
        mate1 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate1.se.fastq.gz")
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
            "{name}_bam_to_fastq_filtered.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq_filtered.se.stderr.log"),

    shell:
        "(python {params.script} \
        --fastq_in {input.mate1} \
        --fastq_out {output.mate1};) 1> {log.stdout} 2> {log.stderr}"


rule TR_tr_quantification_se:
    """
        Quantification of transcript isoforms using Salmon.
    """
    input:
        transcriptome = os.path.join(
            config["output_dir"],
            "TR_salmon_index"),
        reads = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate1.se.fastq.gz"),

    output:
        quantification = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_quantification_se_{name}",
                "quant.sf")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        library = lambda wildcards: get_library_type(
            config[wildcards.experiment][wildcards.name][0],
            config[config[wildcards.experiment][wildcards.name][0]]['sense']),
        outdir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_quantification_se_{name}"),

    shadow: "full"
    singularity:
        "docker://zavolab/salmon:0.11.0"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_tr_quantification.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_tr_quantification.se.stderr.log"),

    shell:
        "(salmon quant \
        -i {input.transcriptome} \
        -l {params.library} \
        -r {input.reads} \
        -o {params.outdir} \
        -p {threads} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_map_to_transcriptome_se:
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
            "{experiment}_{name}.mate1.se.fastq.gz"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map",
                "{name}",
                "transcriptome.se.bam")
                ),
        logfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map",
                "{name}",
                "transcriptome.se.Log.final.out")
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
            "transcriptome.se."),
        multimappers = config['multimappers'] + 10,
        mismatches = lambda wildcards: get_mismatches(
            config[wildcards.experiment][wildcards.name][0]),
    
    shadow: "full"
    
    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_transcriptome.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_transcriptome.se.stderr.log"),

    shell:
        "(STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.transcriptome} \
        --readFilesIn {input.mate1} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 1 \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outFilterScoreMinOverLread 0.9 \
        --outFilterMatchNminOverLread 0.9 \
        --outFilterMismatchNoverLmax {params.mismatches} \
        --outSAMattributes All \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.experiment_id} \
        --alignIntronMax 1 \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"



rule TR_map_to_genome_se:
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
            "{experiment}_{name}.mate1.se.fastq.gz"),

    output:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome.se.bam"),
        logfile = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome.se.Log.final.out"),
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
            "genome.se."),
        multimappers = config["multimappers"],
        mismatches = lambda wildcards: get_mismatches(
            config[wildcards.experiment][wildcards.name][0]),
    
    shadow: "full"
    
    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_genome.se.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_genome.se.stderr.log"),

    shell:
        "(STAR \
        --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.genome} \
        --readFilesIn {input.mate1} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMinOverLread 0.9 \
        --outFilterMatchNminOverLread 0.9 \
        --outFilterMismatchNoverLmax {params.mismatches} \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outSAMattributes All \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.experiment_id} \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"