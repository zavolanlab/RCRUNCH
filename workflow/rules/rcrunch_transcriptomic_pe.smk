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

rule TR_initial_map_to_genome_pe:
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
                    "{sample}_mate1.pe.cutadapt_{format}.fastq.gz"),
                sample=wildcards.sample,
                format=config[wildcards.sample]['format'],
                ),
        reads2 = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "cutadapt",
                    "{sample}_mate2.pe.cutadapt_{format}.fastq.gz"),
                sample=wildcards.sample,
                format=config[wildcards.sample]['format'],
                ),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}",
                "{sample}.pe.bam")),
        logfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "initial_map",
                "{sample}",
                "{sample}.pe.Log.final.out")),

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
            "{sample}.pe."),
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
            "{sample}_initial_map_to_genome.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_initial_map_to_genome.pe.stderr.log"),

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
        --outFilterMismatchNoverLmax {params.mismatches} \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rcrunch SM:{params.sample_id} \
        --alignEndsType EndToEnd > {output.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"

rule TR_flag_duplicates_pe:
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
                "{sample}.duplicates.pe.Processed.out.bam")
                ),

    params:
        cluster_log_path = config["cluster_log"],
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.pe."),

    singularity:
        "docker://zavolab/star:2.6.0a"

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_duplicates.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_flag_duplicates.pe.stderr.log"),

    shell:
        "(STAR \
        --inputBAMfile {input.bam} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --runMode inputAlignmentsFromBAM \
        --outFileNamePrefix {params.outFileNamePrefix} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_index_dupl_mappings_pe:
    """
        Sort and index alignments
        according to coordinates
    """
    input:  
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.pe.Processed.out.bam")

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.pe.bam"
            )),

        bai = temp(os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.pe.bam.bai"
            )),
    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.pe_temp"
            ),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{sample}",
            "sort_dupl_mappings.pe.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{sample}",
            "sort_dupl_mappings.pe.stderr.log"
            ),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_remove_duplicates_pe:
    """
        Remove reads flagged as duplicates by STAR.
        Custom script.
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "flag_duplicates",
            "{sample}.duplicates.pe.bam"),
    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.duplicates.pe.bam")
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
            "{sample}_remove_duplicates.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_remove_duplicates.pe.stderr.log"),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --outfile {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_no_duplicate_removal_pe:
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
                "{sample}.with_duplicates.pe.bam")
                ),
        bai = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "remove_duplicates",
                "{sample}.with_duplicates.pe.bam.bai")
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
            "{sample}_no_duplicate_removal.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_no_duplicate_removal.pe.stderr.log"),

    shell:
        "(cp {input.bam} {output.bam}; \
        cp {input.bai} {output.bai}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_umi_collapse_pe:
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
            "{sample}.umis.pe.bam")
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
            "{sample}_umi_collapse.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "preprocessing",
            "{sample}_umi_collapse.pe.stderr.log"),

    shell:
        "(umi_tools \
        dedup \
        -I {input.bam} \
        --paired \
        -S {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_fastq_pe:
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
                "{experiment}_{name}.mate1.pe.unfiltered.fastq")
                ),
        fastq_mate2 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate2.pe.unfiltered.fastq")
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
            "{name}_bam_to_fastq.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq.pe.stderr.log"),

    shell:
        "(bedtools bamtofastq \
        -i {input.bam} \
        -fq {output.fastq_mate1} \
        -fq2 {output.fastq_mate2}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_bam_to_fastq_filtered_pe:
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
            "{experiment}_{name}.mate1.pe.unfiltered.fastq"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.pe.unfiltered.fastq"),

    output:
        mate1 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate1.pe.fastq.gz")
                ),
        mate2 = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "{experiment}_{name}.mate2.pe.fastq.gz")
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
            "{name}_bam_to_fastq_filtered.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_bam_to_fastq_filtered.pe.stderr.log"),

    shell:
        "(python {params.script} \
        --fastq_in {input.mate1} \
        --fastq_out {output.mate1} \
        --fastq_in2 {input.mate2} \
        --fastq_out2 {output.mate2}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_tr_quantification_pe:
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
            "{experiment}_{name}.mate1.pe.fastq.gz"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.pe.fastq.gz"),

    output:
        quantification = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_quantification_pe_{name}",
                "quant.sf")
                ),
    params:
        cluster_log_path = config["cluster_log"],
        library = lambda wildcards: get_library_type(
            config[config[wildcards.experiment][wildcards.name][0]]['mate2'],
            config[config[wildcards.experiment][wildcards.name][0]]['sense']),
        outdir = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "tr_quantification_pe_{name}"),

    shadow: "full"

    singularity:
        "docker://zavolab/salmon:0.11.0"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_tr_quantification.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_tr_quantification.pe.stderr.log"),

    shell:
        "(salmon quant \
        -i {input.transcriptome} \
        -l {params.library} \
        -1 {input.mate1} \
        -2 {input.mate2} \
        -o {params.outdir} \
        -p {threads} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule TR_map_to_transcriptome_pe:
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
            "{experiment}_{name}.mate1.pe.fastq.gz"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.pe.fastq.gz"),

    output:
        bam = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map",
                "{name}",
                "transcriptome.pe.bam")
                ),
        logfile = temp(
            os.path.join(
                config["output_dir"],
                "TR",
                "{experiment}",
                "tr_map",
                "{name}",
                "transcriptome.pe.Log.final.out")
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
            "transcriptome.pe."),
        multimappers = config['multimappers'] + 10,
        mismatches = lambda wildcards: get_mismatches(
            config[wildcards.experiment][wildcards.name][0])

        
    shadow: "full"
    
    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_transcriptome.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_transcriptome.pe.stderr.log"),

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



rule TR_map_to_genome_pe:
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
            "{experiment}_{name}.mate1.pe.fastq.gz"),
        mate2 = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "{experiment}_{name}.mate2.pe.fastq.gz"),

    output:
        bam = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome.pe.bam"),
        logfile = os.path.join(
            config["output_dir"],
            "TR",
            "{experiment}",
            "gn_map",
            "{name}",
            "genome.pe.Log.final.out"),
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
            "genome.pe."),
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
            "{name}_map_to_genome.pe.stdout.log"),
        stderr = os.path.join(
            config["local_log"],
            "TR",
            "{experiment}",
            "{name}_map_to_genome.pe.stderr.log"),

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