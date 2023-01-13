import os

rule GN_map_STAR_se:
    """ map reads to genome"""
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
        bamfile = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}",
            "{sample}.se.bam"
            ),
        logfile = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}",
            "{sample}.se.Log.final.out"),
    params:
        cluster_log_path = config["cluster_log"],
        sample_id = "{sample}",
        genome = os.path.join(
            config["output_dir"],
            "STAR_index"
            ),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}",
            "{sample}.se.",
            ),
        multimappers = config['multimappers']
    shadow: "full"

    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "map_STAR.se.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "map_STAR.se.stderr.log"
            ),

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
        --outFilterScoreMinOverLread 0.2 \
        --outFilterMatchNminOverLread 0.2 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterType BySJout \
        --outReadsUnmapped Fastx \
        --outSAMattrRGline ID:rcrunch SM:{params.sample_id} \
        --alignEndsType EndToEnd > {output.bamfile} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_umi_collapse_se:
    """
        Duplicate removal when there are UMIs
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.filtered.sorted.bam"
            ),
        bai = os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.filtered.sorted.bam.bai"
            ),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.umis.se.bam"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        metrics_file = os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.umis.se.metrics"
            ),

    singularity:
        "docker://zavolab/umi-tools:0.5.4"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "umi_collapse.se.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "umi_collapse.se.stderr.log"
            ),

    shell:
        "(umi_tools \
        dedup \
        -I {input.bam} \
        -S {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"
