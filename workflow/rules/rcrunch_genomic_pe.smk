# ________________________________________________________________________________
# Genomic rcrunch for paired end CLIP data
# ________________________________________________________________________________
import os

rule GN_map_STAR_pe:
    """ map reads to genome"""
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
                format=config[wildcards.sample]['format']
                ),
        reads2 = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "cutadapt",
                    "{sample}_mate2.pe.cutadapt_{format}.fastq.gz"),
                sample=wildcards.sample,
                format=config[wildcards.sample]['format']
                ),
    
    output:
        bamfile = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}",
            "{sample}.pe.bam"
            ),
        logfile = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}",
            "{sample}.pe.Log.final.out"),
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
            "{sample}.pe.",
            ),
        multimappers = config['multimappers']
    
    singularity:
        "docker://zavolab/star:2.6.0a"

    shadow: "full"
    
    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "map_STAR.pe.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "map_STAR.pe.stderr.log"
            ),

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
        --outFilterScoreMinOverLread 0.2 \
        --outFilterMatchNminOverLread 0.2 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterType BySJout \
        --outReadsUnmapped Fastx \
        --outSAMattrRGline ID:rcrunch SM:{params.sample_id} \
        --alignEndsType EndToEnd > {output.bamfile} \
        ) 1> {log.stdout} 2> {log.stderr}"

rule GN_flag_duplicates_pe:
    """
        Duplicate removal in the absence of UMIs
        - detection step
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
            "flag_duplicates",
            "{sample}.duplicates.pe.Processed.out.bam"))

    params:
        cluster_log_path = config["cluster_log"],
        output_dir = os.path.join(
            config["output_dir"],
            "GN",
            "flag_duplicates"
            ),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "GN",
            "flag_duplicates",
            "{sample}.duplicates.pe."
            ),

    singularity:
        "docker://zavolab/star:2.6.0a"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "flag_duplicates.pe.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "flag_duplicates.pe.stderr.log"
            ),

    shell:
        "(STAR \
        --inputBAMfile {input.bam} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --runMode inputAlignmentsFromBAM \
        --outFileNamePrefix {params.outFileNamePrefix} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_remove_duplicates_pe:
    """
        Duplicate removal in the absence of UMIs
        - removal step
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "flag_duplicates",
            "{sample}.duplicates.pe.Processed.out.bam"
            ),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.duplicates.pe.bam"
            )),
    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_filter_duplicates.py"
            ),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "remove_duplicates.pe.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "remove_duplicates.pe.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --outfile {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_no_duplicate_removal_pe:
    """ 
        No duplicate removal
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
            "{sample}.with_duplicates.pe.bam"
            )),
        bai = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.with_duplicates.pe.bam.bai"
            )),
    params:
        cluster_log_path = config["cluster_log"],
    singularity:
        "docker://bash:5.0.16"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "no_duplicate_collapse.pe.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "no_duplicate_collapse.pe.stderr.log"
            ),

    shell:
        "(cp {input.bam} {output.bam}; \
          cp {input.bai} {output.bai}; \
          ) 1> {log.stdout} 2> {log.stderr}"


rule GN_umi_collapse_pe:
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
            "{sample}.umis.pe.bam"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        metrics_file = os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.umis.pe.metrics"
            ),

    singularity:
        "docker://zavolab/umi-tools:0.5.4"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "umi_collapse.pe.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "umi_collapse.pe.stderr.log"
            ),

    shell:
        "(umi_tools \
        dedup \
        -I {input.bam} \
        --paired \
        -S {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"
