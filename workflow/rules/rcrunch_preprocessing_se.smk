# ________________________________________________________________________________
# Preprocessing subworkflow for paired-end data. First step of analysis
# ________________________________________________________________________________
import os
rule create_index:
    """
        Create an index of the human genome
    """
    input:
        genome = config["genome"],
        gtf = config["gtf"]

    output:
        chromosome_info = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"),

        chromosomes_names = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrName.txt")

    params:
        cluster_log_path = config["cluster_log"],
        output_dir = os.path.join(
            config["output_dir"],
            "STAR_index"),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "STAR_index/STAR_"),
        sjdbOverhang = config["sjdbOverhang"]

    singularity:
        "docker://zavolab/star:2.6.0a"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "STAR_index.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "STAR_index.stderr.log"
            ),

    shell:
        "(mkdir -p {params.output_dir}; \
        chmod -R 777 {params.output_dir}; \
        STAR \
        --runMode genomeGenerate \
        --sjdbOverhang {params.sjdbOverhang} \
        --genomeDir {params.output_dir} \
        --genomeFastaFiles {input.genome} \
        --runThreadN {threads} \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --sjdbGTFfile {input.gtf}) 1> {log.stdout} 2> {log.stderr}"


rule create_csv_from_ensembl_gtf:
    """
        Create a csv based on the gtf file of ensembl
    """
    input:
        gtf = config["gtf"]

    output:
        output = os.path.join(
            config["output_dir"],
            "ensembl_csv",
            "ensembl.csv"),

        flag = os.path.join(
            config["output_dir"],
            "ensembl_csv",
            "done"),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_ensembl__gtf_to_csv.py"
            ),

    singularity:
        "docker://zavolab/rcrunch_python:1.0"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "ensembl_csv.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "ensembl_csv.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --gtf_file {input.gtf} \
        --flag {output.flag} \
        --out_file {output.output}) 1> {log.stdout} 2> {log.stderr}"

rule cutadapt:
    """
        Trim 3' and 5' adapters
    """
    input:
        reads = lambda wildcards:
            os.path.join(
                config["input_dir"],
                config[wildcards.sample]['mate1'] + '.fastq.gz'),
        mate1_3p = lambda wildcards: config[wildcards.sample]['mate1_3p'],
        mate1_5p = lambda wildcards: config[wildcards.sample]['mate1_5p'],
    output:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt_replicates",
            "{sample}.cutadapt_standard.fastq.gz")

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/cutadapt:1.16"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "cutadapt__{sample}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "cutadapt__{sample}.stderr.log"
            )

    shell:
        "cutadapt \
        -e 0.1 \
        -O 1 \
        -j {threads} \
        -m 2 \
        -n 1 \
        -a file:{input.mate1_3p} \
        -g file:{input.mate1_5p} \
        -o {output.reads} \
        {input.reads}) 1> {log.stdout} 2> {log.stderr}"


rule umi_tools_format:
    """
        Create a UMI-tools barcode collapse format for any ENCODE samples.
    """
    input:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt_replicates2",
            "{sample}.cutadapt_standard.fastq.gz")

    output:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt_replicates2",
            "{sample}.cutadapt_encode.fastq.gz")
    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_convert_fastq_to_umi_format.py"
            ),
        read_format = lambda wildcards:
            config[wildcards.sample]['format']
        
    singularity:
        "docker://zavolab/rcrunch_python:1.0"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "umi_format__{sample}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "umi_format__{sample}.stderr.log"
            ),
    shell:
        "(python {params.script} \
        --infile {input.reads} \
        --outfile {output.reads} \
        --read_format {params.read_format};) 1> {log.stdout} 2> {log.stderr}"
