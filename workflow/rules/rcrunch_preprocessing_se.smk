# ________________________________________________________________________________
# Preprocessing subworkflow for paired-end data. First step of analysis
# ________________________________________________________________________________
import os

rule cutadapt_se:
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
            "cutadapt",
            "{sample}_mate1.se.cutadapt_standard.fastq.gz")

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/cutadapt:1.16"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "cutadapt__{sample}.se.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "cutadapt__{sample}.se.stderr.log"
            )

    shell:
        "(cutadapt \
        -e 0.1 \
        -O 1 \
        -j {threads} \
        -m 15 \
        -n 1 \
        -a file:{input.mate1_3p} \
        -g file:{input.mate1_5p} \
        -o {output.reads} \
        {input.reads}) 1> {log.stdout} 2> {log.stderr}"


rule extract_umi_se:
    '''
        Extract UMIs and add them in the sequence header 

    '''
    input:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt",
            "{sample}_mate1.se.cutadapt_standard.fastq.gz")
                
    output:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt",
            "{sample}_mate1.se.cutadapt_umi.fastq.gz")

    params:
        cluster_log_path = config["cluster_log"],
        log = os.path.join(
            config["output_dir"],
            "cutadapt",
            "{sample}.extractumis.log"),
        umi = lambda wildcards: get_umi(wildcards.sample),

    singularity:
        "docker://quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "extract_umi__{sample}.se.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "extract_umi__{sample}.se.stderr.log"
            )

    shell:
        "(umi_tools extract \
        --extract-method=string \
        --stdin={input.reads} \
        --bc-pattern={params.umi} \
        --log={params.log} \
        --stdout {output.reads}) 1> {log.stdout} 2> {log.stderr}"



rule umi_tools_format_se:
    """
        Create a UMI-tools barcode collapse format for any ENCODE samples.
    """
    input:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt",
            "{sample}_mate1.se.cutadapt_standard.fastq.gz")

    output:
        reads = os.path.join(
            config["output_dir"],
            "cutadapt",
            "{sample}_mate1.se.cutadapt_encode.fastq.gz")
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
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "preprocessing",
            "umi_format__{sample}.se.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "preprocessing",
            "umi_format__{sample}.se.stderr.log"
            ),
    shell:
        "(python {params.script} \
        --infile {input.reads} \
        --outfile {output.reads} \
        --read_format {params.read_format};) 1> {log.stdout} 2> {log.stderr}"
