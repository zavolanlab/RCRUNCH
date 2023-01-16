# ________________________________________________________________________________
# Genomic rcrunch for paired end CLIP data
# ________________________________________________________________________________
import os


rule GN_index_mappings:
    """
        Sort and index alignments
        according to coordinates
    """
    input:  
        bam = lambda wildcards: expand(
            os.path.join(
                config["output_dir"],
                "GN",
                "alignment",
                "{sample}",
                "{sample}.{mates}.bam"),
            sample=wildcards.sample,
            mates=get_mates(wildcards.sample)),

        logfile = lambda wildcards: expand(
            os.path.join(
                config["output_dir"],
                "GN",
                "alignment",
                "{sample}",
                "{sample}.{mates}.Log.final.out"),
            sample=wildcards.sample,
            mates=get_mates(wildcards.sample)),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}_sorted.bam"
            )),

        bai = temp(os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}_sorted.bam.bai"
            )),
    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}_temp"
            ),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "sort_mappings.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "sort_mappings.stderr.log"
            ),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_flag_ncRNA_reads:
    """
        Remove reads mapping to specific
        nc categories specified in the config
        - identification step
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}_sorted.bam"
            ),
        bai = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}_sorted.bam.bai"
            ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_filter_ncRNAs.py"
            ),
        ncRNAs = config['ncRNAs'],
        ncRNA_biotypes = expand(config['ncRNA_biotypes']),
        paired = lambda wildcards:
            get_mates_number(wildcards.sample),
        sense = lambda wildcards:
            config[wildcards.sample]['sense'],
        
    output:
        outfile = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.ncrna.txt")),
        flag = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}_done"
            )),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "flag_ncRNA_reads.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "flag_ncRNA_reads.stderr.log"
            ),

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


rule GN_remove_ncRNA_reads:
    """
        Removal of the ncreads 
        - removal step
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "alignment",
            "{sample}_sorted.bam"
            ),
        read_names = os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.ncrna.txt"
            ),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.filtered.bam"
            )),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/picard:2.18.9"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "remove_ncRNA_reads.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "remove_ncRNA_reads.stderr.log"
            ),

    shell:
        "(java -jar /usr/local/bin/picard.jar FilterSamReads \
        I={input.bam} \
        O={output.bam} \
        READ_LIST_FILE={input.read_names} \
        FILTER=excludeReadList \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_sort_ncRNA_rm:
   """
        Sort and index alignments
        according to coordinates
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.filtered.bam"
            ),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.filtered.sorted.bam"
            )),
        bai = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}.filtered.sorted.bam.bai"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "GN",
            "remove_ncRNAs",
            "{sample}_temp"
            ),

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "sort_ncrna_rm.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "sort_ncrna_rm.stderr.log"
            ),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_flag_duplicates:
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
            "{sample}.duplicates.Processed.out.bam"))

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
            "{sample}.duplicates."
            ),

    singularity:
        "docker://zavolab/star:2.6.0a"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "flag_duplicates.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "flag_duplicates.stderr.log"
            ),

    shell:
        "(STAR \
        --inputBAMfile {input.bam} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --runMode inputAlignmentsFromBAM \
        --outFileNamePrefix {params.outFileNamePrefix} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_remove_duplicates:
    """
        Duplicate removal in the absence of UMIs
        - removal step
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "flag_duplicates",
            "{sample}.duplicates.Processed.out.bam"
            ),

    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.duplicates.bam"
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
            "remove_duplicates.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "remove_duplicates.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --outfile {output.bam} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_no_duplicate_removal:
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
            "{sample}.with_duplicates.bam"
            )),
        bai = temp(os.path.join(
            config["output_dir"],
            "GN",
            "remove_duplicates",
            "{sample}.with_duplicates.bam.bai"
            )),

    singularity:
        "docker://bash:5.0.16"

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "no_duplicate_collapse.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{sample}",
            "no_duplicate_collapse.stderr.log"
            ),

    shell:
        "(cp {input.bam} {output.bam}; \
          cp {input.bai} {output.bai}; \
          ) 1> {log.stdout} 2> {log.stderr}"


rule GN_merge:
    input:
        bam = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "GN",
                    "remove_duplicates",
                    "{sample}.{dup_type}.{mates}.bam"),
                sample=config[wildcards.experiment][wildcards.name],
                dup_type=config[config[wildcards.experiment][wildcards.name][0]]['dup_type'],
                mates=get_mates(config[wildcards.experiment][wildcards.name][0]),
                ),
    output:
        bam = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.bam"
            )),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "merge_samples__{name}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "merge_samples__{name}.stderr.log"
            ),

    shell:
        "(samtools merge \
        {output.bam} {input.bam}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_sort_merged:
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.bam"
            ),

    output:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam"
            ),
        bai = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam.bai"
            ),

    params:
        cluster_log_path = config["cluster_log"],
        prefix = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples_temp")

    singularity:
        "docker://zavolab/samtools:1.8"

    threads: 8

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "merge_samples__{name}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "merge_samples__{name}.stderr.log"
            ),

    shell:
        "(samtools sort \
        -T {params.prefix} \
        -@ {threads} \
        {input.bam} > {output.bam}; \
        samtools index \
        {output.bam} {output.bai} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_read_frequencies:
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam"
            ),

    output:
        frequencies = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}.frequencies.csv"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        paired = lambda wildcards:
            get_mates_number(config[wildcards.experiment][wildcards.name][0]),
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_bam_get_read_frequencies.py"
            ),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "read_frequencies__{name}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "read_frequencies__{name}.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --bamfile {input.bam} \
        --paired {params.paired} \
        --outfile {output.frequencies}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_bam_to_bed:
    """
        Convert bamfile to bedfile using bedtools
    """
    input:
        bam = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.sorted.bam"
            ),

    output:
        bed = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_{name}_merge_samples.bed"
            )),

    params:
        cluster_log_path = config["cluster_log"],

    singularity:
        "docker://zavolab/bedtools:2.27.0"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "bam_to_bed___{name}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "bam_to_bed__{name}.stderr.log"
            ),

    shell:
        "(bedtools bamtobed \
        -i {input.bam} > {output.bed}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_genome_coverage:
    """
        Obtain windows of read coverage genome
    """
    input:
        chromosome_info = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),
        fg_bed = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_replicates_merge_samples.bed"
            ),
        bg_bed = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_smis_merge_samples.bed"
            ),
        fg_frequencies = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_replicates.frequencies.csv"
            ),
        bg_frequencies = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_smis.frequencies.csv"
            ),

    output:
        coverage = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}.gn_coverage.bed"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "sliding_windows.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "GN",
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
            get_mates_number(config[wildcards.experiment]["replicates"][0]),
        paired_b = lambda wildcards:
            get_mates_number(config[wildcards.experiment]["smis"][0]),
        prefix = "{experiment}.gn_coverage",
        background_type = lambda wildcards:
            config[wildcards.experiment]["background_type"],

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "genome_coverage.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "genome_coverage.stderr.log"
            ),

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


rule GN_enriched_regions:
    """
        Obtain enriched regions for binding based on EM using
        bg for noise estimation. Modified and adapted from:
        Berger, Severin M., et al. "Crunch: Integrated processing and modeling of
        ChIP-seq data in terms of regulatory motifs."
        Genome research (2019): gr-239319
        Custom script.
   """
    input:
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}.gn_coverage.bed"
            ),

    output:
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_enriched_regions.csv"
            ),
        parameters = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_parameters.csv"
            ),
        all_regions = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),

    params:
        cluster_log_path = config["cluster_log"],
        plot_dir = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "enriched_regions_plots"),
        fdr_cutoff = config['FDR'],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_find_enriched_regions.py"
            ),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "enriched_regions.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "enriched_regions.stderr.log"
            ),

    shell:
        "(mkdir -p {params.plot_dir}; \
        python {params.script} \
        --windows_table {input.reads} \
        --out {output.reads} \
        --out_parameters {output.parameters} \
        --all_regions {output.all_regions} \
        --out_folder {params.plot_dir} \
        --fdr_cutoff {params.fdr_cutoff} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_fit_peaks:
    """
        Obtain peaks
    """
    input:
        fg_bam = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_replicates_merge_samples.sorted.bam"
            ),
        bg_bam = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_smis_merge_samples.sorted.bam"
            ),
        fg_frequencies = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_replicates.frequencies.csv"
            ),
        bg_frequencies = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_smis.frequencies.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_enriched_regions.csv"
            ),
        parameters = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_parameters.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),

    output:
        peaks = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_define_peaks.py"
            ),
        plot_dir = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "fit_peaks_plots"),
        fragment_size = config["fragment_size"],
        window_size = lambda wildcards:
            config[wildcards.experiment]['window_f'],
        paired_f = lambda wildcards:
            get_mates_number(config[wildcards.experiment]["replicates"][0]),
        paired_b = lambda wildcards:
            get_mates_number(config[wildcards.experiment]["smis"][0]),
        sense_f = lambda wildcards:
            config[config[wildcards.experiment]["replicates"][0]]['sense'],
        sense_b = lambda wildcards:
            config[config[wildcards.experiment]["smis"][0]]['sense'],
        background_type = lambda wildcards:
            config[wildcards.experiment]["background_type"],        

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 12

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "fit_peaks.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "fit_peaks.stderr.log"
            ),

    shell:
        "(mkdir -p {params.plot_dir}; \
        python {params.script} \
        --enriched_regions {input.reads} \
        --parameters {input.parameters} \
        --fragment_size {params.fragment_size} \
        --window_size {params.window_size} \
        --chromosomes_g {input.chromosomes_g} \
        --bam_foreground_g {input.fg_bam} \
        --bam_background_g {input.bg_bam} \
        --bam_fg_fq {input.fg_frequencies} \
        --bam_bg_fq {input.bg_frequencies} \
        --paired_f {params.paired_f} \
        --paired_b {params.paired_b} \
        --sense_f {params.sense_f} \
        --sense_b {params.sense_b} \
        --background {params.background_type} \
        --out {output.peaks} \
        --out_folder {params.plot_dir} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_peaks_split_sets_crosslink:
    """
        Split peaks into training and test set
    """
    input:
        genome_fasta = config["genome"],
        all_regions = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),

    output:
        training = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "crosslink",
            "training.fasta"
            )),
        test = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "crosslink",
            "test.fasta"
            )),
        training_bg = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "crosslink",
            "training_bg.fasta"
            )),
        test_bg = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "crosslink",
            "test_bg.fasta"
            )),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "GN",
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
            "mk_peaks_split_sets.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "crosslink"
            ),
        genome_tag = config["genome_tag"],
        peak_size = 50,
        
    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "peaks_split_sets__{run}_crosslink.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "peaks_split_sets__{run}_crosslink.stderr.log"
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
        --chromosomes_g {input.chromosomes_g} \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        --crosslink_type 'crosslink' \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_peaks_split_sets:
    """
        Split peaks into training and test set
    """
    input:
        genome_fasta = config["genome"],
        all_regions = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),

    output:
        training = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "peak_center",
            "training.fasta"
            )),
        test = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "peak_center",
            "test.fasta"
            )),
        training_bg = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "peak_center",
            "training_bg.fasta"
            )),
        test_bg = temp(os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "peak_center",
            "test_bg.fasta"
            )),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "GN",
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
            "mk_peaks_split_sets.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "peak_center"
            ),
        genome_tag = config["genome_tag"],
        peak_size = 50
        
    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "peaks_split_sets__{run}_peak_center.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "peaks_split_sets__{run}_peak_center.stderr.log"
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
        --chromosomes_g {input.chromosomes_g} \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        --crosslink_type 'peak_center' \
        ) 1> {log.stdout} 2> {log.stderr}"

rule GN_peaks_split_sets_crosslink_all_motifs:
    """
        Get all peaks crosslink
    """
    input:
        genome_fasta = config["genome"],
        all_regions = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),

    output:
        test = temp(
            os.path.join(
                config["output_dir"],
                "GN",
                "{experiment}",
                "motif_analysis_final",
                "crosslink",
                "test.fasta")
            ),
        test_bg = temp(
            os.path.join(
                config["output_dir"],
                "GN",
                "{experiment}",
                "motif_analysis_final",
                "crosslink",
                "test_bg.fasta")
            ),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "GN",
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
            "GN",
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
            "GN",
            "{experiment}",
            "motif_analysis_final",
            "peaks_split_sets__all_motifs_crosslink.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
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
        --chromosomes_g {input.chromosomes_g} \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        --crosslink_type 'crosslink' \
        ) 1> {log.stdout} 2> {log.stderr}"


rule GN_peaks_split_sets_all_motifs:
    """
        get all peaks crosslink
    """
    input:
        genome_fasta = config["genome"],
        all_regions = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_regions.csv"
            ),
        reads = os.path.join(
            config["output_dir"],
            "GN",
            "{experiment}",
            "{experiment}_total_peaks.csv"
            ),
        chromosomes_g = os.path.join(
            config["output_dir"],
            "STAR_index",
            "chrNameLength.txt"
            ),

    output:
        test = temp(
            os.path.join(
                config["output_dir"],
                "GN",
                "{experiment}",
                "motif_analysis_final",
                "peak_center",
                "test.fasta")
            ),
        test_bg = temp(
            os.path.join(
                config["output_dir"],
                "GN",
                "{experiment}",
                "motif_analysis_final",
                "peak_center",
                "test_bg.fasta")
            ),
        dirtemp = temp(
            directory(
                os.path.join(
                config["output_dir"],
                "GN",
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
            "GN",
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
            "GN",
            "{experiment}",
            "motif_analysis_final",
            "peaks_split_sets__all_motifs_standard.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "GN",
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
        --chromosomes_g {input.chromosomes_g} \
        --peak_size {params.peak_size} \
        --out_folder {params.output_dir} \
        --crosslink_type 'peak_center' \
        ) 1> {log.stdout} 2> {log.stderr}"

