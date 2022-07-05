import os

rule motif_analysis__randomise_sequences:
    """
        Custom: randomise sequences in the test
        and training datasets simple sequences.
    """
    input:
        sequences = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "{set}.fasta"
            ),
        sequences_bg = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "{set}_bg.fasta"
            ),
        
    output:
        reads = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}",
            "{set}.fasta"
            )),
        randomised_bg = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}",
            "{set}_bg.fasta"
            )),
        reads_pool = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}",
            "{set}_pool.fasta"
            )),
        reads_motevo = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}",
            "{set}_motevo.fasta"
            )),
        randomised_bg_motevo = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}",
            "{set}_bg_motevo.fasta"
            )),
        reads_pool_motevo = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}",
            "{set}_pool_motevo.fasta"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_randomise_sequences.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "{set}"
            ),
        random_seq_per_peak = config["random_sequences_per_peak"],
        prefix = "{set}",
        species = config['genome_tag']

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "randomise_sequences__{run}__{crosslink_type}__{set}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "randomise_sequences__{run}__{crosslink_type}__{set}.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --peaks {input.sequences} \
        --out_folder {params.output_dir} \
        --prefix {params.prefix} \
        --species {params.species} \
        --peaks_bg {input.sequences_bg} \
        --random_seq_per_peak_num {params.random_seq_per_peak} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule motif_analysis__PHYLOGIBBS:
    """
        Phylogibbs: find motif based based on the training set
    """
    input:
        peaks = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "training",
            "training.fasta"
            ),
    output:
        motif = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "PHYLOGIBBS_de_novo_{motif_length}"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        motif_length = "{motif_length}"

    singularity:
        "docker://zavolab/phylogibbs:1.2"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "phylogibbs__{run}__{motif_length}_{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "phylogibbs__{run}__{motif_length}_{crosslink_type}.stderr.log"
            ),

    shell:
        "(phylogibbs \
        -r \
        -D 0 \
        -m {params.motif_length} \
        -z 2 \
        -y $(( $(cat {input.peaks} | wc -l) / 4)) \
        -N 1 \
        -q \
        -f {input.peaks} \
        -o {output.motif}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule motif_analysis__phylogibbs_to_transfac:
    """
        Custom: convert the output of phylogibbs into
        the desired TRANSFAC format, which can be 
        used by Motevo.
    """
    input:
        motifs = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "PHYLOGIBBS_de_novo_{motif_length}"
            ),

    output:
        motif1 = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "de_novo_{motif_length}_1_{experiment}_{method_type}.transfac"
            )),
        motif2 = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "de_novo_{motif_length}_2_{experiment}_{method_type}.transfac"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_motif_phylogibbs_to_transfac.py"
            ),
        prefix = "de_novo_{motif_length}"

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "transfac_{run}__{motif_length}_{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "transfac_{run}__{motif_length}_{crosslink_type}.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --phylogibbs_motifs {input.motifs} \
        --out1 {output.motif1} \
        --out2 {output.motif2} \
        --prefix {params.prefix}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule motif_analysis__motevo__motif_trim:
    """
        Custom: trim customised motifs
    """
    input:
        motif = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "de_novo_{motif_length}_{motif_number}_{experiment}_{method_type}.transfac"
            ),

    output:
        trimmed_motif = temp(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "de_novo_{motif_length}_{motif_number}_{experiment}_{method_type}_{run}.trimmed"
            )),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_motevo_trim_motif.py"
            ),

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "motevo__motif_trim__{run}__" +
            "{motif_length}_{motif_number}_{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "motevo__motif_trim__{run}__" +
            "{motif_length}_{motif_number}_{crosslink_type}.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --inmotif {input.motif} \
        --outmotif {output.trimmed_motif}; \
        ) 1> {log.stdout} 2> {log.stderr}"


checkpoint motif_analysis__motevo__format_wlib:
    """
        Custom: motifs provided by Attract 
        -- format into separate transfac files
    """
    input:
        de_novo_motif = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "{method_type}",
                    "{experiment}",
                    "motif_analysis_crossvalidation",
                    "{run}",
                    "{crosslink_type}",
                    "de_novo_{motif_length}_{motif_no}_{experiment}_{method_type}_{run}.trimmed"
                    ),
                method_type=wildcards.method_type,
                experiment=wildcards.experiment,
                run=wildcards.run,
                crosslink_type=wildcards.crosslink_type,
                motif_no = ["1", "2"],
                motif_length=config["motif_lengths"])

    output:
        outdir = directory(os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "wlib")),

    params:
        motifs = config["wlib_pwms"],
        names = config["wlib_names"],
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_known_motifs_format.py"
            ),
        organism = config['wlib_organism'],
        rbp = lambda wildcards:
            config[wildcards.experiment]['rbp']

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "motevo__format__wlib_{run}__{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "motevo__format__wlib_{run}__{crosslink_type}.stderr.log"
            ),
    shell:
        "(mkdir -p {output.outdir}; \
        python {params.script} \
        --pwms {params.motifs} \
        --de_novo {input.de_novo_motif} \
        --names {params.names} \
        --organism {params.organism} \
        --rbp_name {params.rbp} \
        --outdir {output.outdir}; \
        ) 1> {log.stdout} 2> {log.stderr}"


def gather_wlib_motifs(wildcards):
    checkpoint_output = checkpoints.motif_analysis__motevo__format_wlib.get(
        **wildcards).output.outdir
    ivals = glob_wildcards(os.path.join(
        checkpoint_output, "motif_{wlib}.pwm")).wlib
    pathout = os.path.join(
        config["output_dir"],
        wildcards.method_type,
        wildcards.experiment,
        "motif_analysis_crossvalidation",
        wildcards.run,
        wildcards.crosslink_type,
        "wlib")
    return expand(
        os.path.join(
            pathout,
            "motif_{wlib}.pwm"),
        wlib=ivals)


rule motif_analysis__motevo:
    input:
        training_pool = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "training",
            "training_pool_motevo.fasta"
            ),
        training = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "training",
            "training_motevo.fasta"
            ),
        test_pool = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "randomise_sequences",
            "test",
            "test_pool_motevo.fasta"
            ),
        motifs = gather_wlib_motifs,

    output:
        enrichment = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            "motif_enrichment.tsv"
            ),
        tmpdir = temp(
            directory(
                os.path.join(
                    config["output_dir"],
                    "{method_type}",
                    "{experiment}",
                    "motif_analysis_crossvalidation",
                    "{run}",
                    "{crosslink_type}",
                    "motevo"
                    )
                )
        ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_motevo_tfbs.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "{run}",
            "{crosslink_type}",
            ),
        organism = config["genome_tag"]

    singularity:
        "docker://zavolab/motevo:1.12_python3.6.9"

    threads: 4

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "phylogenetic__motevo__tfbs_{run}__{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_crossvalidation",
            "phylogenetic__motevo__tfbs_{run}__{crosslink_type}.stderr.log"
            ),
    shell:
        "(mkdir -p {params.output_dir}; \
        python {params.script} \
        --outpath {params.output_dir} \
        --training {input.training} \
        --training_pool {input.training_pool} \
        --test_pool {input.test_pool} \
        --wms {input.motifs} \
        --genome_tag {params.organism} \
        --outfile {output.enrichment} \
        ) 1> {log.stdout} 2> {log.stderr}"


rule motif_analysis__randomise_sequences_all_motifs:
    """
        Custom: randomise sequences and make motevo compatible
    """
    input:
        sequences = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "{crosslink_type}",
            "test.fasta"
            ),
        sequences_bg = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "{crosslink_type}",
            "test_bg.fasta"
            ),
        
    output:
        reads = temp(
            os.path.join(
                config["output_dir"],
                "{method_type}",
                "{experiment}",
                "motif_analysis_final",
                "{crosslink_type}",
                "test",
                "test.fasta")
            ),
        randomised_bg = temp(
            os.path.join(
                config["output_dir"],
                "{method_type}",
                "{experiment}",
                "motif_analysis_final",
                "{crosslink_type}",
                "test",
                "test_bg.fasta")
            ),
        reads_pool = temp(
            os.path.join(
                config["output_dir"],
                "{method_type}",
                "{experiment}",
                "motif_analysis_final",
                "{crosslink_type}",
                "test",
                "test_pool.fasta")
            ),
        reads_motevo = temp(
            os.path.join(
                config["output_dir"],
                "{method_type}",
                "{experiment}",
                "motif_analysis_final",
                "{crosslink_type}",
                "test",
                "test_motevo.fasta")
            ),
        randomised_bg_motevo = temp(
            os.path.join(
                config["output_dir"],
                "{method_type}",
                "{experiment}",
                "motif_analysis_final",
                "{crosslink_type}",
                "test",
                "test_bg_motevo.fasta")
            ),
        reads_pool_motevo = temp(
            os.path.join(
                config["output_dir"],
                "{method_type}",
                "{experiment}",
                "motif_analysis_final",
                "{crosslink_type}",
                "test",
                "test_pool_motevo.fasta")
            ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "mk_randomise_sequences.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "{crosslink_type}",
            "test",
            ),
        random_seq_per_peak = config["random_sequences_per_peak"],
        prefix = "test",
        species = config['genome_tag']

    singularity:
        "docker://zavolab/rcrunch_python:1.0.5"

    threads: 1

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "randomise_sequences__test_all_motifs__{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "randomise_sequences__test_all_motifs_{crosslink_type}.stderr.log"
            ),

    shell:
        "(python {params.script} \
        --peaks {input.sequences} \
        --out_folder {params.output_dir} \
        --prefix {params.prefix} \
        --species {params.species} \
        --peaks_bg {input.sequences_bg} \
        --random_seq_per_peak_num {params.random_seq_per_peak} \
        ) 1> {log.stdout} 2> {log.stderr}"


def gather_wlib_motifs_all_motifs(wildcards):
    all_motifs = []
    for run in [f'run{i}' for i in np.arange(int(config['runs']))]:
        checkpoint_output = checkpoints.motif_analysis__motevo__format_wlib.get(
            **wildcards, run=run).output.outdir
        ivals = glob_wildcards(os.path.join(
            checkpoint_output, "motif_{wlib}.pwm")).wlib
        all_motifs.extend(
            expand(
                os.path.join(
                    config["output_dir"],
                    wildcards.method_type,
                    wildcards.experiment,
                    "motif_analysis_crossvalidation",
                    "{run}",
                    wildcards.crosslink_type,
                    "wlib",
                    "motif_{wlib}.pwm"),
                wlib=ivals,
                run=run)
            )
    return all_motifs


rule motif_analysis__motevo_all_motifs:
    input:
        motifs = gather_wlib_motifs_all_motifs,
        test_pool = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "{crosslink_type}",
            "test",
            "test_pool_motevo.fasta"
            )

    output:
        enrichment = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "{crosslink_type}",
            "motif_enrichment.tsv"
            ),

    params:
        cluster_log_path = config["cluster_log"],
        script = os.path.join(
            workflow.basedir,
            "scripts",
            "total_enrichment_motevo_tfbs.py"
            ),
        output_dir = os.path.join(
            config["output_dir"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "{crosslink_type}",
            "motevo"
            ),
        organism = config["genome_tag"]

    singularity:
        "docker://zavolab/motevo:1.12_python3.6.9"

    threads: 6

    log:
        stdout = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "phylogenetic__motevo__tfbs_all_motifs_{crosslink_type}.stdout.log"
            ),
        stderr = os.path.join(
            config["local_log"],
            "{method_type}",
            "{experiment}",
            "motif_analysis_final",
            "phylogenetic__motevo__tfbs_all_motifs_{crosslink_type}.stderr.log"
            ),
    shell:
        "(mkdir -p {params.output_dir}; \
        python {params.script} \
        --outpath {params.output_dir} \
        --test_pool {input.test_pool} \
        --wms {input.motifs} \
        --genome_tag {params.organism} \
        --outfile {output.enrichment} \
        ) 1> {log.stdout} 2> {log.stderr}"

# runs = expand(
#             os.path.join(
#                 config["output_dir"],
#                 "{method_type}",
#                 "{experiment}",
#                 "motif_analysis_crossvalidation",
#                 "{run}",
#                 "{crosslink_type}",
#                 "motif_enrichment.tsv"),
#             method_type=config["method_types"],
#             experiment=config["EXPERIMENT_SET"],
#             run=[f'run{i}' for i in np.arange(int(config['runs']))],
#             crosslink_type=config["peak_center"]),