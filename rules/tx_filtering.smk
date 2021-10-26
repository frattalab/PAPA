rule gtf_list_by_condition:
    input:
        lambda wildcards: expand(os.path.join(STRINGTIE_SUBDIR,
                                              "min_jnc_" + wildcards.min_jnc,
                                              "min_frac_" + wildcards.min_frac,
                                              "min_cov_" + wildcards.min_cov,
                                              "{sample}.assembled.gtf"
                                              ),
                                 sample=get_condition_samples(wildcards.condition, OPTIONS))

    output:
        txt = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}",
                           "min_frac_{min_frac}",
                           "min_cov_{min_cov}",
#                           "{condition}",
                           "gtf_list_{condition}.txt")

    group:
        "transcript_filtering"

    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule merge_by_condition:
    '''
    Merge assembled transcripts by condition into a single GTF file of non-redundant transcripts
    '''
    input:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
#                     "{condition}",
                     "gtf_list_{condition}.txt")

    output:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}",
                           "min_frac_{min_frac}",
                           "min_cov_{min_cov}",
                           "{condition}.all_samples.combined.gtf"),
        tracking = os.path.join(STRINGTIE_SUBDIR,
                                "min_jnc_{min_jnc}",
                                "min_frac_{min_frac}",
                                "min_cov_{min_cov}",
                                "{condition}.all_samples.tracking")

    params:
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}",
                                  "{condition}.all_samples"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "{condition}_merge_by_condition.log")

    group:
        "transcript_filtering"

    shell:
        """
        gffcompare \
        -o {params.out_prefix} \
        -p {params.label} \
        -V \
        {input} \
        2> {log}
        """


rule tx_condition_mean_tpm_filter:
    input:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}",
                           "min_frac_{min_frac}",
                           "min_cov_{min_cov}",
                           "{condition}.all_samples.combined.gtf"),
        tracking = os.path.join(STRINGTIE_SUBDIR,
                                "min_jnc_{min_jnc}",
                                "min_frac_{min_frac}",
                                "min_cov_{min_cov}",
                                "{condition}.all_samples.tracking"),
        gtf_list = os.path.join(STRINGTIE_SUBDIR,
                                "min_jnc_{min_jnc}",
                                "min_frac_{min_frac}",
                                "min_cov_{min_cov}",
#                                "{condition}",
                                "gtf_list_{condition}.txt")

    output:
        os.path.join(STRINGTIE_SUBDIR,
                            "min_jnc_{min_jnc}",
                            "min_frac_{min_frac}",
                            "min_cov_{min_cov}",
                            "{condition}.min_mean_tpm_filtered.gtf")

    params:
        script = "scripts/filter_tx_by_condition_tpm.py",
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}",
                                  "{condition}.min_mean_tpm_filtered"),
        min_tpm = config["min_mean_tpm"]


    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "{condition}_tx_condition_mean_tpm_filter.log")

    group:
        "transcript_filtering"

    shell:
        """
        python {params.script} \
        -i {input.gtf} \
        -t {input.tracking} \
        -l {input.gtf_list} \
        -m {params.min_tpm} \
        -o {params.out_prefix} 2> {log}
        """


rule gtf_list_all_tpm_filtered:
    input:
        lambda wildcards: expand(os.path.join(STRINGTIE_SUBDIR,
                                              "min_jnc_" + wildcards.min_jnc,
                                              "min_frac_" + wildcards.min_frac,
                                              "min_cov_" + wildcards.min_cov,
                                              "{condition}.min_mean_tpm_filtered.gtf"
                                              ),
                                 condition=CONDITIONS)

    output:
        txt = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}",
                           "min_frac_{min_frac}",
                           "min_cov_{min_cov}",
                           "all_conditions_tpm_filtered_gtf_list.txt")

    group:
        "transcript_filtering"

    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule merge_all_tpm_filtered:
    input:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "all_conditions_tpm_filtered_gtf_list.txt")
    output:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}",
                           "min_frac_{min_frac}",
                           "min_cov_{min_cov}",
                           "tpm_filtered.all_samples.combined.gtf"),

    params:
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}",
                                  "tpm_filtered.all_samples"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "merge_all_tpm_filtered.log")

    group:
        "transcript_filtering"

    shell:
        """
        gffcompare \
        -o {params.out_prefix} \
        -p {params.label} \
        -V \
        {input} \
        2> {log}
        """


rule intron_chain_filter:
    """
    Filter novel transcripts for those with matching intron chains to reference transcripts up until their penultimate introns (i.e. novel last exons)
    """
    input:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}",
                           "min_frac_{min_frac}",
                           "min_cov_{min_cov}",
                           "tpm_filtered.all_samples.combined.gtf")

    output:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                           "tpm_filtered.intron_chain_filtered.all_samples.combined.gtf"),
        match_stats = os.path.join(STRINGTIE_SUBDIR,
                                   "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                                   "tpm_filtered.intron_chain_filtered.all_samples.combined.match_stats.tsv")

    params:
        script = "scripts/filter_tx_by_intron_chain.py",
        ref_gtf = GTF,
        annot_source = config["annotation_source"],
        match_by = config["intron_chain_filter_mode"],
        max_terminal_non_match = config["max_terminal_non_match"],
        min_ext_length = config["min_extension_length"],
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}",
                                  "tpm_filtered.intron_chain_filtered.all_samples.combined")

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "intron_chain_filter.log")

    resources:
        threads = config["intron_chain_filter_threads"]

    group:
        "transcript_filtering"

    shell:
        """
        python {params.script} \
        -i {input} \
        -r {params.ref_gtf} \
        -s {params.annot_source} \
        -m {params.match_by} \
        -n {params.max_terminal_non_match} \
        -e {params.min_ext_length} \
        -c {resources.threads} \
        -o {params.out_prefix} \
        2> {log}
        """

rule three_end_filter:
    """
    Filter intron-chain filtered transcripts for those with
    atlas polyA sites or polyA signal motifs close to predicted 3'end
    """

    input:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                           "tpm_filtered.intron_chain_filtered.all_samples.combined.gtf"),
        match_stats = os.path.join(STRINGTIE_SUBDIR,
                                   "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                                   "tpm_filtered.intron_chain_filtered.all_samples.combined.match_stats.tsv"),
        fasta = config["genome_fasta"],
        atlas = config["polya_site_atlas"]

    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                     "tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples.combined.gtf")

    params:
        script = "scripts/filter_tx_by_three_end.py",
        motifs = config["polya_signal_motifs"],
        max_atl_dist = config["max_atlas_distance"],
        motif_len = config["motif_search_region_length"],
        output_prefix = os.path.join(STRINGTIE_SUBDIR,
                                     "min_jnc_{min_jnc}",
                                     "min_frac_{min_frac}",
                                     "min_cov_{min_cov}",
                                     "tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples.combined")

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "three_end_filter.log")

    group:
        "transcript_filtering"

    shell:
        """
        python {params.script} \
        -i {input.gtf} \
        -s {input.match_stats} \
        -f {input.fasta} \
        -a {input.atlas} \
        -p {params.motifs} \
        -m {params.max_atl_dist} \
        -u {params.motif_len} \
        -o {params.output_prefix}
        """


rule merge_filtered_with_ref:
    '''
    Merge intron_chain_filtered transcripts across samples with reference GTF files
    Into a single GTF file of non-redundant annotated & novel transcripts
    '''
    input:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples.combined.gtf")

    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "ref_merged.tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples.combined.gtf")

    params:
        ref_gtf = GTF,
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}",
                                  "ref_merged.tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "merge_filtered_with_ref.log")

    group:
        "transcript_filtering"

    shell:
        """
        gffcompare \
        -r {params.ref_gtf} \
        -o {params.out_prefix} \
        -p {params.label} \
        -V \
        {input} \
        2> {log}
        """
