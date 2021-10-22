# wildcard_constraints:
#     sample = "|".join(SAMPLES),
#     min_frac = "min_frac_\\d$",
#     min_jnc = "min_jnc_\\d$",
#     min_cov = "min_cov_\\d$"


rule stringtie:
    input:
        bam = lambda wildcards: get_bam(wildcards.sample, OPTIONS, OUTPUT_DIR)

    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                     "{sample}.assembled.gtf")

    params:
        gtf = GTF,
        point_feats = "--ptf " + config["polya_site_point_features"] if config["use_point_features"] else "",
        strandedness = config["strandedness"],
        label = config["label"],
        min_iso_frac = "{min_frac}",
        min_iso_len = config["min_isoform_length"],
        gene_abund = lambda wildcards: " ".join(["-A", os.path.join(STRINGTIE_SUBDIR, wildcards.sample + config["gene_abundances_suffix"])]) if config["report_gene_abundances"] else "",
        annot_tr = lambda wildcards: " ".join(["-C", os.path.join(STRINGTIE_SUBDIR, wildcards.sample + config["covered_txipts_suffix"])]) if config["report_covered_annot_txipts"] else "",
        min_jnc_ohang = config["min_junction_overhang"],
        min_jnc_reads = "{min_jnc}",
        trimming = "-t" if config["disable_end_trimming"] else "",
        min_cov = "{min_cov}",
        min_se_cov = config["min_single_exon_coverage"],
        conservative = "--conservative" if config["conservative_mode"] else "",
        min_locus_gap = config["min_locus_gap"],
        max_multimap_frac = config["max_fraction_multi_mapped"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "{sample}.stringtie_assemble.log")

    shell:
        """
        stringtie {input.bam} \
        -G {params.gtf} \
        {params.strandedness} \
        {params.point_feats} \
        -l {params.label} \
        -f {params.min_iso_frac} \
        {params.gene_abund} \
        {params.annot_tr} \
        -a {params.min_jnc_ohang} \
        -j {params.min_jnc_reads} \
        {params.trimming} \
        -c {params.min_cov} \
        -s {params.min_se_cov} \
        {params.conservative} \
        -g {params.min_locus_gap} \
        -M {params.max_multimap_frac} \
        -o {output} \
        2> {log}
        """

# rule extract_novel_stringtie:
#     """
#     Filter for assembled transcripts that do not match reference transcripts (i.e. extract novel isoforms)
#     """
#     input:
#         os.path.join(STRINGTIE_SUBDIR, "{sample}.assembled.gtf")
#
#     output:
#         os.path.join(STRINGTIE_SUBDIR, "{sample}.no_ref_id.assembled.gtf")
#
#     params:
#         ref_string = config["stringtie_ref_string"]
#
#     log:
#         os.path.join(LOG_SUBDIR, "{sample}.extract_novel_stringtie.log")
#
#     shell:
#         """
#         grep -v '{params.ref_string}' {input} > {output} 2> {log}
#         """

rule intron_chain_filter:
    """
    Filter novel transcripts for those with matching intron chains to reference transcripts up until their penultimate introns (i.e. novel last exons)
    """
    input:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                     "{sample}.assembled.gtf")

    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
                     "{sample}.intron_chain_filtered.assembled.gtf")

    params:
        script = "scripts/filter_tx_by_intron_chain.py",
        ref_gtf = GTF,
        annot_source = config["annotation_source"],
        match_by = config["intron_chain_filter_mode"],
        max_terminal_non_match = config["max_terminal_non_match"],
        min_ext_length = config["min_extension_length"],
        out_prefix = os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "{sample}.intron_chain_filtered.assembled")

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "{sample}.intron_chain_filter.log")

    resources:
        threads = config["intron_chain_filter_threads"]

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


rule compose_gtf_list_stringtie:
    input:
        lambda wildcards: expand(os.path.join(STRINGTIE_SUBDIR,
                                              "min_jnc_" + wildcards.min_jnc,
                                              "min_frac_" + wildcards.min_frac,
                                              "min_cov_" + wildcards.min_cov,
                                              "{sample}.intron_chain_filtered.assembled.gtf"
                                              ),
                                 sample=SAMPLES)
        #
        # expand(os.path.join(STRINGTIE_SUBDIR, "{{min_jnc}}_{{min_frac}}_{{min_cov}}", "{sample}.intron_chain_filtered.assembled.gtf"),
        #        sample=SAMPLES,
        #        min_jnc=min_jnc_vals,
        #        min_frac=min_frac_vals,
        #        min_cov=min_cov_vals)
    output:
        txt = os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}_min_frac_{min_frac}_min_cov_{min_cov}_gtf_list.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule gtf_merge_novel:
    '''
    Merge intron_chain_filtered transcripts only into a single GTF file of non-redundant NOVEL transcripts
    This set does not include any reference transcripts
    '''
    input:
        os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}_min_frac_{min_frac}_min_cov_{min_cov}_gtf_list.txt")

    output:
        os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "all_samples.intron_chain_filtered.novel.combined.gtf")

    params:
        out_prefix = os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "all_samples.intron_chain_filtered.novel"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "stringtie_merge_novel.log")

    shell:
        """
        gffcompare \
        -o {params.out_prefix} \
        -p {params.label} \
        -V \
        {input} \
        2> {log}
        """

rule gtf_merge_ref:
    '''
    Merge intron_chain_filtered transcripts across samples with reference GTF files
    Into a single GTF file of non-redundant annotated & novel transcripts
    '''
    input:
        os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}_min_frac_{min_frac}_min_cov_{min_cov}_gtf_list.txt")

    output:
        os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "all_samples.intron_chain_filtered.ref_merged.combined.gtf")

    params:
        gtf = GTF,
        out_prefix = os.path.join(STRINGTIE_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "all_samples.intron_chain_filtered.ref_merged"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "stringtie_merge_novel.log")

    shell:
        """
        gffcompare \
        -r {params.gtf} \
        -o {params.out_prefix} \
        -p {params.label} \
        -V \
        {input} \
        2> {log}
        """


# rule stringtie_merge_novel:
#     '''
#     Merge intron_chain_filtered transcripts only into a single GTF file of non-redundant NOVEL transcripts
#     This set does not include any reference transcripts
#     '''
#         input:
#             os.path.join(STRINGTIE_SUBDIR, "gtf_list.txt")
#
#         output:
#             os.path.join(STRINGTIE_SUBDIR, "all_samples.intron_chain_filtered.merged.gtf")
#
#         params:
#             min_len = config["min_length_merge"],
#             min_cov = config["min_cov_merge"],
#             min_fpkm = config["min_fpkm_merge"],
#             min_tpm = config["min_tpm_merge"],
#             min_frac = config["min_iso_frac_merge"],
#             keep_ri = "-i" if config["keep_retained_introns_merge"] else "",
#             label = config["label"]
#
#         conda:
#             "../envs/papa.yaml"
#
#         log:
#             os.path.join(LOG_SUBDIR, "stringtie_merge_novel.log")
#
#         shell:
#             """
#             stringtie --merge \
#             -m {params.min_len} \
#             -c {params.min_cov} \
#             -F {params.min_fpkm} \
#             -T {params.min_tpm} \
#             -f {params.min_frac} \
#             {params.keep_ri} \
#             -l {params.label} \
#             -o {output} \
#             {input} \
#             2> {log}
#             """
#
#
# rule stringtie_merge_ref:
#     input:
#         os.path.join(STRINGTIE_SUBDIR, "gtf_list.txt")
#
#     output:
#         os.path.join(STRINGTIE_SUBDIR, "all_samples.intron_chain_filtered.ref_merged.gtf")
#
#     params:
#         gtf = GTF,
#         min_len = config["min_length_merge"],
#         min_cov = config["min_cov_merge"],
#         min_fpkm = config["min_fpkm_merge"],
#         min_tpm = config["min_tpm_merge"],
#         min_frac = config["min_iso_frac_merge"],
#         keep_ri = "-i" if config["keep_retained_introns_merge"] else "",
#         label = config["label"]
#
#     conda:
#         "../envs/papa.yaml"
#
#     log:
#         os.path.join(LOG_SUBDIR, "stringtie_merge_ref.log")
#
#     shell:
#         """
#         stringtie --merge \
#         -G {params.gtf} \
#         -m {params.min_len} \
#         -c {params.min_cov} \
#         -F {params.min_fpkm} \
#         -T {params.min_tpm} \
#         -f {params.min_frac} \
#         {params.keep_ri} \
#         -l {params.label} \
#         -o {output} \
#         {input} \
#         2> {log}
#         """
