rule gtf_list_by_condition:
    input:
        lambda wildcards: expand(os.path.join(STRINGTIE_SUBDIR,
                                              "{sample}.assembled.gtf"
                                              ),
                                 sample=get_condition_samples(wildcards.condition, OPTIONS))

    output:
        txt = os.path.join(STRINGTIE_SUBDIR,
#                           "{condition}",
                           "gtf_list_{condition}.txt")

    group:
        "transcript_filtering_tpm"

    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule merge_by_condition:
    '''
    Merge assembled transcripts by condition into a single GTF file of non-redundant transcripts
    '''
    input:
        rules.gtf_list_by_condition.output.txt
        # os.path.join(STRINGTIE_SUBDIR,
        #              "min_jnc_{min_jnc}",
        #              "min_frac_{min_frac}",
        #              "min_cov_{min_cov}",
        #              "gtf_list_{condition}.txt")

    output:
        gtf = temp(os.path.join(STRINGTIE_SUBDIR,
                           "{condition}.all_samples.combined.gtf")),
        tracking = os.path.join(STRINGTIE_SUBDIR,
                                "{condition}.all_samples.tracking")

    params:
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "{condition}.all_samples"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "{condition}_merge_by_condition.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "{condition}_merge_by_condition.log")

    group:
        "transcript_filtering_tpm"

    shell:
        """
        gffcompare \
        -o {params.out_prefix} \
        -p {params.label} \
        -V \
        -i {input} \
        &> {log}
        """


rule tx_condition_mean_tpm_filter:
    '''
    Calculate mean TPM for transcripts across samples of the same condition
    Apply a minimum threshold to retain transcripts
    The script works with a .tracking' file output by GFFcompare to map tx IDs to individual samples
    Then goes back to individual sample GTFs and filters them for tx_ids passing min threshold (stripping only the '.gtf' from the filename in input list)
    Merged GTFs take the longest predicted 3'end for structurally matched transcripts across replicates, but want to consider as many putative 3'ends as possible

    Note: the '.done' dummy file is a hacky solution to ensure that the rule runs once per condition,
    as I couldn't find a way to propagate the 'sample' wildcard (but to only expect 'sample' belonging to the condition)
    The script will output per-sample files, but a downstream dummy rule will check this
    '''
    input:
        tracking = os.path.join(STRINGTIE_SUBDIR,
                                "{condition}.all_samples.tracking"),
        gtf_list = os.path.join(STRINGTIE_SUBDIR,
                                "gtf_list_{condition}.txt"),

    output:
        dummy = temp(os.path.join(STRINGTIE_SUBDIR, "{condition}_mean_tpm.done")),
        # expand(os.path.join(STRINGTIE_SUBDIR,
        #                     "{condition}",
        #                     "{sample}.mean_tpm_filtered.gtf"
        #                     ),
        #        sample=get_condition_samples("{condition}", OPTIONS),
        #        )

    # wildcard_constraints:
    #     sample = get_condition_samples("{condition}", OPTIONS)

    params:
        script = "scripts/filter_tx_by_condition_tpm.py",
        out_suffix = ".mean_tpm_filtered.gtf",
        min_tpm = config["min_mean_tpm"]


    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "{condition}_tx_condition_mean_tpm_filter.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "{condition}_tx_condition_mean_tpm_filter.log")

    group:
        "transcript_filtering_tpm"

    shell:
        """
        python {params.script} \
        -t {input.tracking} \
        -l {input.gtf_list} \
        -m {params.min_tpm} \
        -o {params.out_suffix} \
        &> {log} && touch {output.dummy}
        """

rule check_per_sample_mean_tpm_filtered:
    '''
    Dummy rule to ensure per-sample mean tpm filtered GTF files have been generated
    Simply exits with error code if expected file doesn't exist
    '''
    input:
        lambda wildcards: os.path.join(STRINGTIE_SUBDIR,
                                       OPTIONS[wildcards.sample]["condition"] + "_mean_tpm.done")
    output:
        os.path.join(STRINGTIE_SUBDIR, "mv.{sample}.assembled.mean_tpm_filtered.gtf")

    log:
        os.path.join(LOG_SUBDIR,
                     "{sample}_check_per_sample_mean_tpm_filtered.log")
    shell:
        """
        # Remove the '.mv' from output name - this is file output by previous rule
        filt_fname=$(echo $(dirname {output})/$(basename {output} | cut -d'.' -f2-))

        if [ ! -e $filt_fname ]; then
            echo Expected mean TPM filtered sample file was not generated - $filt_fname - check rule tx_condition_mean_tpm_filter
            exit 1
        else
            echo Success! Expected mean TPM filtered sample file was generated by rule tx_condition_mean_tpm_filter - $filt_fname
            mv $filt_fname {output}
            exit 0
        fi &> {log}
        """

rule get_novel_last_exons:
    '''
    Extract candidate novel last exons for each sample of TPM filtered predicted transcripts
    '''
    input:
        novel_gtf = os.path.join(STRINGTIE_SUBDIR,
                                 "mv.{sample}.assembled.mean_tpm_filtered.gtf"
                                 ),
        ref_gtf = rules.filter_ref_gtf.output if config["filter_ref_gtf"] else GTF

    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "{condition}",
                     "{sample}.last_exons.gtf")

    params:
        script = "scripts/get_novel_last_exons.py",
        min_ext_length = config["min_extension_length"],
        first_ex_5p_tolerance = config["first_exon_5p_tolerance"],
        other_ex_5p_tolerance = config["other_exon_5p_tolerance"],
        trust_input_exon_number = "",
        trust_ref_exon_number = "",
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "{condition}",
                                  "{sample}")

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "{condition}",
                     "{sample}.get_novel_last_exons.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "{condition}",
                     "{sample}.get_novel_last_exons.log")

    group:
        "transcript_filtering_tpm"

    shell:
        """
        python {params.script} \
        -i {input.novel_gtf} \
        -r {input.ref_gtf} \
        -m {params.min_ext_length} \
        -f {params.first_ex_5p_tolerance} \
        -t {params.other_ex_5p_tolerance} \
        {params.trust_input_exon_number} \
        {params.trust_ref_exon_number} \
        -o {params.out_prefix} &> {log}
        """


rule combine_novel_by_condition:
    '''
    Concatenate last exons across samples of the same condition, grouping overlapping last exons across samples with a common identifier
    '''
    input:
        gtfs = lambda wildcards: expand(os.path.join(STRINGTIE_SUBDIR,
                                                     wildcards.condition,
                                                     "{sample}.last_exons.gtf"
                                                     ),
                                        sample = get_condition_samples(wildcards.condition, OPTIONS))

    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "{condition}.merged_last_exons.gtf"
                     )

    params:
        script = "scripts/combine_novel_last_exons.py",
        gtf_suffix = ".last_exons.gtf",
        sample_id_key = "sample_id",
        last_exon_id_key = "last_exon_id"

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "{condition}_combine_novel_by_condition.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "{condition}_combine_novel_by_condition.log")

    group:
        "transcript_filtering_tpm"

    shell:
        """
        python {params.script} \
        -i {input.gtfs} \
        -s {params.gtf_suffix} \
        --sample-id-attribute-name {params.sample_id_key} \
        --last-exon-id-attribute-name {params.last_exon_id_key} \
        -o {output} &> {log}
        """


rule three_end_filter:
    """
    Filter candidate novel last exons for those with
    atlas polyA sites or polyA signal motifs close to predicted 3'end
    """

    input:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "{condition}.merged_last_exons.gtf"
                           ),
        fasta = config["genome_fasta"],
        atlas = config["polya_site_atlas"]

    output:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "{condition}.merged_last_exons.3p_end_filtered.gtf"),

    params:
        script = "scripts/filter_tx_by_three_end.py",
        motifs = config["polya_signal_motifs"],
        max_atl_dist = config["max_atlas_distance"],
        motif_len = config["motif_search_region_length"],
        motif_exp_dist = config["motif_expected_distance"],
        output_prefix = os.path.join(STRINGTIE_SUBDIR,
                                     "{condition}.merged_last_exons.3p_end_filtered")

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "{condition}_three_end_filter.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "{condition}_three_end_filter.log")

    group:
        "transcript_filtering_chain_3p"

    shell:
        """
        python {params.script} \
        -i {input.gtf} \
        -f {input.fasta} \
        -a {input.atlas} \
        -m {params.motifs} \
        -d {params.max_atl_dist} \
        -e {params.motif_exp_dist} \
        -u {params.motif_len} \
        -o {params.output_prefix} \
        &> {log}
        """


rule combine_novel_filtered_by_condition:
    '''
    Combine 3'end filtered last exons across conditions into a single GTF file
    '''
    input:
        expand(os.path.join(STRINGTIE_SUBDIR,
                            "{condition}.merged_last_exons.3p_end_filtered.gtf"),
                            condition=CONDITIONS)

    output:
        gtf = os.path.join(STRINGTIE_SUBDIR,
                           "all_conditions.merged_last_exons.3p_end_filtered.gtf")

    params:
        script = "scripts/combine_novel_last_exons.py",
        gtf_suffix = ".merged_last_exons.3p_end_filtered.gtf",
        sample_id_col = "condition_id",
        no_le_id = "--no-last-exon-id"

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "combine_novel_filtered_by_condition.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "combine_novel_filtered_by_condition.log")

    group:
        "transcript_filtering_chain_3p"

    shell:
        """
        python {params.script} \
        -i {input} \
        -s {params.gtf_suffix} \
        --sample-id-attribute-name {params.sample_id_col} \
        {params.no_le_id} \
        -o {output.gtf} \
        &> {log}
        """


rule get_combined_quant_gtf:
    '''
    Merge novel last exons with reference last exons
    Group transcripts together according to shared last exon
    Extract unique regions for last exons overlapping first/internal reference exons
    '''
    input:
        novel_gtf = rules.combine_novel_filtered_by_condition.output.gtf,
        ref_gtf = rules.filter_ref_gtf.output if config["filter_ref_gtf"] else GTF

    output:
        quant_gtf = os.path.join(STRINGTIE_SUBDIR,
                                 "novel_ref_combined.quant.last_exons.gtf"),
        le_gtf = os.path.join(STRINGTIE_SUBDIR, "novel_ref_combined.last_exons.gtf"),
        tx2le = os.path.join(STRINGTIE_SUBDIR, "novel_ref_combined.tx2le.tsv"),
        tx2gene = os.path.join(STRINGTIE_SUBDIR, "novel_ref_combined.tx2gene.tsv"),
        le2gene = os.path.join(STRINGTIE_SUBDIR, "novel_ref_combined.le2gene.tsv")

    params:
        script = "scripts/get_combined_quant_gtf.py",
        output_prefix = os.path.join(STRINGTIE_SUBDIR,
                                     "novel_ref_combined"),
        trust_ref_exon_number = ""

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     "get_combined_quant_gtf.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     "get_combined_quant_gtf.log")

    shell:
        """
        python {params.script} \
        -i {input.novel_gtf} \
        -r {input.ref_gtf} \
        {params.trust_ref_exon_number} \
        -o {params.output_prefix} \
        &> {log}
        """



# rule intron_chain_filter:
#     """
#     Filter novel transcripts for those with matching intron chains to reference transcripts up until their penultimate introns (i.e. novel last exons)
#     """
#     input:
#         gtf = os.path.join(STRINGTIE_SUBDIR,
#                            "min_jnc_{min_jnc}",
#                            "min_frac_{min_frac}",
#                            "min_cov_{min_cov}",
#                            "tpm_filtered.all_samples.combined.gtf")
#
#     output:
#         gtf = os.path.join(STRINGTIE_SUBDIR,
#                            "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
#                            "tpm_filtered.intron_chain_filtered.all_samples.combined.gtf"),
#         match_stats = os.path.join(STRINGTIE_SUBDIR,
#                                    "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}",
#                                    "tpm_filtered.intron_chain_filtered.all_samples.combined.match_stats.tsv")
#
#     params:
#         script = "scripts/filter_tx_by_intron_chain.py",
#         ref_gtf = GTF,
#         annot_source = config["annotation_source"],
#         match_by = config["intron_chain_filter_mode"],
#         max_terminal_non_match = config["max_terminal_non_match"],
#         min_ext_length = config["min_extension_length"],
#         out_prefix = os.path.join(STRINGTIE_SUBDIR,
#                                   "min_jnc_{min_jnc}",
#                                   "min_frac_{min_frac}",
#                                   "min_cov_{min_cov}",
#                                   "tpm_filtered.intron_chain_filtered.all_samples.combined")
#
#     conda:
#         "../envs/papa.yaml"
#
#     log:
#         os.path.join(LOG_SUBDIR,
#                      "min_jnc_{min_jnc}",
#                      "min_frac_{min_frac}",
#                      "min_cov_{min_cov}",
#                      "intron_chain_filter.log")
#
#     resources:
#         threads = config["intron_chain_filter_threads"]
#
#     group:
#         "transcript_filtering_chain_3p"
#
#     shell:
#         """
#         # remove undefined strand rows from novel GTF
#         awk '{{if ($7!=".") {{print $0}} }}' {input} > {input}.tmp
#
#         python {params.script} \
#         -i {input}.tmp \
#         -r {params.ref_gtf} \
#         -s {params.annot_source} \
#         -m {params.match_by} \
#         -n {params.max_terminal_non_match} \
#         -e {params.min_ext_length} \
#         -c {resources.threads} \
#         -o {params.out_prefix} \
#         2> {log}
#
#         rm {input}.tmp
#         """
