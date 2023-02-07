rule gtf_list_by_condition:
    input:
        lambda wildcards: expand(os.path.join(STRINGTIE_SUBDIR,
                                              "{sample}.assembled.gtf"
                                              ),
                                 sample=get_condition_samples(wildcards.condition, OPTIONS))

    output:
        txt = os.path.join(TX_FILT_SUBDIR,
                           "{condition}",
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

    output:
        gtf = temp(os.path.join(TX_FILT_SUBDIR,
                                "{condition}",
                                "{condition}.all_samples.combined.gtf")),
        tracking = os.path.join(TX_FILT_SUBDIR,
                                "{condition}",
                                "{condition}.all_samples.tracking")

    params:
        out_prefix = os.path.join(TX_FILT_SUBDIR,
                                  "{condition}",
                                  "{condition}.all_samples"),
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{condition}_merge_by_condition.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
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
    The script will output per-sample files, but a downstream dummy rule will check these have been created
    '''
    input:
        tracking = os.path.join(TX_FILT_SUBDIR,
                                "{condition}",
                                "{condition}.all_samples.tracking"),
        gtf_list = os.path.join(TX_FILT_SUBDIR,
                                "{condition}",
                                "gtf_list_{condition}.txt"),

    output:
        dummy = temp(os.path.join(TX_FILT_SUBDIR,
                                  "{condition}",
                                  "{condition}_mean_tpm.done")),

    params:
        script = "scripts/filter_tx_by_condition_tpm.py",
        out_suffix = ".mean_tpm_filtered.gtf",
        min_tpm = config["min_mean_tpm"]


    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{condition}_tx_condition_mean_tpm_filter.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
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
        lambda wildcards: os.path.join(TX_FILT_SUBDIR,
                                       OPTIONS[wildcards.sample]["condition"],
                                       OPTIONS[wildcards.sample]["condition"] + "_mean_tpm.done")
    output:
        os.path.join(STRINGTIE_SUBDIR,
                     "mv.{sample}.assembled.mean_tpm_filtered.gtf")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
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
        os.path.join(TX_FILT_SUBDIR,
                     "{condition}",
                     "{sample}.last_exons.gtf")

    params:
        script = "scripts/get_novel_last_exons.py",
        min_ext_length = config["min_extension_length"],
        first_ex_5p_tolerance = config["first_exon_5p_tolerance"],
        other_ex_5p_tolerance = config["other_exon_5p_tolerance"],
        trust_input_exon_number = "",
        trust_ref_exon_number = "",
        ignore_ext_tol = "" if config["extension_tolerance_filter"] else "--ignore-extension-tolerance",
        ignore_spl_tol = "" if config["spliced_tolerance_filter"] else "--ignore-spliced-tolerance",
        out_prefix = os.path.join(TX_FILT_SUBDIR,
                                  "{condition}",
                                  "{sample}")

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{sample}.get_novel_last_exons.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{sample}.get_novel_last_exons.log")

    # group:
    #     "transcript_filtering_events"

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
        {params.ignore_ext_tol} \
        {params.ignore_spl_tol} \
        -o {params.out_prefix} &> {log}
        """


rule combine_novel_by_condition:
    '''
    Concatenate last exons across samples of the same condition, grouping overlapping last exons across samples with a common identifier
    '''
    input:
        gtfs = lambda wildcards: expand(os.path.join(TX_FILT_SUBDIR,
                                                     wildcards.condition,
                                                     "{sample}.last_exons.gtf"
                                                     ),
                                        sample=get_condition_samples(wildcards.condition, OPTIONS))

    output:
        os.path.join(TX_FILT_SUBDIR,
                     "{condition}",
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
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{condition}_combine_novel_by_condition.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{condition}_combine_novel_by_condition.log")

    # group:
    #     "transcript_filtering_3p"

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
        gtf = os.path.join(TX_FILT_SUBDIR,
                           "{condition}",
                           "{condition}.merged_last_exons.gtf"
                           ),
        fasta = config["genome_fasta"],
        atlas = config["polya_site_atlas"]

    output:
        gtf = os.path.join(TX_FILT_SUBDIR,
                           "{condition}",
                           "{condition}.merged_last_exons.3p_end_filtered.gtf"),

    params:
        script = "scripts/filter_tx_by_three_end.py",
        motifs = config["polya_signal_motifs"],
        max_atl_dist = config["max_atlas_distance"],
        motif_len = config["motif_search_region_length"],
        motif_exp_dist = config["motif_expected_distance"],
        output_prefix = os.path.join(TX_FILT_SUBDIR,
                                     "{condition}",
                                     "{condition}.merged_last_exons.3p_end_filtered")

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{condition}_three_end_filter.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "{condition}",
                     "{condition}_three_end_filter.log")

    # group:
    #     "transcript_filtering_3p"

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


def filtered_novel_targets(conditions, three_end_filter):
    '''
    Returns target filtered novel last exon GTFs depending on whether applying 3'end filter or not
    '''

    assert isinstance(three_end_filter, bool)

    if three_end_filter:
        return expand(os.path.join(TX_FILT_SUBDIR,
                                   "{condition}",
                                   "{condition}.merged_last_exons.3p_end_filtered.gtf"
                                   ),
                      condition=conditions)

    else:
        # No 3'end filtering
        return expand(os.path.join(TX_FILT_SUBDIR,
                                   "{condition}",
                                   "{condition}.merged_last_exons.gtf"
                                   ),
                      condition=conditions
                      )


rule combine_novel_filtered_by_condition:
    '''
    Combine 3'end filtered last exons across conditions into a single GTF file
    '''
    input:
        filtered_novel_targets(CONDITIONS, config["three_end_filter"])

    output:
        gtf = os.path.join(TX_FILT_SUBDIR,
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
                     config["tx_filtering_subdir_name"],
                     "combine_novel_filtered_by_condition.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "combine_novel_filtered_by_condition.log")

    # group:
    #     "transcript_filtering_3p"

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
        novel_gtf = rules.combine_novel_filtered_by_condition.output.gtf if not config["use_provided_novel_les"] else config["input_novel_gtf"],
        ref_gtf = rules.filter_ref_gtf.output if config["filter_ref_gtf"] else GTF

    output:
        quant_gtf = os.path.join(TX_FILT_SUBDIR,
                                 "novel_ref_combined.quant.last_exons.gtf"),
        le_gtf = os.path.join(TX_FILT_SUBDIR, "novel_ref_combined.last_exons.gtf"),
        tx2le = os.path.join(TX_FILT_SUBDIR, "novel_ref_combined.tx2le.tsv"),
        tx2gene = os.path.join(TX_FILT_SUBDIR, "novel_ref_combined.tx2gene.tsv"),
        le2gene = os.path.join(TX_FILT_SUBDIR, "novel_ref_combined.le2gene.tsv"),
        le2genename = os.path.join(TX_FILT_SUBDIR, "novel_ref_combined.le2genename.tsv"),
        info_tbl = os.path.join(TX_FILT_SUBDIR, "novel_ref_combined.info.tsv")

    params:
        script = "scripts/get_combined_quant_gtf.py",
        output_prefix = os.path.join(TX_FILT_SUBDIR,
                                     "novel_ref_combined"),
        trust_ref_exon_number = "",
        ref_extensions_string = "" if len(config["ref_extensions_string"]) == 0 else "--ref-extensions-string " + config["ref_extensions_string"]

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "get_combined_quant_gtf.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "get_combined_quant_gtf.log")

    shell:
        """
        python {params.script} \
        -i {input.novel_gtf} \
        -r {input.ref_gtf} \
        {params.trust_ref_exon_number} \
        {params.ref_extensions_string} \
        -o {params.output_prefix} \
        &> {log}
        """

rule get_ref_quant_gtf:
    '''
    '''
    input:
        gtf = rules.filter_ref_gtf.output if config["filter_ref_gtf"] else GTF

    output:
        quant_gtf = os.path.join(TX_FILT_SUBDIR,
                                 "ref.quant.last_exons.gtf"),
        le_gtf = os.path.join(TX_FILT_SUBDIR, "ref.last_exons.gtf"),
        tx2le = os.path.join(TX_FILT_SUBDIR, "ref.tx2le.tsv"),
        tx2gene = os.path.join(TX_FILT_SUBDIR, "ref.tx2gene.tsv"),
        le2gene = os.path.join(TX_FILT_SUBDIR, "ref.le2gene.tsv"),
        le2genename = os.path.join(TX_FILT_SUBDIR, "ref.le2genename.tsv"),
        info_tbl = os.path.join(TX_FILT_SUBDIR, "ref.info.tsv")

    params:
        script = "scripts/get_ref_quant_gtf.py",
        output_prefix = os.path.join(TX_FILT_SUBDIR,
                                     "ref"),
        trust_ref_exon_number = "",
        ref_extensions_string = "" if len(config["ref_extensions_string"]) == 0 else "--ref-extensions-string " + config["ref_extensions_string"]

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "get_ref_quant_gtf.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["tx_filtering_subdir_name"],
                     "get_ref_quant_gtf.log")

    shell:
        """
        python {params.script} \
        -i {input.gtf} \
        {params.trust_ref_exon_number} \
        {params.ref_extensions_string} \
        -o {params.output_prefix} \
        &> {log}
        """
