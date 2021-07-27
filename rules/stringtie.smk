

rule stringtie:
    input:
        bam = lambda wildcards: get_bam(wildcards.sample, OPTIONS, OUTPUT_DIR)

    output:
        os.path.join(STRINGTIE_SUBDIR, "{sample}.assembled.gtf")

    params:
        gtf = GTF,
        point_feats = "--ptf " + config["polya_site_point_features"] if config["use_point_features"] else "",
        strandedness = config["strandedness"],
        label = config["label"],
        min_iso_frac = config["min_isoform_fraction_abundance"],
        min_iso_len = config["min_isoform_length"],
        gene_abund = lambda wildcards: " ".join(["-A", os.path.join(STRINGTIE_SUBDIR, wildcards.sample + config["gene_abundances_suffix"])]) if config["report_gene_abundances"] else "",
        annot_tr = lambda wildcards: " ".join(["-C", os.path.join(STRINGTIE_SUBDIR, wildcards.sample + config["covered_txipts_suffix"])]) if config["report_covered_annot_txipts"] else "",
        min_jnc_ohang = config["min_junction_overhang"],
        min_jnc_reads = config["min_junction_reads"],
        trimming = "-t" if config["disable_end_trimming"] else "",
        min_cov = config["min_txipt_coverage"],
        min_se_cov = config["min_single_exon_coverage"],
        conservative = "--conservative" if config["conservative_mode"] else "",
        min_locus_gap = config["min_locus_gap"],
        max_multimap_frac = config["max_fraction_multi_mapped"]


    conda:
        "../envs/papa.yaml"

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
        -o {output}
        """

rule extract_novel_stringtie:
    """
    Filter for assembled transcripts that do not match reference transcripts (i.e. extract novel isoforms)
    """
    input:
        os.path.join(STRINGTIE_SUBDIR, "{sample}.assembled.gtf")

    output:
        os.path.join(STRINGTIE_SUBDIR, "{sample}.no_ref_id.assembled.gtf")

    params:
        ref_string = config["stringtie_ref_string"]
    shell:
        """
        grep -v '{params.ref_string}' {input} > {output}
        """

rule intron_chain_filter:
    """
    Filter novel transcripts for those with matching intron chains to reference transcripts up until their penultimate introns (i.e. novel last exons)
    """
    input:
        os.path.join(STRINGTIE_SUBDIR, "{sample}.no_ref_id.assembled.gtf")

    output:
        os.path.join(STRINGTIE_SUBDIR, "{sample}.intron_chain_filtered.no_ref_id.assembled.gtf")

    params:
        script = "scripts/filter_tx_by_intron_chain.py",
        ref_gtf = GTF,
        match_by = config["intron_chain_filter_mode"],
        max_terminal_non_match = config["max_terminal_non_match"]

    conda:
        "../envs/papa.yaml"

    resources:
        threads = 4

    shell:
        """
        python {params.script} \
        -i {input} \
        -r {params.ref_gtf} \
        -m {params.match_by} \
        -n {params.max_terminal_non_match} \
        -c {resources.threads} \
        -o {output}
        """


rule compose_gtf_list_stringtie:
    input:
        expand(os.path.join(STRINGTIE_SUBDIR, "{sample}.intron_chain_filtered.no_ref_id.assembled.gtf"), sample = SAMPLES)
    output:
        txt = os.path.join(STRINGTIE_SUBDIR,"gtf_list.txt")
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule stringtie_merge:
    input:
        os.path.join(STRINGTIE_SUBDIR, "gtf_list.txt")

    output:
        os.path.join(STRINGTIE_SUBDIR, "all_samples.intron_chain_filtered.merged.gtf")

    params:
        gtf = GTF,
        min_len = config["min_length_merge"],
        min_cov = config["min_cov_merge"],
        min_fpkm = config["min_fpkm_merge"],
        min_tpm = config["min_tpm_merge"],
        min_frac = config["min_iso_frac_merge"],
        keep_ri = "-i" if config["keep_retained_introns_merge"] else "",
        label = config["label"]

    conda:
        "../envs/papa.yaml"

    shell:
        """
        stringtie --merge \
        -G {params.gtf} \
        -m {params.min_len} \
        -c {params.min_cov} \
        -F {params.min_fpkm} \
        -T {params.min_tpm} \
        -f {params.min_frac} \
        {params.keep_ri} \
        -l {params.label} \
        -o {output} \
        {input}
        """
