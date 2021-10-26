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
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
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
