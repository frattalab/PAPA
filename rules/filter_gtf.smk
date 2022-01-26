rule filter_ref_gtf:
    input:
        GTF

    output:
        os.path.join(OUTPUT_DIR, "reference_filtered.gtf")

    params:
        script = "scripts/filter_gtf.py",
        min_tsl = config["min_transcript_support_level"],
        include_flags = parse_filter_flags(config["ref_gtf_include_flags"], "--include-flags"),
        exclude_flags = parse_filter_flags(config["ref_gtf_exclude_flags"], "--exclude-flags"),
        include_tags = parse_filter_tags(config["include_tags"], "--tag-include"),
        exclude_tags = parse_filter_tags(config["exclude_tags"], "--tag-exclude")

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR, "filter_ref_gtf.txt")

    log:
        os.path.join(LOG_SUBDIR, "filter_ref_gtf.log")

    shell:
        """
        python {params.script} \
        -i {input} \
        -t {params.min_tsl} \
        {params.include_flags} \
        {params.exclude_flags} \
        {params.include_tags} \
        {params.exclude_tags} \
        -o {output} &> {log}
        """
