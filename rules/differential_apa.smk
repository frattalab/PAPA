

rule assign_tx_to_pas:
    input:
        os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "ref_merged.tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples.combined.gtf")

    output:
        tx2pas = os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "ref_merged.tx2pas.tsv"),
        tx2gene = os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "ref_merged.tx2gene.tsv"),
        pas_le_df = os.path.join(STRINGTIE_SUBDIR,
                               "min_jnc_{min_jnc}",
                               "min_frac_{min_frac}",
                               "min_cov_{min_cov}",
                               "ref_merged.pas_assignment.tsv")

    params:
        script = "scripts/assign_tx_to_pas.py",
        merge_window = config["pas_merge_window_size"],
        group_key = "gene_id",
        tx_key = "transcript_id",
        out_prefix = os.path.join(STRINGTIE_SUBDIR,
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}",
                                  "ref_merged")

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "assign_tx_to_pas.log")

    group:
        "tx_to_per_polya_site_quant"

    shell:
        """
        python {params.script} -i {input} \
        -w {params.merge_window} \
        -g {params.group_key} \
        -t {params.tx_key} \
        -o {params.out_prefix} \
        2> {log}
        """


rule tx_to_polya_quant:
    input:
        tx2pas = os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "ref_merged.tx2pas.tsv"),
        tx2gene = os.path.join(STRINGTIE_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "ref_merged.tx2gene.tsv"),
        quant = expand(os.path.join(SALMON_SUBDIR,
                                    "quant",
                                    "min_jnc_{{min_jnc}}",
                                    "min_frac_{{min_frac}}",
                                    "min_cov_{{min_cov}}",
                                    "{sample}", "quant.sf"),
                       sample=SAMPLES)

    output:
        os.path.join(SALMON_SUBDIR,
                     "pas_quant",
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "summarised_pas_quantification.tsv")

    params:
        script = "scripts/tx_to_polya_quant.R",
        sample_tbl = config["sample_tbl"],
        salmon_dir = os.path.join(SALMON_SUBDIR,
                                  "quant",
                                  "min_jnc_{min_jnc}",
                                  "min_frac_{min_frac}",
                                  "min_cov_{min_cov}")

    conda:
        "../envs/papa_r.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     "min_jnc_{min_jnc}",
                     "min_frac_{min_frac}",
                     "min_cov_{min_cov}",
                     "tx_to_polya_quant.log")

    group:
        "tx_to_per_polya_site_quant"

    shell:
        """
        Rscript {params.script} \
        -s {params.sample_tbl} \
        -d {params.salmon_dir} \
        -t {input.tx2pas} \
        -g {input.tx2gene} \
        -o {output} \
        2> {log}
        """
