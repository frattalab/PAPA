rule tx_to_le_quant:
    input:
        tx2le = rules.get_combined_quant_gtf.output.tx2le,
        tx2gene = rules.get_combined_quant_gtf.output.tx2gene,
        quant = expand(os.path.join(SALMON_SUBDIR,
                                    "quant",
                                    "{sample}",
                                    "quant.sf"),
                       sample=SAMPLES)

    output:
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.counts.tsv")

    params:
        script = "scripts/tx_to_polya_quant.R",
        sample_tbl = config["sample_tbl"],
        salmon_dir = os.path.join(SALMON_SUBDIR,
                                  "quant"),
        output_prefix = os.path.join(DAPA_SUBDIR,
                                     "summarised_pas_quantification")

    conda:
        "../envs/papa_r.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "tx_to_polya_quant.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "tx_to_le_quant.txt")

    shell:
        """
        Rscript {params.script} \
        -s {params.sample_tbl} \
        -d {params.salmon_dir} \
        -t {input.tx2le} \
        -g {input.tx2gene} \
        -o {params.output_prefix} \
        &> {log}
        """
