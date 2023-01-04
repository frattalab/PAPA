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
        counts = os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.counts.tsv"),
        tpm = os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.tpm.tsv"),
        ppau = os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.ppau.tsv"),
        gene_tpm = os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.gene_tpm.tsv"),

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


rule saturn_apa:
    input:
        counts = rules.tx_to_le_quant.output.counts,
        le2gene = rules.get_combined_quant_gtf.output.le2gene,
        sample_tbl = config["sample_tbl"]

    output:
        tbl = os.path.join(DAPA_SUBDIR,
                            "saturn_apa.results.tsv"),
        rda = os.path.join(DAPA_SUBDIR,
                           "saturn_apa.image.RData")

    params:
        script = "scripts/run_differential_usage.R",
        output_prefix = os.path.join(DAPA_SUBDIR, "saturn_apa"),
        min_mean_count = config["min_mean_count"]
    
    threads:
        config["saturn_threads"]

    conda:
        "../envs/papa_r.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "saturn_apa.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "saturn_apa.txt")   

    shell:
        """
        Rscript {params.script} \
        -i {input.counts} \
        -g {input.le2gene} \
        -s {input.sample_tbl} \
        -c {threads} \
        --min-mean-count {params.min_mean_count} \
        -o {params.output_prefix} \
        &> {log}
        """
