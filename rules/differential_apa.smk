def id2id_target(run_identification: bool, use_precomputed_salmon_index: bool, tx2id: str):
    '''
    Helper function to return path to id2id table/info table given flags to control workflow specification
    '''

    assert tx2id in ["tx2le", "tx2gene", "le2gene", "info_tbl"]

    if use_precomputed_salmon_index:
        # no matter what run_id value is, want to use precomputed index (and therefore precomputed GTF annotation, tx2id etc.)
        if tx2id == "tx2le":
            return config["precomputed_tx2le"]
        
        elif tx2id == "tx2gene":
            return config["precomputed_tx2gene"]
        
        elif tx2id == "le2gene":
            return config["precomputed_le2gene"]
        
        else:
            return config["precomputed_info_tbl"]

    elif run_identification:
        # use output of get_combined_quant_gtf
        if tx2id == "tx2le":
            return rules.get_combined_quant_gtf.output.tx2le
        
        elif tx2id == "tx2gene":
            return rules.get_combined_quant_gtf.output.tx2gene
        
        elif tx2id == "le2gene":
            return rules.get_combined_quant_gtf.output.le2gene
        
        else:
            return rules.get_combined_quant_gtf.output.info_tbl
    
    else:
        # use output of get_ref_quant_gtf (i.e. reference/input GTF only quantified)
        if tx2id == "tx2le":
            return rules.get_ref_quant_gtf.output.tx2le
        
        elif tx2id == "tx2gene":
            return rules.get_ref_quant_gtf.output.tx2gene
        
        elif tx2id == "le2gene":
            return rules.get_ref_quant_gtf.output.le2gene
        
        else:
            return rules.get_ref_quant_gtf.output.info_tbl


rule tx_to_le_quant:
    input:
        tx2le = id2id_target(config["run_identification"], config["use_precomputed_salmon_index"], "tx2le"),
        tx2gene = id2id_target(config["run_identification"], config["use_precomputed_salmon_index"], "tx2gene"),
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


rule make_formulas_txt:
    output:
        os.path.join(DAPA_SUBDIR, "formulas.txt")
    
    params:
        full = config["dexseq_formula_full"],
        reduced = config["dexseq_formula_reduced"]
    
    run:
        with open(output[0], "w") as outfile:
            outfile.write(params.full + "\n")
            outfile.write(params.reduced + "\n")


rule dexseq_apa:
    input:
        counts = rules.tx_to_le_quant.output.counts,
        le2gene =  id2id_target(config["run_identification"], config["use_precomputed_salmon_index"], "le2gene"),
        sample_tbl = config["sample_tbl"],
        formulas = rules.make_formulas_txt.output

    output:
        os.path.join(DAPA_SUBDIR, "dexseq_apa.results.tsv")

    params:
        script = os.path.join("scripts", "run_dexseq.R"),
        output_prefix = os.path.join(DAPA_SUBDIR, "dexseq_apa"),
        min_mean_count = config["min_mean_count"],
        min_rel_usage = config["min_relative_usage"],
        contrast_name = CONTRAST_NAME,
        base_key = BASE_KEY,
        contrast_key = CONTRAST_KEY,
        condition_col = "condition" # set by pipeline

    log:
        stdout = os.path.join(LOG_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "dexseq_apa.stdout.log"),
        stderr = os.path.join(LOG_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "dexseq_apa.stderr.log"),

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "dexseq_apa.txt")

    threads:
        config["dexseq_threads"]

    conda:
        "../envs/papa_r.yaml"

    shell:
        """
        Rscript {params.script} \
        -i {input.counts} \
        -g {input.le2gene} \
        -s {input.sample_tbl} \
        --formulas {input.formulas} \
        --base-key {params.base_key} \
        --contrast-key {params.contrast_key} \
        -n {params.contrast_name} \
        --condition-col {params.condition_col} \
        -m {params.min_mean_count} \
        -r {params.min_rel_usage} \
        -c {threads} \
        -o {params.output_prefix} \
        > {log.stdout} \
        2> {log.stderr}
        """


rule process_dexseq_tbl:
    input:
        dexseq_tbl = rules.dexseq_apa.output,
        info_tbl = id2id_target(config["run_identification"], config["use_precomputed_salmon_index"], "info_tbl"),
        ppau = rules.tx_to_le_quant.output.ppau

    output:
        os.path.join(DAPA_SUBDIR, "dexseq_apa.results.processed.tsv")

    params:
        script = os.path.join("scripts", "process_dexseq_tbl.R"),
        output_prefix = os.path.join(DAPA_SUBDIR, "dexseq_apa.results")

    log:
        stdout = os.path.join(LOG_SUBDIR, config["diff_apa_subdir_name"], "process_dexseq_tbl.stdout.log"),
        stderr = os.path.join(LOG_SUBDIR, config["diff_apa_subdir_name"], "process_dexseq_tbl.stderr.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "process_dexseq_tbl.txt")

    conda:
        "../envs/papa_r.yaml"

    shell:
        """
        Rscript {params.script} \
        -i {input.dexseq_tbl} \
        -a {input.info_tbl} \
        -p {input.ppau} \
        -o {params.output_prefix} \
        > {log.stdout} \
        2> {log.stderr}
        """


rule saturn_apa:
    input:
        counts = rules.tx_to_le_quant.output.counts,
        le2gene =  id2id_target(config["run_identification"], config["use_precomputed_salmon_index"], "le2gene"),
        sample_tbl = config["sample_tbl"]

    output:
        saturn_tbl = os.path.join(DAPA_SUBDIR,
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


rule process_saturn_tbl:
    input:
        saturn_tbl = rules.saturn_apa.output.saturn_tbl,
        info_tbl = id2id_target(config["run_identification"], config["use_precomputed_salmon_index"], "info_tbl"),
        ppau = rules.tx_to_le_quant.output.ppau
    
    output:
        processed_tbl = os.path.join(DAPA_SUBDIR,
                                    "saturn_apa.results.processed.tsv")
    
    params:
        script = "scripts/process_saturn_tbl.R",
        output_prefix = os.path.join(DAPA_SUBDIR,
                                    "saturn_apa.results")

    conda:
        "../envs/papa_r.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "process_saturn_tbl.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["diff_apa_subdir_name"],
                     "process_saturn_tbl.txt")

    shell:
        """
        Rscript {params.script} \
        -i {input.saturn_tbl} \
        -a {input.info_tbl} \
        -p {input.ppau} \
        -o {params.output_prefix} \
        &> {log}
        """   