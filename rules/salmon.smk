

rule custom_txome_fasta:
    '''
    Generate FASTA file of reference and novel filtered transcripts
    For use with Salmon
    '''
    input:
        os.path.join(STRINGTIE_SUBDIR, "all_samples.intron_chain_filtered.ref_merged.gtf")

    output:
        os.path.join(SALMON_SUBDIR, "papa.transcripts.fa")

    params:
        genome_fa = config["genome_fasta"]

    log:
        os.path.join(LOG_SUBDIR, "custom_txome_fasta.gffread.log")

    conda:
        "../envs/papa.yaml"

    shell:
        """
        gffread \
        -w {output} \
        -g {params.genome_fa} \
        {input}
        """

rule generate_full_decoys:
    '''
    Generate combined FASTA of target transcripts and rest of genome
    Used to generate selective-alignment compatible FASTA file for Salmon index
    https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    '''
    input:
        genome_fa = config["genome_fasta"],
        txome_fa = os.path.join(SALMON_SUBDIR, "papa.transcripts.fa"),
    output:
        gentrome_fa = os.path.join(SALMON_SUBDIR, "gentrome.fa"),
        decoys = os.path.join(SALMON_SUBDIR, "decoys.txt")

    log:
        os.path.join(LOG_SUBDIR, "generate_full_decoys.log")

    shell:
        """
        grep "^>" {input.genome_fa} | cut -d " " -f 1 > {output.decoys} && \
        sed -i.bak -e 's/>//g' {output.decoys} && \
        cat {input.txome_fa} {input.genome_fa} > {output.gentrome_fa} \
        2> {log}
        """

rule salmon_index:
    input:
        gentrome_fa = os.path.join(SALMON_SUBDIR, "gentrome.fa"),
        decoys = os.path.join(SALMON_SUBDIR, "decoys.txt")

    output:
        os.path.join(SALMON_SUBDIR, "index", "seq.bin"),
        os.path.join(SALMON_SUBDIR, "index", "pos.bin")

    params:
        k = config["kmer_size"],
        outdir = os.path.join(SALMON_SUBDIR, "index","")

    threads:
        config["salmon_index_threads"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR, "salmon_index.log")

    shell:
        """
        salmon index \
        -t {input.gentrome_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        -p {threads} \
        2> {log}
        """


rule salmon_quant_pe:
    input:
        fast1 = lambda wildcards: OPTIONS[wildcards.sample]["fastq1"],
        fast2 = lambda wildcards: OPTIONS[wildcards.sample]["fastq2"],
        index = os.path.join(SALMON_SUBDIR, "index", "seq.bin")

    output:
        os.path.join(SALMON_SUBDIR, "quant", "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(SALMON_SUBDIR, "index",),
        output_dir = os.path.join(SALMON_SUBDIR, "quant", "{sample}"),
        libtype = "A"

    threads:
        config["salmon_quant_threads"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR, "salmon_quant_pe.log")

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --threads {threads} \
        -o {params.output_dir} \
        2> {log}
        """
