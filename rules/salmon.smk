#     Snakemake rules to perform Salmon index generation and quantification
#     Copyright (C) 2024  Sam Bryce-Smith samuel.bryce-smith.19@ucl.ac.uk

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

if single_end:
    ruleorder: salmon_quant_se > salmon_quant_pe

else:
    ruleorder: salmon_quant_pe > salmon_quant_se


rule custom_txome_fasta:
    '''
    Generate FASTA file of reference and novel filtered transcripts
    For use with Salmon
    '''
    input:
        rules.get_combined_quant_gtf.output.quant_gtf if config["run_identification"] else rules.get_ref_quant_gtf.output.quant_gtf

    output:
        os.path.join(SALMON_SUBDIR, "papa.transcripts.fa")

    params:
        genome_fa = config["genome_fasta"]

    log:
        os.path.join(LOG_SUBDIR,
                     config["salmon_subdir_name"],
                     "custom_txome_fasta.log")

    conda:
        "../envs/papa.yaml"

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["salmon_subdir_name"],
                     "custom_txome_fasta.txt")

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
        txome_fa = rules.custom_txome_fasta.output
        # os.path.join(SALMON_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}","papa.transcripts.fa"),
    output:
        gentrome_fa = os.path.join(SALMON_SUBDIR, "gentrome.fa"),
        decoys = os.path.join(SALMON_SUBDIR, "decoys.txt")

    log:
        os.path.join(LOG_SUBDIR,
                     config["salmon_subdir_name"],
                     "generate_full_decoys.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["salmon_subdir_name"],
                     "generate_full_decoys.txt")

    shell:
        """
        grep "^>" {input.genome_fa} | cut -d " " -f 1 > {output.decoys} && \
        sed -i.bak -e 's/>//g' {output.decoys} && \
        cat {input.txome_fa} {input.genome_fa} > {output.gentrome_fa} \
        2> {log}
        """


rule salmon_index:
    input:
        gentrome_fa = rules.generate_full_decoys.output.gentrome_fa,
        decoys = rules.generate_full_decoys.output.decoys

    output:
        seq = os.path.join(SALMON_SUBDIR, "index", "seq.bin"),
        pos = os.path.join(SALMON_SUBDIR, "index", "pos.bin")

    params:
        k = config["salmon_kmer_size"],
        outdir = os.path.join(SALMON_SUBDIR, "index", "")

    threads:
        config["salmon_index_threads"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     config["salmon_subdir_name"],
                     "salmon_index.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["salmon_subdir_name"],
                     "salmon_index.txt")

    shell:
        """
        salmon index \
        -t {input.gentrome_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        -p {threads} \
        &> {log}
        """


rule salmon_quant_pe:
    input:
        fast1 = lambda wildcards: OPTIONS[wildcards.sample]["fastq1"],
        fast2 = lambda wildcards: OPTIONS[wildcards.sample]["fastq2"],
        index = rules.salmon_index.output.seq if not config["use_precomputed_salmon_index"] else os.path.join(config["precomputed_salmon_index"], "seq.bin")

    output:
        os.path.join(SALMON_SUBDIR, "quant", "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(SALMON_SUBDIR, "index",) if not config["use_precomputed_salmon_index"] else config["precomputed_salmon_index"],
        output_dir = os.path.join(SALMON_SUBDIR, "quant", "{sample}"),
        extra_params = " ".join(config["salmon_quant_extra_params"]),
        libtype = "A"

    threads:
        config["salmon_quant_threads"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     config["salmon_subdir_name"],
                     "salmon_quant_pe.{sample}.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["salmon_subdir_name"],
                     "salmon_quant_pe.{sample}.txt")

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --threads {threads} \
        -o {params.output_dir} \
        {params.extra_params} \
        &> {log}
        """


rule salmon_quant_se:
    input:
        fast1 = lambda wildcards: OPTIONS[wildcards.sample]["fastq1"],
        index = rules.salmon_index.output.seq if not config["use_precomputed_salmon_index"] else os.path.join(config["precomputed_salmon_index"], "seq.bin")

    output:
        os.path.join(SALMON_SUBDIR, "quant", "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(SALMON_SUBDIR, "index",) if not config["use_precomputed_salmon_index"] else config["precomputed_salmon_index"],
        output_dir = os.path.join(SALMON_SUBDIR, "quant", "{sample}"),
        extra_params = " ".join(config["salmon_quant_extra_params"]),
        libtype = "A"

    threads:
        config["salmon_quant_threads"]

    conda:
        "../envs/papa.yaml"

    log:
        os.path.join(LOG_SUBDIR,
                     config["salmon_subdir_name"],
                     "salmon_quant_pe.{sample}.log")

    benchmark:
        os.path.join(BMARK_SUBDIR,
                     config["salmon_subdir_name"],
                     "salmon_quant_pe.{sample}.txt")

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        -r {input.fast1} \
        --threads {threads} \
        -o {params.output_dir} \
        {params.extra_params} \
        &> {log}
        """
