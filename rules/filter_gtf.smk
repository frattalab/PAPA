#     Snakemake rule to filter a reference GTF based on reference flags
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
