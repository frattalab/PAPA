#!/usr/bin/env bash

snakemake -p \
	--configfile="config/test_data_config.yaml" \
	--use-conda \
	--cores 1
