#!/usr/bin/env bash

snakemake -p \
	--configfile $1 \
	--use-conda \
	--cores 2
