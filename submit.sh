#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


if [ "$2" != "" ]; then
    RUN_NAME="PAPA"
else
    RUN_NAME=$2
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp $1 ${FOLDER}/${RUN_NAME}_config.yaml

snakemake \
-p \
--use-conda \
--configfile $1 \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 25 \
--nolock \
--rerun-incomplete \
--latency-wait 50
