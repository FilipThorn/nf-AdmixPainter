#!/bin/bash -l

module load bioinfo-tools Nextflow/22.10.1 ANGSD/0.940-stable gcc/10.3.0 cairo/1.17.4 texinfo/6.8 texlive/2022-09-18 libcurl/7.85.0 readline/6.2-11 libicu/5.2-4 xz/5.2.6 bzip2/1.0.8 zlib/1.2.12 R/4.2.1 java/OracleJDK_11.0.9 R_packages/4.2.1 #environment

NXF_HOME=/PATH/nf-admixPainter

export NXF_DEFAULT_DSL=1 

nextflow run main.nf -profile rackham -resume
