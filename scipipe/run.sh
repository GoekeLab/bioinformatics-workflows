#!/bin/bash
go run example.go \
    -ref ../test_data/transcriptome.fa \
    -left ../test_data/reads_1.fq.gz \
    -right ../test_data/reads_2.fq.gz \
    -outdir results
