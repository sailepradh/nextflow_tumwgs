#!/usr/bin/env nextflow

/*
* 'TumWgs' - A tumor whole genome sequencing pipeline that takes use senteion genomic distributed workflow  mode.   
*/ 


nextflow.enable.dsl = 2

log.info """\
======================================================================
TUMWGS -NF v3.0
======================================================================
outdir                  :       $params.outdir
subdir                  :       $params.subdir
crondir                 :       $params.crondir
genome                  :       $params.genome_file
csv                     :       $params.csv
=====================================================================
"""

include { SWGP_COMMON               } from './workflows/common.nf'

println(params.genome_file)
csv = file(params.csv)

/*
    Different workflows for the analysis; current entry point is SPP i.e -entry SPP
*/

workflow SWGP {
    SWGP_COMMON ()
}
