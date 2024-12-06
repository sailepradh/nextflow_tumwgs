#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { VARDICT                  } from '../../modules/local/vardict/main'
include { TNSCOPE                  } from '../../modules/local/sentieon/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS           } from '../../modules/local/concatenate_vcfs/main'


workflow SNV_CALLING {
    take: 
        bam_bqsr         // channel: [mandatory] [ val(group), val(meta), file("umi.bam"), file("umi.bam.bai"), file(bqsr) ]
        cram_bqsr       // channel: [mandatory] [ val(group), val(meta), file("dedup.bam"), file("dedup.bam.bai") ]
        beds            // channel: [mandatory] [ file(bed) ]
        meta            // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
    main:
        ch_versions = Channel.empty()

        // Variantcallers //
        // split by bed-file to speed up calling //
        FREEBAYES ( bam_bqsr, beds)
        ch_versions         = ch_versions.mix(FREEBAYES.out.versions.first())

        VARDICT ( bam_bqsr, beds)
        ch_versions         = ch_versions.mix(VARDICT.out.versions.first())

       // TNSCOPE ( cram_bqsr, beds)
       // ch_versions         = ch_versions.mix(TNSCOPE.out.versions.first())

        // Prepare vcf parts for concatenation //
        vcfparts_freebayes  = FREEBAYES.out.vcfparts_freebayes.groupTuple(by:[0,1])
        vcfparts_vardict    = VARDICT.out.vcfparts_vardict.groupTuple(by:[0,1])
        // vcfparts_tnscope    = TNSCOPE.out.vcfparts_tnscope.groupTuple(by:[0,1])

        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        // vcfs_to_concat      = vcfparts_freebayes.mix(vcfparts_vardict,vcfparts_tnscope)
        vcfs_to_concat      = vcfparts_freebayes.mix(vcfparts_vardict)

        // Join vcfs split by bedparts //
        CONCATENATE_VCFS { vcfs_to_concat }
        ch_versions         = ch_versions.mix(CONCATENATE_VCFS.out.versions.first())

        // Aggregate all callers to one VCF
        AGGREGATE_VCFS { CONCATENATE_VCFS.out.concatenated_vcfs.groupTuple().join(meta.groupTuple())  }
        AGGREGATE_VCFS.out.vcf_concat.view()
        ch_versions         = ch_versions.mix(AGGREGATE_VCFS.out.versions.first())

    emit:
        concat_vcfs =   CONCATENATE_VCFS.out.concatenated_vcfs  // channel: [ val(group), val(vc), file(vcf.gz) ]
        agg_vcf     =   AGGREGATE_VCFS.out.vcf_concat           // channel: [ val(group), val(meta), file(agg.vcf) ]
        versions    =   ch_versions                             // channel: [ file(versions) ]

}