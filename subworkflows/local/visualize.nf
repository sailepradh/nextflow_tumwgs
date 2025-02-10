include { GENS_VIZ              } from '../../modules/local/gens/main'

workflow VISUALIZE {
    take: 
        dnascope_vcf           
        gatk_count    

    main:
        ch_versions = Channel.empty()

        GENS_VIZ ( dnascope_vcf, gatk_count  )
        ch_versions = ch_versions.mix(GENS_VIZ.out.versions)

    emit:
        gens        =   GENS_VIZ.out.gens
        versions    =   ch_versions                     // channel: [ file(versions) ]
}