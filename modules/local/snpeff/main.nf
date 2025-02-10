process SNPEFF {
    label 'process_medium'
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.merged.annotated.vcf"),            emit: snpeff_vcf
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        snpEff -Xmx${avail_mem}M $args ${vcf} > ${prefix}.merged.annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.merged.annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """ 
}


process SNPEFF_SV_ANN {
    label 'process_medium'
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.BND.annotated.vcf"),               emit: snpeff_BND
        tuple val(group), val(meta), file("*.TANDEM.SV_annotated.vcf"),         emit: snpeff_TANDEM
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        snpEff -Xmx${avail_mem}M $args ${vcf} > ${prefix}.SV.annotated.vcf
        grep -e '^#' -e 'MantaBND:' ${prefix}.SV.annotated.vcf >  ${prefix}.BND.annotated.vcf
        grep -v 'MantaBND:' ${prefix}.SV.annotated.vcf  | grep -e '^#' -e 'MantaINV:' >  ${prefix}.INV.annotated.vcf
        grep -v 'MantaBND:' ${prefix}.SV.annotated.vcf | grep -v 'MantaINV:' | grep -e '^#' -e 'TANDEM' >  ${prefix}.TANDEM.SV_annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        """
        echo ${args}
        touch ${prefix}.BND.annotated.vcf
        touch ${prefix}.TANDEM.SV_annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """ 
}

