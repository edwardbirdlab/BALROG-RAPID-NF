process AMRFINDER {

    label 'amrfinder'
    container 'ncbi/amr:latest'

    input:
        tuple val(sample), path(contigs)

    output:
        tuple val(sample), path("${sample}_amrfinder.tsv"), emit: results
        path("versions.yml"),                                emit: versions

    script:
    """
    if [ -s ${contigs} ]; then
        amrfinder \
            -n ${contigs} \
            -o ${sample}_amrfinder.tsv \
            --threads ${task.cpus} \
            --plus
    else
        # Empty contigs - create header-only AMRFinder output
        echo -e "Protein id\\tContig id\\tStart\\tStop\\tStrand\\tElement symbol\\tElement name\\tScope\\tType\\tSubtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference\\t% Identity to reference\\tAlignment length\\tClosest reference accession\\tClosest reference name\\tHMM id\\tHMM description" > ${sample}_amrfinder.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinder: \$(amrfinder -V 2>&1 | grep "Software version" | sed -e "s/Software version: //g" || echo "unknown")
    END_VERSIONS
    """
}
