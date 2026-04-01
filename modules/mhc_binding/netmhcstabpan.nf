process NETMHCSTABPAN {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/netmhc:0.1.0'

    input:
    tuple val(patient_id), path(candidates_fasta), path(hla_alleles)

    output:
    tuple val(patient_id), path("netmhcstabpan_output.txt"), emit: predictions

    script:
    """
    > netmhcstabpan_output.txt

    while IFS= read -r allele; do
        case "\${allele}" in
            HLA-DR*|HLA-DQ*|HLA-DP*) continue ;;
        esac

        netMHCstabpan -p ${candidates_fasta} \\
            -a "\${allele}" \\
            -l 8,9,10,11 \\
            >> netmhcstabpan_output.txt 2>/dev/null || true
    done < ${hla_alleles}
    """
}
