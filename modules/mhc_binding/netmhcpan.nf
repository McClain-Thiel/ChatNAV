process NETMHCPAN {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/netmhc:0.1.0'

    input:
    tuple val(patient_id), path(candidates_fasta), path(hla_alleles)

    output:
    tuple val(patient_id), path("netmhcpan_output.txt"), emit: predictions

    script:
    """
    # Run NetMHCpan 4.1 for each HLA-A/B/C allele
    > netmhcpan_output.txt

    while IFS= read -r allele; do
        # Skip Class II alleles
        case "\${allele}" in
            HLA-DR*|HLA-DQ*|HLA-DP*) continue ;;
        esac

        # NetMHCpan expects allele format: HLA-A02:01
        netmhcpan -p ${candidates_fasta} \\
            -a "\${allele}" \\
            -l 8,9,10,11 \\
            -BA \\
            >> netmhcpan_output.txt 2>/dev/null || true
    done < ${hla_alleles}
    """
}
