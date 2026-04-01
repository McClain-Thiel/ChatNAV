process NETMHCIIPAN {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/netmhc:0.1.0'

    input:
    tuple val(patient_id), path(candidates_fasta), path(hla_alleles)

    output:
    tuple val(patient_id), path("netmhciipan_output.txt"), emit: predictions

    script:
    """
    # Run NetMHCIIpan 4.3 for Class II alleles (DRB1, DQB1)
    > netmhciipan_output.txt

    while IFS= read -r allele; do
        # Only Class II alleles
        case "\${allele}" in
            HLA-DR*|HLA-DQ*|HLA-DP*)
                netmhciipan -f ${candidates_fasta} \\
                    -a "\${allele}" \\
                    -length 13,15,17,19,21 \\
                    >> netmhciipan_output.txt 2>/dev/null || true
                ;;
        esac
    done < ${hla_alleles}
    """
}
