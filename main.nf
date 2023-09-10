process rna_simpleaf_index_expanded {
    input:
        path genome
        path gtf
        
    output:
        path "simpleaf_index", emit: simpleaf_index
    script:
        cmd = "simpleaf index --output simpleaf_index --threads ${params.num_threads} --kmer-length ${params.kmer_length}"
        // options
        if (params.keep_duplicates) {
            cmd += " --keep-duplicates"
        }
        if (params.sparse)  {
            cmd += " --sparse"
        }
        // piscem options
        if (params.use_piscem) {
            cmd += " --use-piscem --minimizer-length ${params.minimizer_length}"
        }

        // expanded reference options
        cmd += " --fasta ${genome} --gtf ${gtf} --ref-type ${params.ref_type} --rlen ${params.read_length}"

        if (params.dedup) {
            cmd += " --dedup"
        }

        if (params.spliced != null) {
            cmd += " --spliced ${params.spliced}"
        }
        if (params.unspliced != null) {
            cmd += " --unspliced ${params.unspliced}"
        }

        """
        ${cmd}
        """
}


process make_template {
    input:
        // rna related 
        val simpleaf_index_path
        path rna_read1
        path rna_read2
        // adt related
        val adt_reference_barcode_csv_gz
        path adt_read1
        path adt_read2
        // hto related
        val hto_reference_barcode_csv_gz
        path hto_read1
        path hto_read2
        // meta parameters
        val sample_name
        val threads
        val use_piscem

    script:
        template 'cite-seq.sh'
}