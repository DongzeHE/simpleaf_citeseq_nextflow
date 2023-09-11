// this process build the expanded index for simpleaf
process rna_simpleaf_index_expanded {
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy'
    input:
        path genome
        path gtf

    output:
        path "simpleaf_index", emit: index_dir
        path "simpleaf_index/index/t2g_3col.tsv", emit: t2g_path
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

// this process build the expanded index for simpleaf
process rna_simpleaf_index_refseq {
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy'
    input:
        path genome
        path gtf

    output:
        path "simpleaf_index", emit: index_dir
        path "simpleaf_index/index/t2g.tsv", emit: t2g_path
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
        cmd += " --ref-seq ${params.ref_seq}"
        
        // execute command

        """
        ${cmd}
        """
}

// this process build the expanded index for simpleaf
process rna_simpleaf_index_existing {
    publishDir "${params.output_dir}", mode: 'copy'
    input:
        path index_dir
        path t2g_path

    output:
        path index_dir, emit: index_dir
        path t2g_path, emit: t2g_path
    script:
        """
        cp $t2g_path $index_dir
        """
}

workflow make_index {
    main:
        ep = params.rna_simpleaf_index.extended_transcriptome
        eid = params.rna_simpleaf_index.existing_index_dir
        if (eid.index_dir != null) {
            simpleaf_index_dir = rna_simpleaf_index_existing(eid.index_dir, eid.t2g_path)
        } else if (params.rna_simpleaf_index.ref_seq != null) {
            ref_seq = Channel.fromPath(params.rna_simpleaf_index.ref_seq, checkIfExists: true)
            simpleaf_index_dir = rna_simpleaf_index_refseq(ref_seq)
        } else if (ep.fasta != null && ep.gtf != null) {
            fasta = Channel.fromPath(ep.fasta, checkIfExists: true)
            gtf = Channel.fromPath(ep.fasta, checkIfExists: true)
            simpleaf_index_dir = rna_simpleaf_index_expanded(fasta, gtf)
        } else {
            error "None of the index options are correctly specified; Cannot proceed"
        }
    emit:
        index_dir = simpleaf_index_dir.index_dir
        t2g_path = simpleaf_index_dir.t2g_path
}


// this process makes the workflow template file
process make_template {
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy'
    input:
        tuple val(sample_name),
            val(rna_read1),
            val(rna_read2),
            val(adt_read1),
            val(adt_read2),
            val(hto_read1),
            val(hto_read2),
            val(num_threads),
            val(use_piscem),
            path(adt_reference_barcode_csv_gz),
            path(hto_reference_barcode_csv_gz)

        path simpleaf_index_path
        path t2g_path

    output:
        val sample_name
        path simpleaf_index_path
        path t2g_path
        path adt_reference_barcode_csv_gz
        path hto_reference_barcode_csv_gz
        path "${sample_name}_instantiated_template.jsonnet"
    script:
        template 'cite-seq.sh'
}

process run_simpleaf_workflow {
    publishDir "${params.output_dir}", mode: 'copy'
    input:
        val sample_name
        path simpleaf_index_path
        path t2g_path
        path adt_reference_barcode_csv_gz
        path hto_reference_barcode_csv_gz
        path template_file
    output:
        path "$sample_name"
    
    script:
        cmd = "simpleaf workflow run --template ${template_file} --output ${sample_name}"

        """
        export ALEVIN_FRY_HOME="af_home"
        simpleaf set-paths

        ${cmd}
        """
    

}

workflow {
    // create output directory
    output_dir = file("${params.output_dir}/log")
    output_dir.mkdirs()
    file("$projectDir/nextflow.config").copyTo(file("${output_dir}/nextflow.config"))

    // make index
    make_index()
    
    // process sample sheet
    sample_sheet = Channel
      .fromPath(params.sample_sheet_tsv)
      .splitCsv(header:true, sep: "\t", strip: true)
      .map{ row-> tuple(row.sample_name,
                        row.rna_read1,
                        row.rna_read2,
                        row.adt_read1,
                        row.adt_read2,
                        row.hto_read1,
                        row.hto_read2)
      }
    
    // parse parameters
    adt_reference_barcode_csv_gz = Channel.fromPath(params.adt_reference_barcode_csv_gz, checkIfExists: true)
    hto_reference_barcode_csv_gz = Channel.fromPath(params.hto_reference_barcode_csv_gz, checkIfExists: true)
    num_threads = Channel.of("$params.num_threads")
    use_piscem = Channel.of("$params.use_piscem")

    // make_template input 
    make_template_input = sample_sheet
        .combine(num_threads)
        .combine(use_piscem)
        .combine(adt_reference_barcode_csv_gz)
        .combine(hto_reference_barcode_csv_gz)

    make_template(make_template_input, make_index.out.index_dir, make_index.out.t2g_path) | run_simpleaf_workflow
}

