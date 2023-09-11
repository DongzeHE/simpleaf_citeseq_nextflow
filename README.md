# A nextflow workflow for processing CITE-seq data using simpleaf 

This workflow is designed to process multiple CITE-seq datasets from the same species in one run using the [simpleaf](https://simpleaf.readthedocs.io/en/latest/) workflow program. The reason we need this nextflow program is because currently simpleaf workflow processes one sample at a time. The simpleaf team is working on providing this functionality. When that is done, this nextflow program will be obsolete.

## Usage
Running the nextflow program is easy: you need to have a nextflow installed on your system by following their [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation). Once you have nextflow installed, you can fill in the parameters in the `nextflow.config` file by replacing the existing values with your own. For example, you can use an existing salmon index instead of building a new one by replacing the `salmon_index` parameter with the path to your salmon index. You can also select to provide a transcriptome fasta file to build a salmon/piscem index directly from it, or providing a genome FASTA file and gene annotaiton GTF file to build a expanded transcriptome reference, for example "spliced+unspliced" or "spliced+intronic". **NOTICE** that you should only choose one set of parameters to build the index.

```
params {
    rna_simpleaf_index {
        existing_index_dir {
            index_dir = "/path/to/salmon/index"
            t2g_path = "/path/to/t2g.tsv"
        }
    }
}
```

This program takes 4 parts of input:
1. A csv file recording the information of the samples you want to process.
2. A salmon/piscem index folder containing the index files, or the files used for generating a new index.
3. The reference barcode CSV file of HTO and ADT barcodes.
4. The parameters used for running simpleaf.

We will talk about each part in detail in the later sections.

Once you have the parameters filled in, you can run the nextflow program by typing the following command in your terminal:
```bash
nextflow main.nf -resume
```

Here we use the `-resume` flag because we want to resume the execution if possible, from where it left off in case it was interrupted. If you want to start the workflow from the beginning, you can omit the `-resume` flag.

# Input files
## Sample information TSV file
The sample information TSV file is a TSV file that contains the information of the samples you want to process. The file should have the following columns:
1. `sample_name`: The name of the sample. This will be used as the prefix of the output files.
3. `RNA_read1`: The path to the read1 FASTQ files of the RNA-seq data, separated by comma if there are multiple files.
4. `RNA_read2`: The path to the read2 FASTQ files of the RNA-seq data, separated by comma if there are multiple files.
5. `HTO_read1`: The path to the read1 FASTQ files of the HTO data, separated by comma if there are multiple files.
6. `HTO_read2`: The path to the read2 FASTQ files of the HTO data, separated by comma if there are multiple files.
7. `ADT_read1`: The path to the read1 FASTQ files of the ADT data, separated by comma if there are multiple files.
8. `ADT_read2`: The path to the read2 FASTQ files of the ADT data, separated by comma if there are multiple files.

```tsv
sample_name     rna_read1       rna_read2       adt_read1       adt_read2       hto_read1       hto_read2
sample1 rna_read1_1.fastq,rna_read1_2.fastq     rna_read2_1.fastq,rna_read2_2.fastq     adt_read1.fastq adt_read2.fastq hto_read1.fastq hto_read2.fastq
sample2 rna_read1_1.fastq,rna_read1_2.fastq     rna_read2_1.fastq,rna_read2_2.fastq     adt_read1.fastq adt_read2.fastq hto_read1.fastq hto_read2.fastq

```

**Please** make sure that there is an empty line at the bottom of the TSV file. Otherwise, nextflow will ignore the last sample. 




