# A nextflow workflow for processing CITE-seq data using simpleaf 

This workflow is designed to process multiple CITE-seq datasets from the same species in one run using the [simpleaf](https://simpleaf.readthedocs.io/en/latest/) workflow program. The reason we need this nextflow program is because currently simpleaf workflow processes one sample at a time. The simpleaf team is working on providing this functionality. When that is done, this nextflow program will be obsolete.

## Usage
Running nextflow program is easy: you need to have a nextflow installed on your system by following their [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation). Once you have nextflow installed, you can fill in the parameters in the `nextflow.config` file by replacing the existing values with your own. For example, you can use an existing salmon index instead of building a new one by replacing the `salmon_index` parameter with the path to your salmon index.
```
params {
    simpleaf_index {
        existing_index_dir = "/path/to/salmon/index"
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

Here we use the `-resume` flag because we want to resume the workflow if possible from where it left off in case it was interrupted. If you want to start the workflow from the beginning, you can omit the `-resume` flag.

# Input files
## Sample information CSV file
The sample information CSV file is a CSV file that contains the information of the samples you want to process. The file should have the following columns:
1. `sample_name`: The name of the sample. This will be used as the prefix of the output files.
2. `chemistry`: The chemistry used for the sample. Currently only `10xv2` and `10xv3` are supported.
3. `RNA_read1`: The path to the read1 FASTQ files of the RNA-seq data, separated by comma if there are multiple files.
4. `RNA_read2`: The path to the read2 FASTQ files of the RNA-seq data, separated by comma if there are multiple files.
5. `HTO_read1`: The path to the read1 FASTQ files of the HTO data, separated by comma if there are multiple files.
6. `HTO_read2`: The path to the read2 FASTQ files of the HTO data, separated by comma if there are multiple files.
7. `ADT_read1`: The path to the read1 FASTQ files of the ADT data, separated by comma if there are multiple files.
8. `ADT_read2`: The path to the read2 FASTQ files of the ADT data, separated by comma if there are multiple files.







