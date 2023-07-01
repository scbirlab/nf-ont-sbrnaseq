# Nanopore sbRNA-seq pipeline

Nextflow pipeline to process Nanopore POD5 files from multiple bacterial samples into a gene $\times$ cell count table.

Duplex basecalling is carried out, which requires v10.4 chemistry to have been used. 

## Processing steps

1. Merge all POD5 files (in case they've been already demultiplexed) into a single file.
2. Basecall the merged POD5 file using `guppy` duplex basecalling, producing a `.fastq.gz`. 
3. Demultiplex the `.fastq.gz` based on the `sample_sheet` barcodes using `guppy` according to the parameter `barcode_kit`.

Then for each demultiplexed sample:

1. Trim reads to adapters using `cutadapt`. Reads from sbRNA-seq will be flanked by sequences containing cell barcodes and UMIs, so this step trims any extra sequences either side of these flanking sequences.
2. Extract cell barcodes and UMIs using `umitools`.
3. Map to genome FASTA using `minimap2`.
4. Deduplicate mapped reads using `umitools`.
5. Count deduplicated reads per gene using `featureCounts`.
6. Count deduplicated reads per gene per cell using `umitools`.

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Map to genome FASTA using `bowtie2` because `minimap2` logs are not compatible with `multiqc`. This way, some kind of alignment metrics are possible.
3. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

You need to have Nextflow and either `conda` or `mamba` installed on your system. If possible, use `mamba` because it will be faster.

You will also need the `guppy` basecaller from Oxford Nanopore. It can be downloaded from their [community site](https://community.nanoporetech.com/downloads). When you've installed it, `guppy_path` is a required parameter of the pipeline.

### Reference genome and genome annotations

You also need the genome FASTA and GFF annotations for the bacteria you are sequencing. These can be obtained from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/):

1. Search for your strain of interest, and open its main page
2. On the right-hand side, click `Customize view`, then `Customize` and check `Show sequence`. Finally, click `Update view`. You may have to wait a few minute while the sequence downloads.
3. Click `Send to: > Complete record > File > [FASTA or GFF3] > Create file`
4. Save the files to directories which you provide as parameters below.

### First time using Nextflow?

If it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet) and, optionally, a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-ont-sbrnaseq
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which will not be automatically updated. If you want to ensure that you're using the very latest version of the pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-ont-sbrnaseq -latest
```
If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so using

```bash 
nextflow run scbirlab/nf-ont-sbrnaseq -r v0.0.1
```

For help, use `nextflow run scbirlab/nf-ont-sbrnaseq --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed
- `data_dir`: Path to root directory to find POD5 files (using pattern `<data_dir>/pod5_*/*/*.pod5`)
- `genome_fasta_dir`: path to directory containing genome FASTA files (for mapping)
- `genome_gff_dir`: path to directory containing genome GFF files (for feature counting)
- `guppy_path`: path to Guppy (`ont-guppy` directory) provided by Oxford Nanopore
- `barcode_kit`: SKU of the barcoding kit used, e.g. `SQK-NBD114-24` for the 24-barcode ligation kit 

The following parameters have default values can be overridden if necessary.

- `model = "dna_r10.4.1_e8.2_400bps"` : For `guppy`, the basecalling model to use.
- `trim_qual = 10` : For `cutadapt`, the minimum Phred score for trimming 3' calls
- `min_length = 11` : For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
- `umitools_error = 6`: For `umitools`, the number of errors allowed to correct cell barcodes
- `strand = 1` : For `featureCounts`, the strandedness of RNA-seq. `1` for forward, `2` for reverse.
- `ann_type = 'gene'` : For `featureCounts`, features from GFF column 3 to use for counting
- `label = 'Name'` : For `featureCounts`, one or more (comma-separated) fields from column 9 of GFF for labeling counts

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
   
    sample_sheet = "/path/to/sample_sheet.csv"
    data_dir = "/path/to/pod5/root"

    guppy_path = "/path/to/ont-guppy"

    genome_fasta_dir = "/path/to/fastas"
    genome_gff_dir = "/path/to/gffs"

    barcode_kit = "SQK-NBD114-24"
}
```

Alternatively, you can provide these on the command line:

```bash
nextflow run scbirlab/nf-ont-sbrnaseq \
    --sample_sheet /path/to/sample_sheet.csv \
    --data_dir /path/to/pod5/root \
    --guppy_path /path/to/ont-guppy \
    --genome_fasta_dir /path/to/fastas \
    --genome_gff_dir /path/to/gffs \
    --barcode_kit SQK-NBD114-24
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which demultiplexed FASTQ files belong to which sample, which genome each sample should be mapped to, and the UMI and cell barcode scheme for each sample.

The file must have a header with the column names below, and one line per sample to be processed.

- `barcode_id`: the ONT barcode name from the barcoding kit you used
- `sample_id`: the unique name of the sample
- `genome_id`: the name of the genome to map to. Each entry must match the name of one file (apart from the extension) in `genome_fasta_dir` and `genome_gff_dir`
- `n_cells`: maximum number of uniquely barcoded cells in the sample (used by `umi_tools`). This can be more than the number of reads, in which case the number of reads is taken instead.
- `adapter`: the adapter sequence to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequences _either side_ will be removed, but the adapters themselves will be retained.
- `umi`: the cell barcode and UMI pattern in [`umitools` regex format](https://umi-tools.readthedocs.io/en/latest/regex.html#regex-regular-expression-mode)

Here is an example of the sample sheet:

| barcode_id | sample_id | genome_id             | n_cells | adapter   | umi |
| ---------- | --------- | --------------------- | ------- | --------- | --- |
| barcode01  | Eco1      | EcoMG1655-NC_000913.3 | 885000 | AGACAGN{6}G{3}...N{7}AGATCG | ^(?P<discard_1>.{6})(?P<cell_1>.{6})(?P<discard_2>.{3}).+(?P<cell_3>.{7})(?P<discard_4>.{6})$ |
| barcode02  | Eco2      | EcoMG1655-NC_000913.3 | 885000 | AGACAGN{6}G{3}...N{7}AGATCG | ^(?P<discard_1>.{6})(?P<cell_1>.{6})(?P<discard_2>.{3}).+(?P<cell_3>.{7})(?P<discard_4>.{6})$ |

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under four directories:

- `fastq`: Demultiplexed FASTQ files.
- `processed`: FASTQ files and logs resulting from trimming and UMI extraction
- `counts`: tables and BAM files corresponding to cell $\times$ gene counts
- `multiqc`: HTML report on processing steps

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-ont-sbrnaseq/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [umitools](https://umi-tools.readthedocs.io/en/latest/index.html)
- [featureCounts](https://subread.sourceforge.net/featureCounts.html)
- [minimap2](https://lh3.github.io/minimap2/minimap2.html)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools](http://www.htslib.org/doc/samtools.html)