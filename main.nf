#!/usr/bin/env nextflow

/*
========================================================================================
   Nanopore RNASeq Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-ont-sbrnaseq
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println """\
         S C B I R   N A N O P O R E  s b R N A - S E Q   P I P E L I N E
         ================================================================
         Nextflow pipeline to process raw Nanopore POD5 files from 
         multiple bacterial samples into a gene x cell count table. Applies 
         duplex basecalling, therefore requires v10.4 chemistry.

         Usage:
            nextflow run sbcirlab/nf-ont-sbrnaseq --sample_sheet <csv> --data_dir <dir> --genome_fasta_dir <dir> --genome_gff_dir <dir> --guppy_path <path>
            nextflow run sbcirlab/nf-ont-sbrnaseq -c <config-file>

         Required parameters:
            sample_sheet         Path to a CSV containing sample IDs matched with ONT barcode IDs, genome information, and adapter sequences.
            data_dir             Path to root directory to find POD5 files (using pattern <data_dir>/pod5_*/*/*.pod5).
            genome_fasta_dir     Path to directory containing genome FASTA files (for mapping)
            genome_gff_dir       Path to directory containing genome GFF files (for feature counting)
            guppy_path           Path to Guppy (`ont-guppy` directory) provided by Oxford Nanopore
            barcode_kit          SKU of the barcoding kit used, e.g. SQK-NBD114-24 for the 24-barcode ligation kit 

         Optional parameters (with defaults):   
            model = "dna_r10.4.1_e8.2_400bps"   For `guppy`, the basecalling model to use.
            trim_qual = 5                       For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = 64                     For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
            umitools_error = 6                  For `umitools`, the number of errors allowed to correct cell barcodes
            strand = 1                          For `featureCounts`, the strandedness of RNA-seq. `1` for forward, `2` for reverse.
            ann_type = 'gene'                   For `featureCounts`, features from GFF column 3 to use for counting
            label = 'Name'                      For `featureCounts`, one or more (comma-separated) fields from column 9 of GFF for labeling counts

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.data_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to fastq_dir")
}
if ( !params.genome_fasta_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to genome_fasta_dir")
}
if ( !params.genome_gff_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to genome_gff_dir")
}
if ( !params.guppy_path ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to guppy_path")
}
if ( !params.barcode_kit ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a barcode_kit")
}

wd = file(params.sample_sheet)
working_dir = wd.getParent()

fastq_o = "${working_dir}/fastq"
processed_o = "${working_dir}/processed"
counts_o = "${working_dir}/counts"
multiqc_o = "${working_dir}/multi_qc"

log.info """\
         S C B I R   N A N O P O R E  s b R N A - S E Q   P I P E L I N E
         ================================================================
         inputs
            sample sheet   : ${params.sample_sheet}
            data directory : ${params.data_dir}
         basecalling
            Guppy path     : ${params.guppy_path}
            model          : ${params.model}
         genome locations
            FASTA          : ${params.genome_fasta_dir}
            GFF            : ${params.genome_gff_dir}
         trimming 
            quality        : ${params.trim_qual}
            minimum length : ${params.min_length}
         UMI detection
            Error number   : ${params.umitools_error}
         output            
            FASTQ          : ${fastq_o}
            Processed      : ${processed_o}
            Counts         : ${counts_o}
            MultiQC        : ${multiqc_o}
         """
         .stripIndent()

dirs_to_make = [fastq_o, processed_o, counts_o, multiqc_o]

log.info  """
         Making directories:
          """.stripIndent()

dirs_to_make.each { 
   log.info "${it}: " 
   log.info file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}

/*
========================================================================================
   Create Channels
========================================================================================
*/

pod5_ch = Channel.fromPath( "${params.data_dir}/pod5_*/*/*.pod5", 
                            checkIfExists: true )
                  .map { tuple( it.getSimpleName(), it ) }

csv_ch = Channel.fromPath( params.sample_sheet, 
                           checkIfExists: true )
                .splitCsv( header: true )

sample_ch = csv_ch.map { tuple( it.barcode_id,
                                it.genome_id,
                                it.sample_id,
                                it.adapter,
                                it.umi, 
                                it.n_cells ) }

genome_ch = csv_ch.map { tuple( it.genome_id,
                                file( "${params.genome_fasta_dir}/${it.genome_id}.fasta",
                                      checkIfExists: true ),
                                file( "${params.genome_gff_dir}/${it.genome_id}.gff3",
                                      checkIfExists: true ) ) }
                  .unique()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   pod5_ch.collect { it[1] } \
      | MERGE_POD5 | BASECALL 
   BASECALL.out.main | DEMUX 
   
   DEMUX.out.main.map { ['_', it] }
                 .transpose()
                 .map { tuple( it[1].getParent().getSimpleName(), 
                               it[1] ) }
                 .join( sample_ch )
                 .set { demuxed }

   demuxed | FASTQC
   demuxed | TRIM_CUTADAPT

   TRIM_CUTADAPT.out.main | UMITOOLS_WHITELIST
   TRIM_CUTADAPT.out.main.join( UMITOOLS_WHITELIST.out.main, 
                                by: 1 ) | UMITOOLS_EXTRACT_WL

   genome_ch | MINIMAP2_INDEX 
   UMITOOLS_EXTRACT_WL.out.main.combine( MINIMAP2_INDEX.out, 
                                         by: 0 ) | MINIMAP2_ALIGN
   MINIMAP2_ALIGN.out.main \
      | UMITOOLS_DEDUP

   genome_ch | BOWTIE2_INDEX 
   UMITOOLS_EXTRACT_WL.out.main.combine( BOWTIE2_INDEX.out, 
                                         by: 0 ) | BOWTIE2_ALIGN

   UMITOOLS_DEDUP.out.main \
      | FEATURECOUNTS

   FEATURECOUNTS.out.main \
      | UMITOOLS_COUNT 

   TRIM_CUTADAPT.out.logs.concat(
         FASTQC.out.logs,
         BOWTIE2_ALIGN.out.logs,
         FEATURECOUNTS.out.logs )
      .flatten()
      .unique()
      .collect() \
      | MULTIQC

}

/*
 * Merge all POD5 files
 */
process MERGE_POD5 {

   time '5h'

   input:
   path( 'input*.pod5' )

   output:
   path( "*.pod5" )

   script:
   """
   pod5 merge input*.pod5 --output ${file(params.data_dir).getBaseName()}.pod5
   """
}

/*
 * Basecall with Guppy (duplex)
 */
process BASECALL {

   memory '128 GB'
   time '24h'
   queue 'gpu'
   clusterOptions '--gres=gpu:4'

   input:
   path( pod5 )

   output:
   path( "*.fastq.gz" ), emit: main
   tuple path( "final/*plex/*.log" ), path( "final/*plex/*.txt" ), path( "final/*plex/*.js" ), emit: logs

   script:
   var guppy = "${params.guppy_path}/bin/guppy_basecaller"
   """
   guppy_duplex --input_path . --save_path . \
      --basecaller_exe ${guppy} \
      --duplex_basecaller_exe "${guppy}"_duplex \
      --call_non_duplex_reads \
      --do_read_splitting \
      --simplex_config ${params.model}_fast.cfg \
      --duplex_config ${params.model}_sup.cfg \
      --duplex_chunks_per_runner 416 \
      --device 'cuda:all:100%'

   cat final/*plex/*/*.fastq | gzip --best > ${file(params.data_dir).getBaseName()}.fastq.gz
   rm final/*plex/*/*.fastq
   """
}

/*
 * Demultiplex basecalled data
 */
process DEMUX {

   memory '128 GB'
   time '24h'
   queue 'gpu'
   clusterOptions '--gres=gpu:4'

   publishDir( fastq_o, 
               mode: 'copy' )

   input:
   path( reads )

   output:
   path( "barcode??/*.fastq.gz" ), emit: main
   tuple path( "*.txt" ),  path( "*.log" ), emit: logs

   script:
   var guppy_barcoder = "${params.guppy_path}/bin/guppy_barcoder"
   """
   ${guppy_barcoder} --input_path . --save_path . \
      --config configuration.cfg \
      --barcode_kits ${params.barcode_kit} \
      --compress_fastq \
      --records_per_fastq 0 \
      --verbose_logs
   """
}

/* 
 * Do quality control checks
 */
process FASTQC {

   memory '32GB'

   tag "${sample_id}"

   input:
   tuple val( barcode_id ), path( reads ), val( genome_id ), val( sample_id ), val( adapters ), val( umis ), val( n_cells )

   output:
   path( "*.zip" ), emit: logs

   script:
   var bc_sample = "${barcode_id}-${sample_id}"
   """
   zcat ${reads} > ${bc_sample}.fastq
   fastqc --noextract --memory 10000 --threads 4 ${bc_sample}.fastq
   rm ${bc_sample}.fastq
   """

}


/*
 * Trim adapters from reads
 * --revcomp does the stranding
 */
process TRIM_CUTADAPT {

   tag "${barcode_id}-${sample_id}"

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( barcode_id ), path( reads ), val( genome_id ), val( sample_id ), val( adapters ), val( umis ), val( n_cells )

   output:
   tuple val( genome_id ), val( "${barcode_id}-${sample_id}" ), val( adapters ), val( umis ), val( n_cells ), path( "*.with-adapters.fastq.gz" ), emit: main
   tuple val( "${barcode_id}-${sample_id}" ), path( "*.without-adapters.fastq.gz" ), emit: no_adapt
   tuple path( "*.log" ),  path( "*.json" ), emit: logs

   script:
   var bc_sample = "${barcode_id}-${sample_id}"
   """
   ln -s ${reads} ${bc_sample}_0.fastq.gz

   cutadapt \
		-g '${adapters}' \
      --overlap 15 \
		-q ${params.trim_qual},${params.trim_qual} \
      --error-rate 0.35 \
      --revcomp \
		--report=full \
      --action=retain \
		--untrimmed-output=${bc_sample}.without-adapters_0.fastq.gz \
		-o ${bc_sample}_1.fastq.gz \
      --json=${bc_sample}_0.cutadapt.json \
		${bc_sample}_0.fastq.gz > ${bc_sample}_0.cutadapt.log

   cutadapt \
		-a '^${adapters}'\$ \
      --error-rate 0.35 \
		--minimum-length ${params.min_length} \
		--report=full \
      --action=retain \
		--untrimmed-output=${bc_sample}.without-adapters_1.fastq.gz \
		-o ${bc_sample}.with-adapters.fastq.gz \
      --json=${bc_sample}.cutadapt.json \
		${bc_sample}_1.fastq.gz > ${bc_sample}.cutadapt.log

   zcat ${bc_sample}.without-adapters_?.fastq.gz | gzip --best > ${bc_sample}.without-adapters.fastq.gz
   rm ${bc_sample}_1.fastq.gz ${bc_sample}.without-adapters_?.fastq.gz
   """
}

/*
 * Identify cell barcodes
 */
process UMITOOLS_WHITELIST {

   tag "${sample_id}"

   time '24h'
   memory '32GB'

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( genome_id ), val( sample_id ), val( adapters ), val( umis ), val( n_cells ), path( reads )

   output:
   tuple path( "*.txt" ), val( sample_id ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   N_LINES=\$(cat "${reads}" | wc -l)
   echo \$N_LINES
   N_READS=\$((\$N_LINES/4))
   N_CELLS=\$(printf \$N_READS'\\n${n_cells}' | sort -g | head -n1)
   echo \$N_LINES \$N_READS \$N_CELLS

   umi_tools whitelist \
      --knee-method=distance \
      --set-cell-number \$N_CELLS \
      --allow-threshold-error \
      --plot-prefix ${sample_id}.whitelist \
		--bc-pattern="${umis}" \
      --extract-method=regex \
      --method=umis \
      --error-correct-threshold=${params.umitools_error} \
      --ed-above-threshold=correct \
		--log ${sample_id}.whitelist.log \
		--stdin ${reads} \
      --stdout ${sample_id}.whitelist.txt
   """
}

/*
 * Extract cell barcodes and UMIs
 */
process UMITOOLS_EXTRACT_WL {

   tag "${sample_id}"

   time '24h'

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), val( genome_id ), val( adapters ), val( umis ), val( n_cells ), path( reads ), path ( whitelist )

   output:
   tuple val( genome_id ), val( sample_id ), path( "*.gz" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   umi_tools extract \
		--whitelist=${whitelist} \
		--bc-pattern="${umis}" \
      --extract-method=regex \
      --error-correct-cell \
      --log ${sample_id}.extract.log \
		--stdin ${reads} \
		--stdout ${sample_id}.extracted.fastq.gz 
   """
}

/*
 * Index the reference genome for use by Bowtie2.
 */
process BOWTIE2_INDEX {

   tag "${genome_id}"
   
   input:
   tuple val( genome_id ), path( fasta ), path( gff )

   output:
   tuple val( genome_id ), path( "*.bt2" ), path( gff )

   script:
   """
   bowtie2-build ${fasta} ${genome_id}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BOWTIE2_ALIGN {

   tag "${sample_id} - ${genome_id}" 

   memory '16GB'
   errorStrategy 'ignore'

   input:
   tuple val( genome_id ), val( sample_id ), path( reads ), path( idx ), path( gff )

   output:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \
      -x ${genome_id} \
      --rdg 10,10 \
      --very-sensitive-local \
      --trim-to 3:2500 \
      -U ${reads} -S ${sample_id}.mapped.sam 2> ${sample_id}.bowtie2.log
   samtools sort ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   samtools index ${sample_id}.mapped.bam
   """
}

/*
 * Index the reference genome for use by Minimap2.
 */
process MINIMAP2_INDEX {

   tag "${genome_id}"
   
   input:
   tuple val( genome_id ), path( fasta ), path( gff )

   output:
   tuple val( genome_id ), path( "*.mmi" ), path( gff )

   script:
   """
   minimap2 -x map-ont -d ${fasta.getBaseName()}-ont.mmi ${fasta}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process MINIMAP2_ALIGN {

   tag "${sample_id}-${genome_id}" 

   publishDir( path: processed_o, 
               mode: 'copy', 
               pattern: "*.mapped.bam" )

   input:
   tuple val( genome_id ), val( sample_id ), path( reads ), path( idx ), path( gff )

   output:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   minimap2 -y --eqx -Y --sam-hit-only --no-end-flt --sr --frag=yes -F 3000 \
      -a ${idx} <(zcat ${reads}) \
      -o ${sample_id}.mapped.sam 2> ${sample_id}.minimap2.log
   samtools sort ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   samtools index ${sample_id}.mapped.bam
   rm ${sample_id}.mapped.sam
   """
}

/*
 * Dedupicate reads based on mapping coordinates and UMIs
 */
process UMITOOLS_DEDUP {

   tag "${sample_id}" 

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.bam" )
   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.log" )

   input:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( bamfile ), path( bam_idx )

   output:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( "*.bam" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   umi_tools dedup \
		--per-cell \
		--stdin ${bamfile} \
		--log ${sample_id}.dedup.log \
      --stdout ${sample_id}.dedup.bam
   """
}


/*
 * Use featureCounts (from the SubReads package) to count transcripts mapping to each gene
 * ## $(ANN_TYPE): which of "CDS", "gene", "mRNA", etc
 * ## $(LABEL): the tag from column 9 to use to label transcript counts, e.g. "Locus", "Name"
 */
process FEATURECOUNTS {

   tag "${sample_id}"

   memory '32GB'

   publishDir( counts_o, 
               mode: 'copy', 
               pattern: "*.tsv" )
   publishDir( counts_o, 
               mode: 'copy', 
               pattern: "*.summary" )

   input:
   tuple val( genome_id ), path( gff ), val( sample_id ), path( bamfile )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   tuple path( "*.summary" ), path( "*.tsv" ), emit: logs

   script:
   """
   featureCounts \
      -s ${params.strand} \
      -d ${params.min_length} -D 3000  \
      -t ${params.ann_type} -g ${params.label} \
      -a ${gff} \
      -R SAM \
      --verbose \
      -o ${sample_id}.featureCounts.tsv ${bamfile}

   samtools sort ${bamfile}.featureCounts.sam -O bam -l 9 -o ${bamfile}.featureCounts.bam
   samtools index ${bamfile}.featureCounts.bam
   rm ${bamfile}.featureCounts.sam
   """
}

/*
 * Count unique UMIs per cell per gene
 */
process UMITOOLS_COUNT {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( bam_idx )

   output:
   tuple val( sample_id ), path( "*.tsv" )

   script:
   """
   umi_tools count \
		--per-gene --per-cell \
		--gene-tag=XT \
		-I ${bamfile} -S ${sample_id}.umitools_count.tsv
   """
}

/*
 * Make log report
 */
process MULTIQC {

   publishDir( multiqc_o, 
               mode: 'copy' )

   input:
   path '*'

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/