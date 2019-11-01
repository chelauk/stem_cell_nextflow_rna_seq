#!/usr/bin/env nextflow
/*
========================================================================================
                        chela's adaptation of nf-core/rnaseq
========================================================================================
 chela_rnaseq/rnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/chelauk
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    in the working directory set up directory structure as follows
    sample_1/
├── replicate_1
│   ├── sample_1-ID-1_qorts
│   ├── dupradar
│   ├── fastqc
│   ├── fastq_raw
│   │   ├── sample_1-ID-1_L001_R1_001.fastq.gz
│   │   ├── sample_1-ID-1_L001_R2_001.fastq.gz
│   │   ├── sample_1-ID-1_L002_R1_001.fastq.gz
│   │   ├── sample_1-ID-1_L002_R2_001.fastq.gz
│   │   ├── sample_1-ID-1_L003_R1_001.fastq.gz
│   │   ├── sample_1-ID-1_L003_R2_001.fastq.gz
│   │   ├── sample_1-ID-1_L004_R1_001.fastq.gz
│   │   └── sample_1-ID-1_L004_R2_001.fastq.gz
│   ├── mark_duplicates
│   ├── STAR
│   └── trim_galore
├── replicate_2
│   ├── sample_1-ID-2_qorts
│   ├── dupradar
│   ├── fastqc
│   ├── fastq_raw
│   │   ├── sample_1-ID-2_L001_R1_001.fastq.gz
│   │   ├── sample_1-ID-2_L001_R2_001.fastq.gz
│   │   ├── sample_1-ID-2_L002_R1_001.fastq.gz
│   │   ├── sample_1-ID-2_L002_R2_001.fastq.gz
│   │   ├── sample_1-ID-2_L003_R1_001.fastq.gz
│   │   ├── sample_1-ID-2_L003_R2_001.fastq.gz
│   │   ├── sample_1-ID-2_L004_R1_001.fastq.gz
│   │   └── sample_1-ID-2_L004_R2_001.fastq.gz
│   ├── mark_duplicates
│   ├── STAR
│   └── trim_galore
└── replicate_3
    ├── sample_1-ID-3_qorts
    ├── dupradar
    ├── fastqc
    ├── fastq_raw
    │   ├── sample_1-ID-3_L001_R1_001.fastq.gz
    │   ├── sample_1-ID-3_L001_R2_001.fastq.gz
    │   ├── sample_1-ID-3_L002_R1_001.fastq.gz
    │   ├── sample_1-ID-3_L002_R2_001.fastq.gz
    │   ├── sample_1-ID-3_L003_R1_001.fastq.gz
    │   ├── sample_1-ID-3_L003_R2_001.fastq.gz
    │   ├── sample_1-ID-3_L004_R1_001.fastq.gz
    │   └── sample_1-ID-3_L004_R2_001.fastq.gz
    ├── mark_duplicates
    ├── STAR
    └── trim_galore
sample_2/etc

* note you only need to create and populate the sample/replicate/fastq_raw directories
the pipeline will create the rest (sample_ID_qorts/dupradar/fastqc/etc)
    nextflow run chela_rnaseq samples --genome GRCh38 -c sge.config
    Mandatory arguments:
      --genome                      Path to input data (must be surrounded with quotes)
    Optional arguments:
      --fusion                      Check RNA-seq for fusion
      --kallisto                    Run Kallisto
      --Hisat                       Run Hisat
      --fromBam                     runs from bam onward

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index                  Path to STAR index
      --hisat2_index                Path to HiSAT2 index
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF file
      --gff                         Path to GFF3 file
      --bed12                       Path to bed12 file
    
    Other options:
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --seqCenter                   Add sequencing center in @RG line of output BAM header
    """.stripIndent()
}

params.help = false
params.genome = false
params.skipFusion = false
params.skipFastqc = false
params.skipTrimFastq = false
params.skipMakeBam = false
params.star_index = params.genome ? params.genomes[ params.genome  ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome  ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome  ].gtf ?: false : false
params.ref_dir = params.genome ? params.genomes[ params.genome  ].ref_dir ?: false : false

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// def emptyMap = [:] the colon inside the square brackets denotes that the
// object is a map and not a list
def summary = [:]
 if(params.aligner == 'star'){
     summary['Aligner'] = "STAR"  // maps "Aligner" to "STAR"
     if(params.star_index) summary['STAR Index'] = params.star_index
     if(params.fasta) summary['FASTA'] = params.fasta
     if(params.gtf) summary['GTF'] = params.gtf
     if(params.ref_dir) summary['resource directory'] = params.ref_dir
     if(params.skipFusion) summary['skipFusion'] = params.skipFusion
     if(params.skipFastqc) summary['skipFastqc'] = params.skipFastqc
     if(params.skipTrimFastq) summary['skipTrimFastq'] = params.skipTrimFastq
     if(params.skipMakeBam) summary['skipMakeBam'] = params.skipMakeBam
 }

 // Validate inputs

 if( params.star_index && params.aligner == 'star' ){
     star_index = file(params.star_index)
     gtf = file(params.gtf)
     fasta = file(params.fasta)
     ref_dir = file(params.ref_dir)
 }

 log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
 log.info "========================================="

if (params.fastqs == 4){
    Channel
      .fromFilePairs("*/replicate_*/fastq_raw/*L00[1-4]_R{1,2}*fastq.gz",size:8 )
      .map { prefix, file -> tuple(prefix, getSampleID(file[0]), getReplicateID(file[0]),file) }
      .into {fastqc_ch;raw_reads_trimgalore;file_locations }
} else {
    Channel
      .fromFilePairs("*/replicate_*/fastq_raw/*L00[1-4]_R{1,2}*fastq.gz",size:2 )
      .map { prefix, file -> tuple(prefix, getSampleID(file[0]), getReplicateID(file[0]),file) }
      .into {fastqc_ch;raw_reads_trimgalore;file_locations }

}
 // Added channel to check if bam is already created
 
 Channel
      .fromPath("*/replicate_*/bam/*bam")
      .map { file -> tuple(getPrefix_bam(file), getSampleID_bam(file), getReplicateID_bam(file),file) }
      .into { ready_fusion; ready_markDuplicates }
 
def getPrefix_bam( file ){
    //using RegEx to extract prefix
     regexpPE = /([\w_\-]+)\/(\w+_[1-6])\/bam\/([\w\-]+S\d+)/
     (file =~ regexpPE)[0][3]
    }
 
 def getSampleID_bam( file ){
     // using RegEx to extract the SampleID
     regexpPE = /([\w_\-]+)\/(\w+_[1-6])\/\w+\/[\w\-]+/
     (file =~ regexpPE)[0][1]
 }
 
 def getReplicateID_bam( file ){
     // using RegEx to extract the SampleID
     regexpPE = /([\w_\-]+)\/(\w+_[1-6])\/\w+\/[\w\-]+/
     (file =~ regexpPE)[0][2]
 }
 
 def getSampleID( file ){
     // using RegEx to extract the SampleID
     regexpPE = /([\w_\-]+)\/(\w+_[1-6])\/\w+\/[\w\-]+(L[0-4]{3})/
     (file =~ regexpPE)[0][1]
 }
 
 def getReplicateID( file ){
     regexpPE = /([\w_\-]+)\/(\w+_[1-6])\/\w+\/[\w\-]+(L[0-4]{3})/
     (file =~ regexpPE)[0][2]
 }

process fastqc {
  tag "${sample_id},${replicate}"
  echo true
  label 'small_mem'

  publishDir "${sample_id}/${replicate}/fastqc", mode: 'copy',
     saveAs: {filename -> filename.indexOf(".zip") > 0 ? "$filename" : "$filename"}
  
  when:
  !params.skipFastqc

  input:
  set sample_prefix, sample_id, replicate, file(reads) from fastqc_ch

  output:
  file "${sample_prefix}*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc -q $reads
  """
  }

process trim_fastq {
  tag "${sample_id},${replicate}"
  echo true
  label 'small_mem'

  publishDir "${sample_id}/${replicate}/trim_galore", mode: 'copy',
    saveAs: {filename ->
       if (filename.indexOf("trimming_report.txt") > 0) "$filename"
       else null
  }
  
  when:
  !params.skipTrimFastq

  input:
  set val(sample_prefix), val(sample_id), val(replicate), file(reads) from raw_reads_trimgalore

  output:
  set val(sample_prefix), val(sample_id), val(replicate), file("${sample_prefix}*fq.gz") into trimmed_ch
  file("*txt") into trimmed_report

  script:
  """
  trim_galore --gzip --paired $reads
  """
}

process make_bam {
  tag "${sample_prefix}, ${replicate}"
  label 'big_mem'

  publishDir "${sample_id}/${replicate}/STAR", mode: 'copy',
    saveAs: { filename ->
    if (filename.indexOf("Log.out") > 0 ) "$filename"
    else if (filename.indexOf("Log.final.out") > 0 ) "$filename"
    else if (filename.indexOf("bam") > 0 ) "$filename"
    else if (filename.indexOf("bai") > 0 ) "$filename"
    else if (filename.indexOf("idxstats") > 0 ) "$filename"
    else null
    }
  
  when:
  !params.skipMakeBam

  input:
  set val(sample_prefix), val(sample_id), val(replicate), file(reads) from trimmed_ch
  file index from star_index

  output:
  set val(sample_prefix), val(sample_id), val(replicate), file("${sample_prefix}*bam") into bam_res, mark_duplicates
  file("${sample_prefix}*Log*") into bam_log
  file("${sample_prefix}*bai") into bam_idx
  file("${sample_prefix}*idxstats") into samtools_stats
  set val(sample_prefix), val(sample_id), val(replicate), file("${sample_prefix}*bam") into junction_ch

  script:
  R1 = reads[0,2,4,6]
  R2 = reads[1,3,5,7]
  """
  STAR --genomeDir ${index} \\
       --readFilesIn ${R1.join(',')} ${R2.join(',')} \\
       --readFilesCommand zcat \\
       --twopassMode Basic \\
       --outReadsUnmapped None \\
       --chimSegmentMin 12 \\
       --chimJunctionOverhangMin 12 \\
       --alignSJDBoverhangMin 10 \\
       --alignMatesGapMax 100000 \\
       --alignIntronMax 100000 \\
       --chimSegmentReadGapMax 3 \\
       --alignSJstitchMismatchNmax 5 -1 5 5 \\
       --runThreadN 2 \\
       --outSAMstrandField intronMotif \\
       --outSAMunmapped Within \\
       --outSAMtype BAM SortedByCoordinate \\
       --chimMultimapScoreRange 10 \\
       --chimNonchimScoreDropMin 10 \\
       --outFileNamePrefix ${sample_prefix}_ \\
       --outSAMattrRGline ID:${sample_prefix} CN:NTRGL DS:RNA-seq SM:${sample_id} PL:illumina \\
       --chimOutType WithinBAM \\
       --chimOutJunctionFormat 1

  samtools index ${sample_prefix}_Aligned.sortedByCoord.out.bam
  samtools idxstats ${sample_prefix}_Aligned.sortedByCoord.out.bam > ${sample_prefix}_idxstats
  """
 }

process fusion {

    tag "${sample_id}.fusion"
    errorStrategy 'ignore'
    label 'medium_mem'
 
    publishDir "${sample_id}/${replicate}/fusion", pattern: '*fusion*', mode: 'copy'
    
    when:
    !params.skipFusion
    
    input:
    set val(sample_prefix), val(sample_id), val(replicate), file(bam) from junction_ch.mix(ready_fusion)
    file gtf from gtf
    file fasta from fasta

    output:
    file("*fusion*") into fus_out
 
    script:
    """
    arriba \\
    -x ${bam} \\
    -g ${gtf} \\
    -o ${sample_prefix}.fusions.tsv \\
    -a ${fasta} \\
    -b ~/applications/arriba_v1.1.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv \\
    -T \\
    -P
    """
}
 
process markDuplicates {

    label 'small_mem'

    tag "${sample_prefix}.md"

    publishDir "${sample_id}/${replicate}/mark_duplicates", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

    input:
    set val(sample_prefix), val(sample_id), val(replicate), file(bam) from mark_duplicates.mix(ready_markDuplicates)

    output:
    set val(sample_prefix), val(sample_id), val(replicate), file("${sample_prefix}*markDups.bam") into bam_quali,bam_qc,bam_md,bam_stringtie
    file "${sample_prefix}*markDups_metrics.txt" into picard_results
    file "${sample_prefix}*markDups.bam.bai"

    """
    java -Xmx16g -jar /home/sejjctj/applications/picard/build/libs/picard.jar MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${bam.baseName}.markDups.bam \\
        METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${bam.baseName}.markDups.bam
    """
}


process stringtie {
    errorStrategy 'ignore'
    tag "${sample_prefix}, ${replicate}"
    publishDir "${sample_id}/${replicate}", mode: 'copy'
 
    input:
    set val(sample_prefix), val(sample_id), val(replicate), file(bam) from bam_stringtie

    output:
    set val(sample_prefix), val(sample_id), val(replicate), file("stringtie") into stringtie_out

    """
    mkdir stringtie
    ~/applications/stringtie-2.0.4.Linux_x86_64/stringtie -p 4 \\
    -G  ${gtf}  \\
    -A  ./stringtie/${sample_id}.gene_abund.tab \\
    -o  ./stringtie/${sample_id}.gtf \\
    -B -e \\
    -v ${bam}
    """
}

process qorts {
    tag "${sample_prefix}, ${replicate}"
    errorStrategy 'ignore'
    publishDir "${sample_id}/${replicate}/", mode: 'copy'

    input:
    set val(sample_prefix), val(sample_id), val(replicate), file(bam) from bam_qc

    output:
    file("${sample_prefix}_qorts") into qorts_res

    script:
    """
    java -Xmx32G -jar ~/applications/qorts/QoRTs-STABLE.jar QC \\
    --minMAPQ 60 \\
    --maxReadLength 150 \\
    --prefilterImproperPairs \\
    ${bam} \\
    ${gtf} \\
    ${sample_prefix}_qorts
    """
}

 process qualimap {
     tag "${sample_prefix}, ${replicate}"
     errorStrategy 'ignore'
     publishDir "${sample_id}/${replicate}/", mode: 'copy'
 
     input:
     set val(sample_prefix), val(sample_id), val(replicate), file(bam) from bam_quali
 
     output:
     file("qualimap") into qualimap_res
 
     script:
     """
     mkdir java_temp
     export JAVA_OPTS="-Djava.io.tmpdir=./java_temp"
     mkdir qualimap
     ~/applications/qualimap_v2.2.1/qualimap rnaseq  \\
     --java-mem-size=16G \\
     -bam ${bam} \\
     -gtf ${gtf} \\
     -pe \\
     -outdir ./qualimap
     """
 }

process QoRTsR {
    scratch 'false'
    module 'r/new'
    echo true
    errorStrategy 'ignore'
    
    publishDir ".", mode: 'copy'

    input:
    file(w) from qorts_res.collect()
    file(x) from Channel.fromPath('samples')
    file(y) from Channel.fromPath('/home/sejjctj/useful/get_qorts_dataver2')

    output:
    file('*pdf')

    script:
    """
    bash ${y} ${x} 1 star > decoder
    Rscript /home/sejjctj/pipelines/rna_seq/qorts.R decoder
    """
}

process dupradar {
    
    errorStrategy 'retry'
    maxErrors 5
    module 'r/new'
    tag "${sample_prefix}.dupRadar"
    publishDir "${sample_id}/${replicate}/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    input:
    set val(sample_prefix), val(sample_id), val(replicate), file(bam) from bam_md
    file gtf

    output:
    file "*.{pdf,txt}" into dupradar_results

    """
    Rscript  /home/sejjctj/pipelines/rna_seq/dupRadar.R ${bam} ${gtf} 0 TRUE 4
    """
}
