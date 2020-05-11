#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

--------------------------------------------------------------------
    WGBS pipeline
--------------------------------------------------------------------------

input: SAM files, whose paths are specified in rules.txt
output: bsseq objects, separated by cytosine context

processes:
    1. convert sam to bam
    2. methylation extraction via MethylDackel
    3. HDF5-backed bsseq object creation
*/

def helpMessage() {
    log.info"""
    =========================
        WGBS pipeline
    =========================
    
    Usage:
        nextflow main.nf [options]
    
    Typical use case:
        nextflow main.nf --sample "paired" -profile jhpce
        
    Required flags:
        --sample:      "single" or "paired", depending on your FASTQ reads
        --reference:   "hg38", "hg19", or "mm10". The reference genome to be
                       used for alignment and methylation extraction
    
    Optional flags:
        --input [path]:   the path to the directory containing the rules.txt file
        --use_bme:        include this flag to perform methylation extraction-
                          related processes with Bismark utilities, rather than
                          the default of MethylDackel
    """.stripIndent()
}


// -------------------------------------
//   Define default values for params
// -------------------------------------

params.reference = ""
params.sample = ""
params.input = "${workflow.projectDir}/test"
params.output = "${workflow.projectDir}/out"
params.work = "${workflow.projectDir}/work"
params.use_bme = false

// -------------------------------------
//   Validate Inputs
// -------------------------------------

if (params.sample != "single" && params.sample != "paired") {
    exit 1, "Sample type not provided or invalid. Valid options for --sample are 'single' or 'paired'."
}

if (params.reference != "hg19" && params.reference != "hg38" && params.reference != "mm10") {
    exit 1, "Reference not provided or invalid. Valid options for --reference are 'hg19', 'hg38', or 'mm10'."
}

// ------------------------------------------------------------
//   Construct convenience variables (dependent on reference)
// ------------------------------------------------------------

if (params.reference == "hg38") {
    params.ref_fasta_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz"
} else if (params.reference == "hg19") {
    params.ref_fasta_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
} else {  // mm10
    params.ref_fasta_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz"
}

// ------------------------------------------------------------
//   Utilities for retrieving info from filenames
// ------------------------------------------------------------

def get_prefix(f) {
    //  Remove these regardless of position in the string (note blackListAny is a regular expression)
    blackListAny = ~/_[12]_summary|_[12]_fastqc_data|_success_token|_(trimmed|untrimmed)|_(reverse|forward)_(paired|unpaired)|_R[12]\$(a21|raw|sqm|sqq)|CH[GH]_*O[BT]_|CpG_*O[BT]_|_bedgraph_merged/
    
    //  Remove these if at the end of the file (before the file extension)
    String blackListEnd = "_[12]\\.|_R[12]\\.|_(encode|align)_reads\\.|\\.(c|cfu|txt|gz|sorted|)"

    f.name.toString()
        .replaceAll(blackListAny, "")
        .replaceAll(blackListEnd, ".")
        .tokenize('.')[0]
}

def get_chromosome_name(f) {
    f.name.toString()
        .tokenize('.')[1]
        .replaceAll("chrchr", "chr")
}


// ######################################################
//    Pre-processing steps 
// ######################################################

// Pull the reference fasta for the given reference; build the bisulfite-converted
// genome required by Bismark tools
process PrepareReference {
    storeDir "${workflow.projectDir}/ref/${params.reference}"
    
    input:
        file split_fasta_script from file("${workflow.projectDir}/scripts/split_fasta.sh")
        file encode_ref_script from file("${workflow.projectDir}/scripts/write_configs_encode_ref.R")
    
    output:
        file "$baseName" into BME_genome, MD_genome  // the original reference fasta
        
    shell:
        baseName = file("${params.ref_fasta_link}").getName() - ".gz"
        '''
        #  Pull fasta from GENCODE
        wget !{params.ref_fasta_link}
        gunzip !{baseName}.gz
        
        #  Build the bisulfite genome, needed for bismark (and copy it to publishDir,
        #  circumventing Nextflow's inability to recursively copy)
        !{params.bismark_genome_preparation} --hisat2 --path_to_aligner !{params.hisat2} ./
        mkdir -p !{workflow.projectDir}/ref/!{params.reference}
        cp -R Bisulfite_Genome !{workflow.projectDir}/ref/!{params.reference}/
        '''
}


// ######################################################
//    Begin pipeline
// ######################################################


concordant_sams_out
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .ifEmpty{ error "Concordant sams missing from input to 'FilterSam' process." }
    .set{ concordant_sams_in }
    

//  Extracting alignment information from the AriocP or AriocU logs
process ParseAriocLogs {

    publishDir "${params.output}/Arioc/", mode:'copy'
    
    input:
        file parse_script from file("${workflow.projectDir}/scripts/parse_arioc_logs.R")
        file arioc_logs_in from arioc_logs_out.collect()
        
    output:
        file "alignment_results.rda"
        
    shell:
        '''
        Rscript !{parse_script}
        '''
}
    

//  To be specific, this process takes the sam file of concordant reads from
//  Arioc, filters by mapping quality, and removes duplicate mappings.
process FilterAlignments {

    publishDir "${params.output}/logs/", mode:'copy', pattern:'*.log'
    tag "$prefix"
    
    input:
        set val(prefix), file(sam_file) from concordant_sams_in
        
    output:
        file "${prefix}.cfu.sam" optional true into processed_sams_out
        file "${prefix}.cfu.sorted.bam*" optional true into processed_alignments_out
        file "filter_sam_${prefix}.log"
        
    shell:
        '''
        if [ !{params.use_bme} == "true" ]; then
            #  Quality-filter and deduplicate
            !{params.samtools} view -q 5 -F 0x100 -h !{sam_file} \
                | !{params.samblaster} -r -o !{prefix}.cfu.sam
        else
            #  Quality-filter, deduplicate, compress, and sort
            !{params.samtools} view -q 5 -F 0x100 -h !{sam_file} \
                | !{params.samblaster} -r \
                | !{params.samtools} sort -@ !{task.cpus} -o !{prefix}.cfu.sorted.bam -
                
            #  Index the sorted BAM
            !{params.samtools} index !{prefix}.cfu.sorted.bam
            
        cp .command.log filter_sam_!{prefix}.log
        '''
}

// ############################################################################
//  The following are methylation extraction processes when Bismark Methylation
//  Extractor is selected (via --use_bme)
// ############################################################################

if (params.use_bme) {
    processed_sams_out
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .ifEmpty{ error "Filtered/ deduplicated sams missing from input to 'BME' process." }
        .set{ processed_sams_in }
    
    //  Bismark Methylation Extractor on the quality-filtered, deduplicated sams
    process BME {
    
        publishDir "${params.output}/BME/", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(sam_file) from processed_sams_in
            
        output:
            file "${prefix}/" into BME_outputs
            file "BME_${prefix}.log" into BME_logs_out
            
        shell:
            // BME needs path to samtools if samtools isn't on the PATH
            if (params.samtools != "samtools") {
                flags = "--samtools_path=${params.samtools}"
            } else {
                flags = ""
            }
            
            // paired vs single-end flag
            if (params.sample == "paired") {
                flags += " --paired-end"
            } else {
                flags += " --single-end"
            }
            
            // on multiple cores? BME runs N *additional* threads with the flag
            // "--multicore N", hence the subtraction by 1
            if (task.cpus > 1) {
                flags += " --multicore " + (task.cpus - 1)
            }
            '''
            mkdir !{prefix}
            !{params.BME} !{flags} --gzip -o !{prefix}/ !{sam_file}
            
            cp .command.log BME_!{prefix}.log
            '''
    }
    
    
    //  Aggregate methylation info from BME logs for all samples, into a single .rda file
    process ParseBMELogs {
    
        publishDir "${params.output}/", mode:'copy'
        
        input:
            file BME_logs_in from BME_logs_out.collect()
            file parse_BME_script from file("${workflow.projectDir}/scripts/parse_bme_logs.R")
            
        output:
            file "BME_metrics.rda"
            
        shell:
            '''
            Rscript !{parse_BME_script}
            '''
    }
    
    
    BME_outputs
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .ifEmpty{ error "BME outputs missing from input to 'Bismark2Bedgraph' process." }
        .set{ bedgraph_in }
    
    //  bismark2bedraph on BME outputs. This is done separately because this utility,
    //  BME, and coverage2cytosine all have significantly different hardware resource
    //  requirements. Thus to run the pipeline at maximum speed, while minimizing
    //  wasted resources, the 3 utilities should be separated.
    process Bismark2Bedgraph {
    
        publishDir "${params.output}/Reports/$prefix", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(BME_dir) from bedgraph_in
            
        output:
            file "*_bedgraph_merged*"
            file "${prefix}_bedgraph_merged.gz.bismark.cov.gz" into bedgraph_outputs
            file "bismark2bedgraph_${prefix}.log"
            
        shell:
            '''
            !{params.bismark2bedGraph} -o !{prefix}_bedgraph_merged ./!{BME_dir}/*_!{prefix}.cfu.txt.gz
            
            cp .command.log bismark2bedgraph_!{prefix}.log
            '''
    }
    
    
    bedgraph_outputs
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .ifEmpty{ error "Bedgraphs missing from input to 'Coverage2Cytosine' process." }
        .set{ c2c_in }
    
    process Coverage2Cytosine {
    
        publishDir "${params.output}/Reports/$prefix", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(bedgraph_file) from c2c_in
            
        output:
            file "*.CX_report.txt" into cytosine_reports
            file "C2C_${prefix}.log"
            
        shell:
            '''
            !{params.coverage2cytosine} \
                --split_by_chromosome \
                --CX \
                --genome_folder !{workflow.projectDir}/ref/!{params.reference}/ \
                -o !{prefix} \
                !{bedgraph_file}
                
            cp .command.log C2C_!{prefix}.log
            '''
    }

// ############################################################################
//  The alternative is MethylDackel for methylation extraction
// ############################################################################

} else {
    processed_alignments_out
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .ifEmpty{ error "Processed/ sorted BAMs missing from input to 'MethylationExtraction' process." }
        .set{ processed_alignments_in }
        
    process MethylationExtraction {
        publishDir "${params.output}/MethylDackel/", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(bam_file), file(bam_index) from processed_alignments_in
            file MD_genome
            file meth_count_script from file("${workflow.projectDir}/scripts/meth_count.R")
            
        output:
            file "methyl_extraction_${prefix}.log"
            file "${prefix}*.CX_report.txt" into cytosine_reports
            
        shell:
            '''
            #  Run methylation extraction
            echo "Running 'MethylDackel extract' on the sorted bam..."
            !{params.MethylDackel} extract --cytosine_report !{MD_genome} !{bam_file}
            
            echo "Summary stats for !{prefix}:"
            Rscript !{meth_count_script} -t !{params.data_table_threads}
            
            
            #  Split reports by sequence (pulled from BAM header)
            echo "Splitting cytosine report by sequence..."
            for SN in $(!{params.samtools} view -H !{bam_file} | cut -f 2 | grep "SN:" | cut -d ":" -f 2); do
                awk -v sn=$SN '$1 == sn' !{prefix}.cytosine_report.txt > !{prefix}.$SN.CX_report.txt
            done
            
            echo "Done."
            cp .command.log methyl_extraction_!{prefix}.log
            '''
    }
}

//  Group reports for each chromosome into one channel (each)
cytosine_reports
    .flatten()
    .filter{ get_chromosome_name(it).length() <= 5 } // take only canonical seqs
    .map{ file -> tuple(get_chromosome_name(file), file) }
    .groupTuple()
    .set{ cytosine_reports_by_chr }


process FormBsseqObjects {

    publishDir "${params.output}/BSobjects/logs", mode:'copy', pattern:'*.log'
    tag "$chr"
        
    input:
        set val(chr), file(reports) from cytosine_reports_by_chr
        file bs_creation_script from file("${workflow.projectDir}/scripts/bs_create.R")
        
    output:
        file "assays_${chr}.h5" into assays_out
        file "bs_${chr}_*.rda" into bs_objs_out
        file "create_bs_${chr}.log"
        
    shell:
        '''
        Rscript !{bs_creation_script} -s !{chr} -c !{task.cpus}
        
        cp .command.log create_bs_!{chr}.log
        '''
}


//  Combine Bsseq objects and their HDF5-backed assays into two .rda files
//  (one for CpG context, the other for CpH) and a single .h5 file
process MergeBsseqObjects {

    publishDir "${params.output}/BSobjects/logs", mode:'copy'
    
    input:
        file assays_in from assays_out.collect()
        file bs_objs_in from bs_objs_out.collect()
        file combine_script from file("${workflow.projectDir}/scripts/bs_merge.R")
        
    output:
        file "merge_objects.log"
        
    shell:
        '''
        #  This script actually writes result files directly to publishDir
        Rscript !{combine_script} -d !{params.output}/BSobjects
        
        cp .command.log merge_objects.log
        '''
}