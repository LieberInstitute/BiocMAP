#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

--------------------------------------------------------------------
    WGBS pipeline
--------------------------------------------------------------------------

input: fastq or fastq.gz files, which should be present in ./in/
output: bsseq objects, separated by cytosine context, to be written to ./out/ by default

processes:
    A. Pre-processing on input FASTQ files- this involves decompressing any gzipped
       files, renaming file extensions or paired-end read suffices (if necessary),
       and merging files with the same sample ID
    
    1. FastQC on FASTQ files
    2. Trimming FASTQ files (by default based on adapter content FastQC metric)
    3. FastQC on any trimmed files, to confirm trimming achieved its goal
    4. Alignment with Arioc
    5. Filter to high quality, unique reads, and convert sam to bam
    6. Bismark methylation extraction
    7. bismark2bedgraph on BME outputs
    8. coverage2cytosine on bismark2bedgraph outputs
    9. HDF5-backed bsseq object creation 
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
        --annotation [path]: the path to the directory to store annotation-
                          related files
        --custom_anno [name]: use the FASTA present in the annotation
                          directory, and associate it with a name for future
                          runs
        --input [path]:   the path to the directory containing the 
                          samples.manifest. Defaults to "./test"
        --output [path]:  the directory into which to place pipeline results 
                          and outputs. Defaults to "./out"
        --trim_mode [mode]: determines the conditions under which trimming occurs:
                          "skip": do not perform trimming on samples
                          "adaptive": [default] perform trimming on samples
                              that have failed the FastQC "Adapter content" 
                              metric
                          "force": perform trimming on all samples
        --all_alignments: include this flag to signal Arioc to also write
                          outputs for discondant, rejected, and unmapped reads.
                          Sam files for each outcome are kept as pipeline
                          outputs. By default, only concordant reads are used
                          for later processing (methylation extraction and
                          beyond)
    """.stripIndent()
}


// -------------------------------------
//   Define default values for params
// -------------------------------------

params.all_alignments = false
params.annotation = "${workflow.projectDir}/ref"
params.custom_anno = ""
params.force_trim = false
params.input = "${workflow.projectDir}/test"
params.output = "${workflow.projectDir}/out"
params.reference = ""
params.sample = ""
params.trim_mode = "adaptive"
params.use_bme = false
params.with_lambda = false
params.work = "${workflow.projectDir}/work"

// -------------------------------------
//   Validate Inputs
// -------------------------------------

if (params.sample != "single" && params.sample != "paired") {
    exit 1, "Sample type not provided or invalid. Valid options for --sample are 'single' or 'paired'."
}

if (params.reference != "hg19" && params.reference != "hg38" && params.reference != "mm10") {
    exit 1, "Reference not provided or invalid. Valid options for --reference are 'hg19', 'hg38', or 'mm10'."
}

// Trim mode
if (params.trim_mode != "skip" && params.trim_mode != "adaptive" && params.trim_mode != "force") {
    exit 1, "'--trim_mode' accepts one of three possible arguments: 'skip', 'adaptive', or 'force'."
}

// ------------------------------------------------------------
//   Construct convenience variables (dependent on reference)
// ------------------------------------------------------------

if (params.custom_anno != "") {
    params.anno_version = "custom"
    params.anno_suffix = params.custom_anno + "_custom_build"
} else if (params.reference == "hg38") {
    params.anno_version = params.gencode_version_human
    params.anno_suffix = params.reference + '_gencode_v' + params.gencode_version_human + '_' + params.anno_build
    params.ref_fasta_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.gencode_version_human}/GRCh38.primary_assembly.genome.fa.gz"
} else if (params.reference == "hg19") {
    params.anno_version = params.gencode_version_human
    params.anno_suffix = params.reference + '_gencode_v' + params.gencode_version_human + 'lift37_' + params.anno_build
    params.ref_fasta_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.gencode_version_human}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
} else {  // mm10
    params.anno_version = params.gencode_version_mouse
    params.anno_suffix = params.reference + '_gencode_' + params.gencode_version_mouse + '_' + params.anno_build
    params.ref_fasta_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.gencode_version_mouse}/GRCm38.primary_assembly.genome.fa.gz"
}


// ------------------------------------------------------------
//   Utilities for retrieving info from filenames
// ------------------------------------------------------------

def replace_listed(x, pattern_list, replacement) {
    for (pattern in pattern_list) {
        x = x.replaceAll(pattern, replacement)
    }
    return x
}

def get_prefix(f) {
    //  Remove these regardless of position in the string (note blackListAny is a regular expression)
    blackListAny = [~/_[12]_(summary|fastqc_data)/, ~/_success_token/, ~/_(trimmed|untrimmed)/, ~/_(reverse|forward)/, ~/_(paired|unpaired)/, ~/_R[12]\$(a21|raw|sqm|sqq)/, ~/CH[GH]_*O[BT]_|CpG_*O[BT]_/, ~/|_bedgraph_merged/]
    
    //  Replace these with a dot
    blackListDot = [~/(_[12]|_R[12])\./, ~/_(encode|align)_reads\./, ~/\.(c|cfu)/, ~/\.txt/, ~/\.gz/, ~/\.sorted/, ~/_val/]

    f = replace_listed(f.name.toString(), blackListAny, "")
    f = replace_listed(f, blackListDot, ".")
    
    return f.tokenize('.')[0]
}

def get_chromosome_name(f) {
    f.name.toString()
        .tokenize('.')[1]
        .replaceAll("chrchr", "chr")
}

def get_file_ext(f) {
    if (f.name.toString().tokenize(".")[-1] == "gz") {
        return('.fastq.gz')
    } else {
        return('.fastq')
    }
}

def get_context(f) {
    f.name.toString()
        .tokenize('.')[0][-3..-1]
}


// ------------------------------------------------------------
//  Print all parameters to log
// ------------------------------------------------------------

log.info "=================================="
log.info " WGBS Pipeline"
log.info "=================================="
def summary = [:]
summary['All alignments'] = params.all_alignments
summary['Annotation dir'] = params.annotation
summary['Annotation release'] = params.anno_version
summary['Annotation build'] = params.anno_build
summary['Input dir'] = params.input
summary['Output dir'] = params.output
summary['Reference'] = params.reference
summary['Sample']	= params.sample
summary['Trim mode'] = params.trim_mode
summary['Working dir'] = workflow.workDir
summary['Current user']		= "$USER"
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "==========================================="

// ######################################################
//    Pre-processing steps 
// ######################################################

if (params.custom_anno == "") {
    // Pull the reference fasta for the given reference; subset to "main" sequences
    // if necessary
    process PullReference {
        storeDir "${params.annotation}/${params.anno_suffix}"

        input:
            file split_fasta_script from file("${workflow.projectDir}/scripts/split_fasta.sh")
        output:
            file "$out_fasta" into raw_genome

        shell:
            //  Name of the primary assembly fasta after being downloaded and unzipped
            baseName = file("${params.ref_fasta_link}").getName() - ".gz"

            //  Name the pipeline will use for the primary and main assembly fastas, respectively
            primaryName = "assembly_${params.anno_suffix}.fa".replaceAll("main", "primary")
            mainName = "assembly_${params.anno_suffix}.fa".replaceAll("primary", "main")

            //  Name of fasta to use for this pipeline execution instance
            out_fasta = "assembly_${params.anno_suffix}.fa"
            '''
            mkdir -p !{params.annotation}/!{params.anno_suffix}
            
            #  Pull fasta from GENCODE
            wget !{params.ref_fasta_link}
            gunzip !{baseName}.gz
            mv !{baseName} !{primaryName} # rename for consistency with pipeline naming conventions
            
            #######################################################################
            #  Create the "main" fasta of canonical seqs only
            #######################################################################
            
            if [ !{params.anno_build} == "main" ]; then
                #  Determine how many chromosomes/seqs to keep
                if [ !{params.reference} == "mm10" ]; then
                    num_chrs=22
                else
                    num_chrs=25
                fi
                #  Find the line of the header for the first extra contig (to not
                #  include in the "main" annotation fasta
                first_bad_line=$(grep -n ">" !{primaryName} | cut -d : -f 1 | paste -s | cut -f $(($num_chrs + 1)))
                
                #  Make a new file out of all the lines up and not including that
                sed -n "1,$(($first_bad_line - 1))p;${first_bad_line}q" !{primaryName} > !{mainName}
            fi
            '''
    }
} else {
    Channel.fromPath("${params.annotation}/*.fa")
        .ifEmpty{ error "Cannot find FASTA in annotation directory (and --custom_anno was specified)" }
        .first()  // This proves to nextflow that the channel will always hold one value/file
        .into{ raw_genome }
}

//  Split the fasta into individual files (1 per canonical sequence) and write
//  the configs for encoding the reference with AriocE
process PrepareReference {
    storeDir "${params.annotation}/${params.anno_suffix}"
    
    input:
        file split_fasta_script from file("${workflow.projectDir}/scripts/split_fasta.sh")
        file encode_ref_script from file("${workflow.projectDir}/scripts/write_configs_encode_ref.R")
        file raw_genome
    
    output:
        file "chr_names_${params.anno_suffix}"
        file "encode_ref_gap.cfg" into encode_ref_gap_cfg
        file "encode_ref_nongap.cfg" into encode_ref_nongap_cfg
        file "prepare_ref.log"
        
    shell:
        if (params.custom_anno != "") {
            genome_dirname = params.custom_anno
        } else {
            genome_dirname = params.anno_build
        }
        '''
        #  Make a file containing a list of seqnames
        grep ">" !{raw_genome} | cut -d " " -f 1 | cut -d ">" -f 2 > chr_names_!{params.anno_suffix}
        
        #  Split by sequence, to prepare for encoding reference with Arioc
        bash !{split_fasta_script} !{raw_genome}
        
        #  Write the Arioc configs for encoding the reference
        out_dir=!{params.annotation}/!{params.anno_suffix}
        mkdir -p $out_dir
        Rscript !{encode_ref_script} -r !{params.reference} -d $out_dir
        
        #  Keep a log of what happened so far
        cp .command.log prepare_ref.log
        '''
}


// Arioc requires an encoded reference sequence. This process builds that within the repo,
// if the encoded sequence hasn't been built before.
process EncodeReference {
    storeDir "${workflow.projectDir}/ref/${params.reference}"
    
    input:
        file encode_ref_gap_cfg
        file encode_ref_nongap_cfg
        
    output:
        file ".success" into success_token_ref
        
    shell:
        '''
        !{params.AriocE} !{encode_ref_gap_cfg}
        !{params.AriocE} !{encode_ref_nongap_cfg}
        touch .success
        '''
}


process Merging {

    publishDir "${params.output}/logs", mode:'copy', pattern:'*.log'
    tag "Performing merging if/where necessary"

    input:
        file original_manifest from file("${params.input}/samples.manifest")
        file merge_script from file("${workflow.projectDir}/scripts/preprocess_inputs_first.R")
        
    output:
        file "*.f*q*" into merged_inputs_flat
        file "arioc_samples.manifest" into arioc_manifest
        file "preprocess_inputs_first_half.log" into rules_input

    shell:
        '''
        Rscript !{merge_script}
        
        cp .command.log preprocess_inputs_first_half.log
        '''
}


//  Group both reads together for each sample, if paired-end, and assign each sample a prefix
if (params.sample == "single") {
    merged_inputs_flat
        .flatten()
        .map{file -> tuple(get_prefix(file), file) }
        .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
        .into{ fastqc_untrimmed_inputs; trimming_inputs; lambda_inputs }
} else {
    merged_inputs_flat
        .flatten()
        .map{file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
        .into{ fastqc_untrimmed_inputs; trimming_inputs; lambda_inputs }
}

// ######################################################
//    Begin pipeline
// ######################################################


//  -----------------------------------------
//   Step 1: Run FastQC on each sample
//  -----------------------------------------

process FastQC_Untrimmed {
    tag "$fq_prefix"
    publishDir "${params.output}/FastQC/Untrimmed", mode:'copy'

    input:
        set val(fq_prefix), file(fq_file) from fastqc_untrimmed_inputs 

    output:
        file "*"
        file "*_summary.txt" into fastq_summaries_untrimmed

    shell:
        if (params.sample == "single") {
            copy_command = "cp ${fq_prefix}_fastqc/summary.txt ${fq_prefix}_summary.txt"
            data_command = "cp ${fq_prefix}_fastqc/fastqc_data.txt ${fq_prefix}_fastqc_data.txt"
        } else {
            copy_command = "cp ${fq_prefix}_1_fastqc/summary.txt ${fq_prefix}_1_summary.txt && cp ${fq_prefix}_2_fastqc/summary.txt ${fq_prefix}_2_summary.txt"
            data_command = "cp ${fq_prefix}_1_fastqc/fastqc_data.txt ${fq_prefix}_1_fastqc_data.txt && cp ${fq_prefix}_2_fastqc/fastqc_data.txt ${fq_prefix}_2_fastqc_data.txt"
        }
        '''
        !{params.fastqc} -t !{task.cpus} *.f*q* --extract
        !{copy_command}
        !{data_command}
        '''
}

//  Combine FASTQ files and FastQC result summaries for each sample, to form the input channel for Trimming
if (params.sample == "single") {

    fastq_summaries_untrimmed
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .join(trimming_inputs)
        .ifEmpty{ error "All files (fastQC summaries on untrimmed inputs, and the FASTQs themselves) missing from input to trimming channel." }
        .set{ trimming_inputs }
        
} else { // paired

    fastq_summaries_untrimmed
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .join(trimming_inputs)
        .ifEmpty{ error "All files (fastQC summaries on untrimmed inputs, and the FASTQs themselves) missing from input to trimming channel." }
        .set{ trimming_inputs }
}    


//  -----------------------------------------------------------------------
//   Step 2: Trim FASTQ files if required (or requested via --force_trim)
//  -----------------------------------------------------------------------

process Trimming {

    tag "Prefix: $fq_prefix"
    publishDir "${params.output}/Trimming",mode:'copy'

    input:
        set val(fq_prefix), file(fq_summary), file(fq_file) from trimming_inputs

    output:
        file "${fq_prefix}*_fastqc.{html,zip}" optional true
        file "${fq_prefix}*.f*q*_trimming_report.txt" optional true into trim_reports_out
        file "${fq_prefix}*.fq" into trimming_outputs

    shell:
        file_ext = get_file_ext(fq_file[0])
        trim_args = "--illumina --fastqc --dont_gzip --basename ${fq_prefix}"
        if (params.sample == "paired") {
            trim_args = trim_args + " --paired"
        }
        '''
        #  Determine whether to trim the FASTQ file(s). This is ultimately
        #  controlled by the '--trim_mode' command flag.
        if [ "!{params.trim_mode}" == "force" ]; then
            do_trim=true
        elif [ "!{params.trim_mode}" == "skip" ]; then
            do_trim=false
        elif [ "!{params.sample}" == "single" ]; then
            #  Then '--trim_mode "adaptive"' was selected, and data is single-end
            #  (was fq_summary)
            if [ $(grep "Adapter Content" *summary.txt | cut -f 1)  == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        else
            #  Then '--trim_mode "adaptive"' was selected, and data is paired-end
            result1=$(grep "Adapter Content" !{fq_prefix}_1_summary.txt | cut -c1-4)
            result2=$(grep "Adapter Content" !{fq_prefix}_2_summary.txt | cut -c1-4)
            if [ $result1 == "FAIL" ] || [ $result2 == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        fi
        
        #  Run trimming if required
        if [ "$do_trim" == true ]; then
            trim_galore !{trim_args} *.f*q*
        else
            #  Otherwise rename files (for compatibility downnstream, and to signal to
            #  nextflow to output these files) and decompress as necessary
            if [ ![file_ext] == '.fastq.gz' ]; then
                if [ "!{params.sample}" == "single" ]; then
                    gunzip -c !{fq_prefix}!{file_ext} > !{fq_prefix}_trimmed.fq
                else
                    gunzip -c !{fq_prefix}_1!{file_ext} > !{fq_prefix}_val_1.fq
                    gunzip -c !{fq_prefix}_2!{file_ext} > !{fq_prefix}_val_2.fq
                fi
            else
                if [ "!{params.sample}" == "single" ]; then
                    mv !{fq_prefix}!{file_ext} !{fq_prefix}_untrimmed.fq
                else
                    mv !{fq_prefix}_1!{file_ext} !{fq_prefix}_untrimmed_1.fq
                    mv !{fq_prefix}_2!{file_ext} !{fq_prefix}_untrimmed_2.fq
                fi
            fi
        fi
        
        cp .command.log trimming_!{fq_prefix}.log
        '''
}

//  Pair trimming output FASTQs into a channel, grouped by sample name ("prefix")
if (params.sample == "single") {

    trimming_outputs
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .ifEmpty{ error "Single-end trimming output channel is empty" }
        .into{ ariocE_inputs1; ariocE_inputs2 }
} else {

    trimming_outputs
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .ifEmpty{ error "Paired-end trimming output channel is empty" }
        .into{ ariocE_inputs1; ariocE_inputs2 }
}


//  -------------------------------------------------------------------------
//   Step 4: Alignment with Arioc
//  -------------------------------------------------------------------------

//  This is done separately from the encode/ align steps, as writing the configs
//  uses R (possibly not available on a GPU node, where encoding/alignment is done)
process WriteAriocConfigs {

    publishDir "${params.output}/Arioc/configs/",mode:'copy'
    tag "$fq_prefix"
    
    input:
        set val(fq_prefix), file(fq_file) from ariocE_inputs1
        file encode_reads_script from file("${workflow.projectDir}/scripts/write_config_encode_reads.R")
        file align_reads_script from file("${workflow.projectDir}/scripts/write_config_align.R")
        file arioc_manifest
        
    output:
        file "write_configs_${fq_prefix}.log"
        file "*_encode_reads.cfg" into encode_reads_cfgs
        file "*_align_reads.cfg" into align_reads_cfgs
        
    shell:
        if (params.sample == "paired") {
            exec_name = "AriocP"
        } else {
            exec_name = "AriocU"
        }
        
        encoded_dir = "${workflow.projectDir}/Arioc/temp_encoded_reads"
        
        //  Form strings containing lines to write in the Arioc configs, that
        //  are dependent on parameters from the first-half config
        arioc_opts = '<' + exec_name + ' gpuMask="' + params.gpu_mask + \
                     '" batchSize="' + params.batch_size + \
                     '" verboseMask="0xE0000007">'
        r_opts = "  <R>${workflow.projectDir}/ref/${params.reference}/${params.anno_suffix}</R>"
        nongapped_opts = '  <nongapped seed="' + params.nongapped_seed + \
                         '" ' + params.nongapped_args + '/>'
        gapped_opts = '  <gapped seed="' + params.gapped_seed + '" ' + \
                      params.gapped_args + '/>'
        x_opts = '  <X ' + params.x_args + '/>'
        q_opts = '  <Q filePath="' + encoded_dir + '">'
        '''
        Rscript !{encode_reads_script} \
            -p !{params.sample} \
            -d !{encoded_dir} \
            -x !{fq_prefix}
            
        Rscript !{align_reads_script} \
            -p !{params.sample} \
            -a !{params.all_alignments} \
            -f !{fq_prefix} \
            -o \'!{arioc_opts}\' \
            -g \'!{gapped_opts}\' \
            -n \'!{nongapped_opts}\' \
            -r \'!{r_opts}\' \
            -x \'!{x_opts}\' \
            -q \'!{q_opts}\'
            
        cp .command.log write_configs_!{fq_prefix}.log
        '''
}

//  Combine (possibly trimmed) FASTQ files into a channel with their associated
//  AriocE.cfg file
encode_reads_cfgs
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .join( ariocE_inputs2 )
    .ifEmpty{ error "Input channel for AriocE is empty" }
    .set{ ariocE_merged_inputs }


//  FASTQs must be encoded with AriocE, before the alignment step
process EncodeReads {

    publishDir "${params.output}/Arioc/logs/", mode:'copy', pattern:'*.log'
    tag "$fq_prefix"
    
    input:
        set val(fq_prefix), file(config), file(fq_file) from ariocE_merged_inputs
        
    output:
        file "${fq_prefix}_success_token" into success_tokens_reads
        file "encode_${fq_prefix}.log"
        
    shell:
        '''
        !{params.AriocE} !{fq_prefix}_encode_reads.cfg
        touch !{fq_prefix}_success_token
        
        cp .command.log encode_!{fq_prefix}.log
        '''
}

//  This channel includes encoded reads and their associated Arioc alignment config
success_tokens_reads
    .mix(align_reads_cfgs)
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .groupTuple()
    .ifEmpty{ error "Encoded reads missing from input to 'AlignReads' process." }
    .set{ align_in }

process AlignReads {   
    
    publishDir "${params.output}/Arioc/sams/", mode:'copy', pattern:'*.sam'
    publishDir "${params.output}/Arioc/logs/", mode:'copy', pattern:'*.log'
    tag "$prefix"
    
    input:
        // This indicates the reference exists/ was properly built
        file success_token_ref
        
        set val(prefix), file(cfg_and_token) from align_in
        
    output:
        file "${prefix}.[dru].sam" optional true
        file "${prefix}.c.sam" into concordant_sams_out
        file "${prefix}_alignment.log" into arioc_reports_out
        
    shell:
        if (params.sample == "paired") {
            exec_name = "${params.AriocP}"
        } else {
            exec_name = "${params.AriocU}"
        }
        '''
        #  Run alignment
        !{exec_name} !{prefix}_align_reads.cfg
        cp .command.log !{prefix}_alignment.log
        
        #  Rename sams to be [sampleName].[alignment_type].sam
        for sam in $(ls Arioc*.sam); do
            sam_type=$(echo $sam | cut -d "." -f 2)
            mv $sam !{prefix}.$sam_type.sam
        done
        '''
}

concordant_sams_out
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .ifEmpty{ error "Concordant sams missing from input to 'FilterSam' process." }
    .set{ concordant_sams_in }
    

//  To be specific, this process takes the sam file of concordant reads from
//  Arioc, filters by mapping quality, and removes duplicate mappings.
process FilterAlignments {

    publishDir "${params.output}/FilteredAlignments/logs/", mode:'copy', pattern:'*.log'
    publishDir "${params.output}/FilteredAlignments/sams/", mode:'copy', pattern:'*.sam'
    tag "$prefix"
    
    input:
        set val(prefix), file(sam_file) from concordant_sams_in
        
    output:
        file "${prefix}.cfu.sam"
        file "filter_sam_${prefix}.log"
        
    shell:
        '''
        #  Quality-filter and deduplicate
        !{params.samtools} view -q 5 -F 0x100 -h !{sam_file} \
            | !{params.samblaster} -r -o !{prefix}.cfu.sam
            
        cp .command.log filter_sam_!{prefix}.log
        '''
}


//  Generate a 'rules.txt' file automatically, for use as input to the second
//  module
process MakeRules {
    publishDir "${params.input}", mode:'copy'
    
    //  The input is just used to "connect" this process to the computational
    //  graph, required for nextflow to know when the workflow is complete
    input:
        file rules_input
        
    output:
        file 'rules.txt'
        
    shell:
        txt = "# Automatically generated input to the second module/half\n" + \
              "manifest = ${params.input}/samples.manifest\n" + \
              "sam = ${params.output}/FilteredAlignments/sams/[id].cfu.sam\n" + \
              "arioc_log = ${params.output}/Arioc/logs/[id]_alignment.log\n" + \
              "trim_report = ${params.output}/Trimming/[id].fastq_trimming_report.txt"
        
        '''
        echo -e "!{txt}" > rules.txt
        '''
}
