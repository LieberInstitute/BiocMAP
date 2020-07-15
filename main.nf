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
        --use_bme:        include this flag to perform methylation extraction-
                          related processes with Bismark utilities, rather than
                          the default of MethylDackel
        --with_lambda:    include this flag if all samples have spike-ins with
                          the lambda bacteriophage genome. Pseudoalignment will
                          then be performed to estimate bisulfite conversion
                          efficiency
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
params.force_trim = false
params.all_alignments = false
params.use_bme = false
params.trim_mode = "adaptive"
params.with_lambda = false

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

if (params.reference == "hg38") {
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

params.lambda_link = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/840/245/GCA_000840245.1_ViralProj14204/GCA_000840245.1_ViralProj14204_genomic.fna.gz"


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
    blackListDot = [~/(_[12]|_R[12])\./, ~/_(encode|align)_reads\./, ~/\.(c|cfu)/, ~/\.txt/, ~/\.gz/, ~/\.sorted/]

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


// ######################################################
//    Pre-processing steps 
// ######################################################

// Pull the reference fasta for the given reference; build the bisulfite-converted
// genome required by Bismark tools; split the fasta into individual files (1 per
// canonical sequence) and write the configs for encoding the reference with AriocE
process PrepareReference {
    storeDir "${workflow.projectDir}/ref/${params.reference}/${params.anno_suffix}"
    
    input:
        file split_fasta_script from file("${workflow.projectDir}/scripts/split_fasta.sh")
        file encode_ref_script from file("${workflow.projectDir}/scripts/write_configs_encode_ref.R")
    
    output:
        file "$out_fasta" into BME_genome, MD_genome
        file "chr_names_${params.anno_suffix}" into chr_names
        file "encode_ref_gap.cfg" into encode_ref_gap_cfg
        file "encode_ref_nongap.cfg" into encode_ref_nongap_cfg
        
    shell:
        //  Name of the primary assembly fasta after being downloaded and unzipped
        baseName = file("${params.ref_fasta_link}").getName() - ".gz"
            
        //  Name the pipeline will use for the primary and main assembly fastas, respectively
        primaryName = "assembly_${params.anno_suffix}.fa".replaceAll("main", "primary")
        mainName = "assembly_${params.anno_suffix}.fa".replaceAll("primary", "main")
            
        //  Name of fasta to use for this pipeline execution instance
        out_fasta = "assembly_${params.anno_suffix}.fa"
        
        '''
        mkdir -p !{workflow.projectDir}/ref/!{params.reference}/!{params.anno_suffix}
        
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
            
            #  Make a new fasta out of all the lines up to and not including that
            sed -n "1,$(($first_bad_line - 1))p;${first_bad_line}q" !{primaryName} > !{mainName}
            
            #  Make a file containing a list of seqnames
            grep ">" !{mainName} | cut -d " " -f 1 | cut -d ">" -f 2 > chr_names_!{params.anno_suffix}
            
            #  Create a link in an isolated directory for compatibility with bismark genome preparation
            mkdir main
            cd main
            ln -s  ../!{mainName} !{mainName}
            cd ..
        else
            #  Create a link in an isolated directory for compatibility with bismark genome preparation
            mkdir primary
            cd primary
            ln -s ../!{primaryName} !{primaryName}
            cd ..
            
            #  Make a file containing a list of seqnames
            grep ">" !{primaryName} | cut -d " " -f 1 | cut -d ">" -f 2 > chr_names_!{params.anno_suffix}
        fi
                    
        #  Build the bisulfite genome, needed for bismark (and copy it to publishDir,
        #  circumventing Nextflow's inability to recursively copy)
        !{params.bismark_genome_preparation} --hisat2 --path_to_aligner !{params.hisat2} ./!{params.anno_build}
        cp -R !{params.anno_build}/Bisulfite_Genome !{workflow.projectDir}/ref/!{params.reference}/!{params.anno_suffix}/
        
        #  Split by sequence, to prepare for encoding reference with Arioc
        bash !{split_fasta_script} !{baseName}
        
        #  Write the Arioc configs for encoding the reference
        out_dir=!{workflow.projectDir}/ref/!{params.reference}/!{params.anno_suffix}/encoded_ref
        mkdir -p $out_dir
        Rscript !{encode_ref_script} -r !{params.reference} -d $out_dir
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


if (params.with_lambda) {
    process PrepareLambda {
        storeDir "${workflow.projectDir}/ref/lambda"
        
        input:
            file convert_script from file("${workflow.projectDir}/scripts/bisulfite_convert.R")
            
        output:
            file "lambda*.idx" into lambda_indices_out
            file "lambda*.fa" into lambda_genomes_out
            file "prepare_lambda.log"
            
        shell:
            baseName = file("${params.lambda_link}").getName() - ".gz"
            '''
            #  Pull, unzip, and simplify sequence name to 'lambda'
            wget !{params.lambda_link}
            gunzip -c !{baseName}.gz | sed '1s/>.*Es/>lambda Es/' > lambda.fa
            
            #  Artificially create a "perfectly bisulfite-converted" fasta
            Rscript !{convert_script}
            
            #  Create a kallisto index for each version of the genome
            kallisto index -i lambda_normal.idx lambda.fa
            kallisto index -i lambda_bs_artificial.idx lambda_bs_artificial.fa
            
            cp .command.log prepare_lambda.log
            '''
    }
}


process Merging {

    publishDir "${params.output}", mode:'copy', pattern:'*.log'
    tag "Performing merging if/where necessary"

    input:
        file original_manifest from file("${params.input}/samples.manifest")
        file merge_script from file("${workflow.projectDir}/scripts/preprocess_inputs.R")
        
    output:
        file "*.f*q*" into merged_inputs_flat
        file "arioc_samples.manifest" into arioc_manifest
        file "preprocess_inputs.log"

    shell:
        '''
        Rscript !{merge_script}
        
        cp .command.log preprocess_inputs.log
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

if (params.with_lambda) {
    process LambdaPseudo {
        
        tag "$prefix"
        publishDir "${params.output}/lambda", mode:'copy'
        
        input:
            file lambda_indices_in from lambda_indices_out.collect()
            file lambda_genomes_in from lambda_genomes_out.collect()
            set val(prefix), file(fq_file) from lambda_inputs
            
        output:
            file "${prefix}_lambda_pseudo.log" into lambda_reports_out
            
        shell:
            '''
            #  This assumes paired-end samples!! (change later)
            fq1=!{prefix}_1.f*q*
            fq2=!{prefix}_2.f*q*
            
            #  Perform pseudoalignment to original and bisulfite-converted genomes
            kallisto quant -i lambda_normal.idx -o ./orig $fq1 $fq2
            kallisto quant -i lambda_bs_artificial.idx -o ./bs $fq1 $fq2
            
            #  Get the counts for number of successfully "aligned" reads
            orig_count=$(sed -n '2p' orig/abundance.tsv | cut -f 4)
            bs_count=$(sed -n '2p' bs/abundance.tsv | cut -f 4)

            #  Write the estimated conversion efficiency to the log
            conv_eff=$(Rscript -e "100 * $bs_count/($orig_count + $bs_count)" | cut -d ']' -f 2)
            echo "Conversion efficiency:${conv_eff}%."
            cp .command.log !{prefix}_lambda_pseudo.log
            '''
    }
}


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
        file "${fq_prefix}*.f*q*_trimming_report.txt" optional true
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
        '''
        refDir=!{workflow.projectDir}/Arioc/encoded_ref/!{params.reference}
        mkdir -p $refDir/temp_encoded_reads
        Rscript !{encode_reads_script} \
            -p !{params.sample} \
            -d !{workflow.projectDir}/Arioc/temp_encoded_reads \
            -x !{fq_prefix}
        Rscript !{align_reads_script} \
            -p !{params.sample} \
            -d !{workflow.projectDir} \
            -r !{params.reference} \
            -b !{params.AriocBatchSize} \
            -a !{params.all_alignments} \
            -x !{fq_prefix}
            
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
        file "${prefix}_alignment.log" into arioc_logs_out
        
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
        fi
            
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
            !{params.MethylDackel} extract --cytosine_report --CHG --CHH !{MD_genome} !{bam_file}
            
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
        file "${chr}_Cp*.success" into bs_tokens_out
        file "create_bs_${chr}.log"
        
    shell:
        '''
        mkdir -p !{params.output}/BSobjects/objects
        Rscript !{bs_creation_script} \
            -s !{chr} \
            -c !{task.cpus} \
            -d !{params.output}/BSobjects/objects
        if [ "$?" == "0" ]; then
            touch !{chr}_CpG.success
            touch !{chr}_CpH.success
        fi
        cp .command.log create_bs_!{chr}.log
        '''
}


//  Group tokens into two elements ("CpG" context and "CpH" context)
bs_tokens_out
    .flatten()
    .map{ file -> tuple(get_context(file), file) }
    .groupTuple()
    .set{ bs_tokens_in }

//  Combine Bsseq objects and their HDF5-backed assays into two .rda files
//  (one for CpG context, the other for CpH) and a single .h5 file
process MergeBsseqObjects {

    publishDir "${params.output}/BSobjects/logs", mode:'copy'
    
    input:
        set val(context), file(token) from bs_tokens_in
        file chr_names
        file combine_script from file("${workflow.projectDir}/scripts/bs_merge.R")
        
    output:
        file "merge_objects_${context}.log"
        
    shell:
        '''
        #  This script actually writes result files directly to publishDir
        Rscript !{combine_script} -d !{params.output}/BSobjects/objects -c !{context}
        
        cp .command.log merge_objects_!{context}.log
        '''
}
