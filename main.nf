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
        --output [path]:  the path to the directory to contain pipeline outputs
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
params.use_bme = false
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

def get_context(f) {
    f.name.toString()
        .tokenize('.')[0][-3..-1]
}

//  Write run info to output
log.info "=================================="
log.info " WGBS Pipeline"
log.info "=================================="
def summary = [:]
summary['Sample']	= params.sample
summary['Reference'] = params.reference
summary['Annotation release'] = params.anno_version
summary['Annotation build'] = params.anno_build
summary['Input dir'] = params.input
summary['Output dir'] = params.output
summary['Working dir'] = workflow.workDir
summary['Use BME'] = params.use_bme
summary['Has lambda spike-ins'] = params.with_lambda
summary['Current user']		= "$USER"
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "==========================================="


// ######################################################
//    Pre-processing steps 
// ######################################################

// Pull the reference fasta for the given reference; build the bisulfite-converted
// genome required by Bismark tools
process PrepareReference {
    storeDir "${workflow.projectDir}/ref/${params.reference}/${params.anno_suffix}"
    
    input:
        file split_fasta_script from file("${workflow.projectDir}/scripts/split_fasta.sh")
        file encode_ref_script from file("${workflow.projectDir}/scripts/write_configs_encode_ref.R")
    
    output:
        file "$out_fasta" into BME_genome, MD_genome
        file "chr_names_${params.anno_suffix}" into chr_names
        
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
            
            #  Make a new file out of all the lines up and not including that
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
            !{params.kallisto} index -i lambda_normal.idx lambda.fa
            !{params.kallisto} index -i lambda_bs_artificial.idx lambda_bs_artificial.fa
            
            cp .command.log prepare_lambda.log
            '''
    }
}


//  Place SAMs and any reports/logs into channels for use in the pipeline
process PreprocessInputs {
    
    publishDir "${params.output}/logs/", mode:'copy', pattern:'preprocess_inputs.log'
    
    input:
        file rules from file("${params.input}/rules.txt")
        file preprocess_script from file("${workflow.projectDir}/scripts/preprocess_inputs.R")
        
    output:
        file "*.sam" into concordant_sams_out
        file "*_arioc.log" into arioc_reports_out
        file "*_trim_report_r*.txt" optional true into trim_reports_out
        file "*_xmc.log" optional true into xmc_reports_out
        file "*_bme.log" optional true into bme_reports_out
        file "*.f*q*" optional true into fastq_out
        file "preprocess_inputs.log"
        
    shell:
        '''
        Rscript !{preprocess_script}
        cp .command.log preprocess_inputs.log
        '''
}


// ######################################################
//    Begin pipeline
// ######################################################


if (params.with_lambda) {
    if (params.sample == "single") {
        fastq_out
            .flatten()
            .map{file -> tuple(get_prefix(file), file) }
            .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
            .set{ fastq_in }
            
    } else {
        fastq_out
            .flatten()
            .map{file -> tuple(get_prefix(file), file) }
            .groupTuple()
            .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
            .set{ fastq_in }
    }
    
    process LambdaPseudo {
        
        tag "$prefix"
        publishDir "${params.output}/lambda", mode:'copy'
        
        input:
            file lambda_indices_in from lambda_indices_out.collect()
            file lambda_genomes_in from lambda_genomes_out.collect()
            set val(prefix), file(fq_file) from fastq_in
            
        output:
            file "${prefix}_lambda_pseudo.log" into lambda_reports_out
            
        shell:
            '''
            #  This assumes paired-end samples!! (change later)
            if [ !{params.sample} == 'paired' ]; then
                command_args='!{prefix}_1.f*q* !{prefix}_2.f*q*'
            else
                command_args='--single !{prefix}.f*q*'
            fi
            
            #  Perform pseudoalignment to original and bisulfite-converted genomes
            kallisto quant -t !{task.cpus} -i lambda_normal.idx -o ./orig $command_args
            kallisto quant -t !{task.cpus} -i lambda_bs_artificial.idx -o ./bs $command_args
            
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


concordant_sams_out
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .ifEmpty{ error "Concordant sams missing from input to 'SamToBam' and (if applicable) 'BME' processes." }
    .into{ concordant_sams_in_stb; concordant_sams_in_bme }
    

//  Sort and compress Arioc SAMs
process SamToBam {

    publishDir "${params.output}/logs/", mode:'copy', pattern:'*.log'
    tag "$prefix"
    
    input:
        set val(prefix), file(sam_file) from concordant_sams_in_stb
        
    output:
        file "${prefix}.cfu.sorted.bam*" into processed_alignments_out
        file "sam_to_bam_${prefix}.log"
        
    shell:
        '''
        # Sort and compress
        !{params.samtools} sort -@ !{task.cpus} -o !{prefix}.cfu.sorted.bam !{prefix}.sam
                
        #  Index the sorted BAM
        !{params.samtools} index !{prefix}.cfu.sorted.bam
            
        cp .command.log sam_to_bam_!{prefix}.log
        '''
}

// ############################################################################
//  The following are methylation extraction processes when Bismark Methylation
//  Extractor is selected (via --use_bme)
// ############################################################################

if (params.use_bme) {
    
    //  Bismark Methylation Extractor on the quality-filtered, deduplicated sams
    process BME {
    
        publishDir "${params.output}/BME/", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(sam_file) from concordant_sams_in_bme
            
        output:
            file "${prefix}/" into BME_outputs
            file "BME_${prefix}.log" into BME_reports
            
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
        publishDir "${params.output}/Reports/$prefix", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(bam_pair) from processed_alignments_in
            file MD_genome
            
        output:
            file "methyl_extraction_${prefix}.log"
            file "${prefix}*.CX_report.txt" into cytosine_reports
            
        shell:
            '''
            #  Run methylation extraction
            echo "Running 'MethylDackel extract' on the sorted bam..."
            !{params.MethylDackel} extract --cytosine_report -o !{prefix} --CHG --CHH !{MD_genome} *.bam
            
            echo "Summary stats for !{prefix}:"
            c_contexts=("CG" "CHG" "CHH")
            for context in ${c_contexts[@]}; do
                m_perc=$(awk -v ctxt=$context '{if ($6 == ctxt) {U += $4; M += $5}}END{print 100*M/(U+M)}' !{prefix}.cytosine_report.txt)
                echo "C methylated in $context context: ${m_perc}%"
            done
            
            #  Split reports by sequence (pulled from BAM header)
            echo "Splitting cytosine report by sequence..."
            for SN in $(!{params.samtools} view -H *.bam | cut -f 2 | grep "SN:" | cut -d ":" -f 2); do
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


//  Take reports and logs from Arioc and optionally Trim Galore, XMC, BME, and
//  lambda pseudoalignment (if applicable); aggregate relevant metrics/
//  QC-related stats into a data frame.
process ParseReports {

    publishDir "${params.output}/metrics", mode:'copy'
    
    input:
        file trim_reports_in from trim_reports_out.collect()
        file xmc_reports_in from xmc_reports_out.collect()
        file bme_reports_in from bme_reports_out.collect()
        file arioc_reports_in from arioc_reports_out.collect()
        file lambda_reports_in from lambda_reports_out.collect()
        
        file parse_reports_script from file("${workflow.projectDir}/scripts/parse_reports.R")
        file rules from file("${params.input}/rules.txt")
        
    output:
        file "metrics.rda"
        file "metrics.log"
        
    shell:
        '''
        Rscript !{parse_reports_script}
        cp .command.log metrics.log
        '''
}


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

//  Combine Bsseq objects and their HDF5-backed assays into two HDF5-backed
//  summarized experiments (directories)- one for each cytosine context
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
