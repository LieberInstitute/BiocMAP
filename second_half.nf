#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

--------------------------------------------------------------------
    BiocMAP- Second Half
--------------------------------------------------------------------

input: SAM files, whose paths are specified in rules.txt
output: bsseq objects, separated by cytosine context

processes:
    1. sort and compress SAM alignments into BAM format
    2. methylation extraction via MethylDackel/ Bismark
    3. HDF5-backed bsseq object creation
*/

def helpMessage() {
    log.info"""
    ================================================================================
        BiocMAP- Second Half
    ================================================================================
    
    Usage:
        nextflow second_half.nf [options]
    
    Typical use case:
        nextflow second_half.nf --sample "paired" --reference "hg38" \\
                                -profile jhpce
        
    Required flags:
        --sample:      "single" or "paired", depending on your FASTQ reads
        --reference:   "hg38", "hg19", or "mm10". The reference genome to be
                       used for methylation extraction
    
    Optional flags:
        --annotation [path]: the path to the directory to store annotation-
                          related files
        --custom_anno [name]: use the FASTA present in the annotation
                          directory, and associate it with a name for future
                          runs
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

params.annotation = "${workflow.projectDir}/ref"
params.custom_anno = ""
params.output = "${workflow.projectDir}/out"
params.reference = ""
params.sample = ""
params.use_bme = false
params.with_lambda = false
params.work = "${workflow.projectDir}/work"

if (params.reference == "mm10") {
    params.input = "${workflow.projectDir}/test/mouse/${params.sample}"
} else {
    params.input = "${workflow.projectDir}/test/human/${params.sample}"
}

// Ensure output directory is absolute, without '.' or '..'. See
// https://github.com/LieberInstitute/BiocMAP/issues/31#issuecomment-2912688552
params.output_clean = java.nio.file.Paths.get(params.output).toAbsolutePath().normalize().toString()

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
    blackListAny = [~/_[12]_(summary|fastqc_data)/, ~/_success_token/, ~/_(trimmed|untrimmed)/, ~/_(reverse|forward)/, ~/_(paired|unpaired)/, ~/_R[12]\$(a21|raw|sqm|sqq)/, ~/CH[GH]_*O[BT]_|CpG_*O[BT]_/, ~/_bedgraph_merged/]
    
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

//  Given a "row" of the 'samples.manifest' file as a string, return the sample
//  ids
def get_ids(row) {
    if (params.sample == "single") {
        return(row.tokenize('\t')[2])
    } else {
        return(row.tokenize('\t')[4])
    }
}

//  Given a "row" of the 'samples.manifest' file as a string, return the FASTQ
//  files
def get_fastq_names(row) {
    if (params.sample == "single") {
        return(file(row.tokenize('\t')[0]))
    } else {
        return(tuple(file(row.tokenize('\t')[0]), file(row.tokenize('\t')[2])))
    }
}

//  Given a single line of 'rules.txt', return either a string or a list of
//  strings containing all glob expression(s) (note that instances of '[id]'
//  are still returned, and not "evaluated") 
def get_rules_glob(row) {
    if (row[0] == "#") {
        return ""
    } else {
        return row.replaceAll("\\s", "").tokenize('=')[1]
    }
}

// This gets the SHA commit ID of the repository where BiocMAP is installed.
// This associates the pipeline run with a precise "version" of BiocMAP. Note
// that nextflow provides the "workflow.commitId" variable with this intended
// function- during testing this variable appears to be null.
params.commitId = "git --git-dir=${workflow.projectDir}/.git rev-parse HEAD".execute().text.trim()

def summary_main = [:]
summary_main['BiocMAP version'] = params.commitId
summary_main['Config profile'] = workflow.profile
summary_main['Annotation dir'] = params.annotation
summary_main['Annotation release'] = params.anno_version
summary_main['Annotation build'] = params.anno_build
summary_main['Custom anno label'] = params.custom_anno
summary_main['Input dir'] = params.input
summary_main['Output dir'] = params.output_clean
summary_main['Reference'] = params.reference
summary_main['Sample']	= params.sample
summary_main['Use BME'] = params.use_bme
summary_main['Has lambda spike-ins'] = params.with_lambda
summary_main['Working dir'] = workflow.workDir
summary_main['Current user']		= "$USER"

def summary_args = [:]
summary_args['Using containers'] = params.using_containers
summary_args['Single-end Kallisto args'] = params.kallisto_single_args

//  Write run info to log
log.info "================================================================================"
log.info "    BiocMAP- Second Module"
log.info "================================================================================"
log.info "---- Main options:"
log.info summary_main.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "---- Software arguments:"
log.info summary_args.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "================================================================================"


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
            file "assembly_${params.anno_suffix}.fa" into raw_genome, MD_genome, BME_genome, C2C_genome
            
        shell:
            //  Name of the primary assembly fasta after being downloaded and unzipped
            baseName = file("${params.ref_fasta_link}").getName() - ".gz"
                
            //  Name the pipeline will use for the primary and main assembly fastas, respectively
            primaryName = "assembly_${params.anno_suffix}.fa".replaceAll("main", "primary")
            mainName = "assembly_${params.anno_suffix}.fa".replaceAll("primary", "main")

            '''
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
        .into{ raw_genome; MD_genome; BME_genome; C2C_genome }
}


// Create 'chr_names' file for the current FASTA
process PrepareReference {
    storeDir "${params.annotation}/${params.anno_suffix}"
    
    input:
        file raw_genome
        
    output:
        file "chr_names_${params.anno_suffix}" into chr_names
        
    shell:
        if (params.custom_anno != "") {
            genome_dirname = params.custom_anno
        } else {
            genome_dirname = params.anno_build
        }
        '''
        #  Make a file containing a list of seqnames
        grep ">" !{raw_genome} | cut -d " " -f 1 | cut -d ">" -f 2 > chr_names_!{params.anno_suffix}
        '''
}
        
        
if (params.with_lambda) {
    process PrepareLambda {
        storeDir "${params.annotation}/lambda"
        
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


//  Read rules.txt, producing a list of globs (one glob per row in rules.txt)
rules_file = file("${params.input}/rules.txt")
rules_globs = rules_file
    .readLines()
    .collect{ get_rules_glob(it) }
    .flatten()
    
rules_globs.removeAll('')

//  Read 'samples.manifest' to produce a list of sample IDs
manifest_index = rules_globs
    .collect{ it.tokenize('/')[-1] }
    .indexOf('samples.manifest')
    
manifest_file = file(rules_globs[manifest_index])
ids = manifest_file
    .readLines()
    .collect{ get_ids(it) }

//  Evaluate instances of '[id]' in rules.txt, to produce a list of pure globs
//  matching all required files for the preprocessing step
path_list = []
for (glob in rules_globs) {
    if (glob.contains('[id]')) {
        for (id in ids) {
            path_list.add(glob.replaceAll('\\[id\\]', id))
        }
    } else {
        path_list.add(glob)
    }
}
    
//  Place all files matched by any of the globs into a channel
Channel
    .fromPath(path_list)
    .collect()
    .set{ glob_channel }
    
// Extract FASTQ file paths from the manifest and place in a channel to pass to
// PreprocessInputs
Channel
    .fromPath(manifest_file)
    .splitText()
    .map{ row -> get_fastq_names(row) }
    .flatten()
    .collect()
    .set{ raw_fastqs }

//  Place SAMs and any reports/logs into channels for use in the pipeline
process PreprocessInputs {
    
    publishDir "${params.output_clean}/preprocessing/", mode:'copy', pattern:'preprocess_inputs_second_half.log'
    
    input:
        //  Stage files with unique names, since we are not guaranteed the
        //  basename of each file is unique as it is. See pditommaso's solution
        //  at https://github.com/nextflow-io/nextflow/issues/516
        file "*." from glob_channel
        
        file raw_fastqs
        file rules_file
        file manifest_file
        file preprocess_script from file("${workflow.projectDir}/scripts/preprocess_inputs_second.R")
        
    output:
        file "*.{sam,bam}" into concordant_bams_out
        file "*.bam.bai" optional true into concordant_indices_out
        file "*_{arioc,trim_report,xmc,fastqc}.log" into misc_reports_out
        file "*.f*q*" includeInputs true into fastq_out
        file "preprocess_inputs_second_half.log"
        
    shell:
        '''
        Rscript !{preprocess_script}
        cp .command.log preprocess_inputs_second_half.log
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
        publishDir "${params.output_clean}/lambda", mode:'copy'
        
        input:
            file lambda_indices_in from lambda_indices_out.collect()
            file lambda_genomes_in from lambda_genomes_out.collect()
            set val(prefix), file(fq_file) from fastq_in
            
        output:
            file "${prefix}_lambda_pseudo.log" into lambda_reports_temp
            
        shell:
            '''
            if [ !{params.sample} == 'paired' ]; then
                command_args='!{prefix}_1.f*q* !{prefix}_2.f*q*'
            else
                command_args='!{params.kallisto_single_args} !{prefix}.f*q*'
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
} else {
    lambda_reports_temp = Channel.empty()
}

lambda_reports_temp.ifEmpty('').set{lambda_reports_out}

concordant_bams_out
    .mix(concordant_indices_out)
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .groupTuple()
    .ifEmpty{ error "Concordant sams/bams missing from input to 'MethylationExtraction' or (if applicable) 'BME' processes." }
    .set{ concordant_bams_in }
    

// ############################################################################
//  The following are methylation extraction processes when Bismark Methylation
//  Extractor is selected (via --use_bme)
// ############################################################################

if (params.use_bme) {
    
    //  Bismark Methylation Extractor on the quality-filtered, deduplicated sams
    process BME {
    
        publishDir "${params.output_clean}/BME/", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(sam_file) from concordant_bams_in
            
        output:
            file "${prefix}/" into BME_outputs
            file "methyl_extraction_${prefix}.log" into methyl_reports_out
            
        shell:            
            // paired vs single-end flag
            if (params.sample == "paired") {
                flags = " --paired-end"
            } else {
                flags = " --single-end"
            }
            
            // on multiple cores? BME runs N *additional* threads with the flag
            // "--multicore N", hence the subtraction by 1
            if (task.cpus > 1) {
                flags += " --multicore " + (task.cpus - 1)
            }
            
            is_test = (params.input == "${workflow.projectDir}/test/mouse/${params.sample}") || (params.input == "${workflow.projectDir}/test/human/${params.sample}")
            '''
            if [[ !{is_test} && !{params.use_bme} ]]; then
                #  "Unsort" the sorted test BAM for compatibility with BME
                samtools sort -n -@ !{task.cpus} -o temp.bam !{prefix}.bam
                rm !{prefix}.bam
                mv temp.bam !{prefix}.bam
            fi
            
            mkdir !{prefix}
            bismark_methylation_extractor !{flags} --gzip -o !{prefix}/ $(ls | grep -E "^*.[bs]am$")
            
            cp .command.log methyl_extraction_!{prefix}.log
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
    
        publishDir "${params.output_clean}/Reports/$prefix", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(BME_dir) from bedgraph_in
            
        output:
            file "*_bedgraph_merged*"
            file "${prefix}_bedgraph_merged.gz.bismark.cov.gz" into bedgraph_outputs
            file "bismark2bedgraph_${prefix}.log"
            
        shell:
            '''
            bismark2bedGraph -o !{prefix}_bedgraph_merged ./!{BME_dir}/*.txt.gz
            
            cp .command.log bismark2bedgraph_!{prefix}.log
            '''
    }
    
    
    bedgraph_outputs
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .ifEmpty{ error "Bedgraphs missing from input to 'Coverage2Cytosine' process." }
        .set{ c2c_in }
    
    process Coverage2Cytosine {
    
        publishDir "${params.output_clean}/Reports/$prefix", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(bedgraph_file) from c2c_in
            file C2C_genome
            
        output:
            file "*.CX_report.txt" into cytosine_reports
            file "C2C_${prefix}.log"
            
        shell:
            '''
            coverage2cytosine \
                --split_by_chromosome \
                --CX \
                --genome_folder . \
                -o !{prefix} \
                !{bedgraph_file}
                
            cp .command.log C2C_!{prefix}.log
            '''
    }

// ############################################################################
//  The alternative is MethylDackel for methylation extraction
// ############################################################################

} else {        
    process MethylationExtraction {
        publishDir "${params.output_clean}/Reports/$prefix", mode:'copy'
        tag "$prefix"
        
        input:
            set val(prefix), file(bam_pair) from concordant_bams_in
            file MD_genome
            
        output:
            file "methyl_extraction_${prefix}.log" into methyl_reports_out
            file "${prefix}*.CX_report.txt" into cytosine_reports
            
        shell:
            '''
            #  Run methylation extraction
            echo "Running 'MethylDackel extract' on the sorted bam..."
            MethylDackel extract --cytosine_report -o !{prefix} --CHG --CHH !{MD_genome} *.bam
            
            echo "Summary stats for !{prefix}:"
            c_contexts=("CG" "CHG" "CHH")
            for context in ${c_contexts[@]}; do
                m_perc=$(awk -v ctxt=$context '{if ($6 == ctxt) {M += $4; U += $5}}END{print 100*M/(U+M)}' !{prefix}.cytosine_report.txt)
                echo "C methylated in $context context: ${m_perc}%"
            done
            
            #  Split reports by sequence (pulled from BAM header)
            echo "Splitting cytosine report by sequence..."
            for SN in $(samtools view -H *.bam | cut -f 2 | grep "SN:" | cut -d ":" -f 2); do
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
    .filter{ get_chromosome_name(it).length() <= 5 && !get_chromosome_name(it).equals('lambda') && !get_chromosome_name(it).equals('phiX') } // take only canonical seqs
    .map{ file -> tuple(get_chromosome_name(file), file) }
    .groupTuple()
    .set{ cytosine_reports_by_chr }


//  Take reports and logs from Arioc and optionally Trim Galore, XMC, BME, and
//  lambda pseudoalignment (if applicable); aggregate relevant metrics/
//  QC-related stats into a data frame.
process ParseReports {

    publishDir "${params.output_clean}/metrics", mode:'copy'
    
    input:
        file misc_reports_in from misc_reports_out.collect()
        file methyl_reports_in from methyl_reports_out.collect()
        file lambda_reports_in from lambda_reports_out.collect()
        
        file parse_reports_script from file("${workflow.projectDir}/scripts/parse_reports.R")
        file rules from file("${params.input}/rules.txt")
        
    output:
        file "metrics.rda" into metrics
        file "metrics.log"
        
    shell:
        '''
        Rscript !{parse_reports_script}
        cp .command.log metrics.log
        '''
}


process FormBsseqObjects {

    publishDir "${params.output_clean}/BSobjects/logs", mode:'copy', pattern:'*.log'
    
    // When using docker, publish the bsseq objects to the output folder
    // normally. Otherwise, use a shortcut where we directly write objects to
    // the output folder during creation. The shortcut saves a significant
    // amount of disk space as well as some I/O strain and time
    publishDir "${params.output_clean}/BSobjects/objects/$chr/CpG", mode:'copy', pattern: '*_CpG.{h5,rds}', saveAs: { filename -> filename.replaceAll("_${chr}_CpG", "") }, enabled: params.using_containers
    publishDir "${params.output_clean}/BSobjects/objects/$chr/CpH", mode:'copy', pattern: '*_CpH.{h5,rds}', saveAs: { filename -> filename.replaceAll("_${chr}_CpH", "") }, enabled: params.using_containers
    
    tag "$chr"
        
    input:
        set val(chr), file(reports) from cytosine_reports_by_chr
        file bs_creation_script from file("${workflow.projectDir}/scripts/bs_create.R")
        
    output:
        file "{*.rds,*.h5,${chr}_Cp*.success}" into bsseq_objects_out
        file "create_bs_${chr}.log"
        
    shell:
        '''
        if [[ !{params.using_containers} == "true" ]]; then
            out_dir=$(pwd)
        else
            out_dir=!{params.output_clean}/BSobjects/objects
            mkdir -p ${out_dir}
        fi
        
        Rscript !{bs_creation_script} \
            -s !{chr} \
            -c !{task.cpus} \
            -d ${out_dir}
           
        if [[ !{params.using_containers} == "true" ]]; then
            #  Give unique filenames so that nextflow knows how to manage the
            #  files
            mv !{chr}/CpG/assays.h5 assays_!{chr}_CpG.h5
            mv !{chr}/CpG/se.rds se_!{chr}_CpG.rds
            
            mv !{chr}/CpH/assays.h5 assays_!{chr}_CpH.h5
            mv !{chr}/CpH/se.rds se_!{chr}_CpH.rds
        else
            touch !{chr}_CpG.success
            touch !{chr}_CpH.success
        fi
        
        cp .command.log create_bs_!{chr}.log
        '''
}


//  Group objects (or "success tokens" if applicable) into two elements ("CpG"
//  context and "CpH" context)
bsseq_objects_out
    .flatten()
    .map{ file -> tuple(get_context(file), file) }
    .groupTuple()
    .set{ bsseq_objects_in }

//  Combine Bsseq objects and their HDF5-backed assays into two HDF5-backed
//  summarized experiments (directories)- one for each cytosine context
process MergeBsseqObjects {

    publishDir "${params.output_clean}/BSobjects/logs", pattern: "merge_objects_${context}.log", mode:'copy'
    
    // When using docker, publish the bsseq objects to the output folder
    // normally. Otherwise, use a shortcut where we directly write objects to
    // the output folder during creation. The shortcut saves a significant
    // amount of disk space as well as some I/O strain and time
    publishDir "${params.output_clean}/BSobjects/objects/combined", pattern: '*.{rds,h5}', mode:'move', enabled: params.using_containers
    
    input:
        set val(context), file(bsobj) from bsseq_objects_in
        file chr_names
        file metrics
        file combine_script from file("${workflow.projectDir}/scripts/bs_merge.R")
        
    output:
        file "merge_objects_${context}.log"
        file "*{.rds,.h5}" optional true
        
    shell:
        '''
        if [[ !{params.using_containers} == "true" ]]; then
            #  Organize bsseq objects into directories as they were
            #  originally produced
            for chr in $(cat !{chr_names}); do
                mkdir -p $chr/!{context}
                mv assays_${chr}_!{context}.h5 $chr/!{context}/assays.h5
                mv se_${chr}_!{context}.rds $chr/!{context}/se.rds
            done
            
            in_dir=$(pwd)
            out_dir=$(pwd)
        else
            #  Write files directly to output dir
            in_dir=!{params.output_clean}/BSobjects/objects
            out_dir=!{params.output_clean}/BSobjects/objects/combined
        fi
        
        mkdir -p $out_dir
        
        Rscript !{combine_script} \
          -i $in_dir \
          -o $out_dir \
          -c !{context}

        
        cp .command.log merge_objects_!{context}.log
        '''
}
