library('here')
library('getopt')

#  Run several instances of the pipeline, testing combinations of major options
#  (such as "--reference"), to verify the pipeline works under a diverse set of
#  potential use cases.
#
#  To keep the number of instances needed to test small, some input options are
#  assumed not to interact with other options (for example, we assume
#  "--trim_mode" will run the same way with or without 'minor' options like
#  "--with_lambda"). This assumption is not made for 'major' options like
#  "--sample".

#  "Major" options
sample_opts = c('single', 'paired')
ref_opts = c('hg38', 'hg19')

#  "Minor" options
trim_opts = c('skip', 'adaptive', 'force')
bme_opts = c('true', 'false')
lambda_opts = c('true', 'false')

spec = matrix(
    c("base_dir", "d", 1, "character", "full path to directory to place grid-test files"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

dir.create(opt$base_dir, showWarnings=FALSE)


first_half_base_command = c(
    '#$ -cwd',
    '#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G',
    '',
    'module load nextflow',
    'export _JAVA_OPTIONS="-Xms8g -Xmx10g"',
    '',
    paste0('nextflow ', here('first_half.nf'), ' \\'),
    '    -profile first_half_jhpce \\'
)
                       
second_half_base_command = c(
    '#$ -cwd',
    '#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G',
    '',
    'module load nextflow',
    'export _JAVA_OPTIONS="-Xms8g -Xmx10g"',
    '',
    paste0('nextflow ', here('second_half.nf'), ' \\'),
    '    -profile second_half_jhpce \\'
)

###############################################################################

#  Shuffle minor options randomly
trim_opts = sample(trim_opts, size=length(trim_opts), replace=FALSE)
bme_opts = sample(bme_opts, size=length(bme_opts), replace=FALSE)
lambda_opts = sample(lambda_opts, size=length(lambda_opts), replace=FALSE)

index = 1
for (sample_opt in sample_opts) {
    for (ref_opt in ref_opts) {
        #  Determine value of minor options
        trim_opt = trim_opts[(index %% length(trim_opts)) + 1]
        bme_opt = bme_opts[(index %% length(bme_opts)) + 1]
        lambda_opt = lambda_opts[(index %% length(lambda_opts)) + 1]
        
        #  Determine command for first half
        log_path = file.path(
            opt$base_dir,
            paste0('grid_test_first_', index, '.log'
        )
        
        first_half_command = c(
            paste('#$ -o', log_path),
            paste('#$ -e', log_path),
            first_half_base_command,
            paste0('    --trim_mode "', trim_opt, '" \\'),
            paste0('    --sample "', sample_opt, '" \\'),
            paste0('    --reference "', ref_opt, '" \\'),
            paste0('    --annotation "', here('ref'), '" \\'),
            paste0('    -w "', opt$base_dir, '/work_first_', index, '" \\'),
            paste0('    --output "', opt$base_dir, '/out_first_', index, '"')
        )
        
        #  Write into shell script and submit as a job
        shell_path = file.path(opt$base_dir, paste0('run_first_', index, '.sh'))
        writeLines(first_half_command, con=shell_path)
        #system(paste('qsub', shell_path))
        
        #  Determine command for second half
        log_path = file.path(
            opt$base_dir,
            paste0('grid_test_second_', index, '.log'
        )
        
        second_half_command = c(
            paste('#$ -o', log_path),
            paste('#$ -e', log_path),
            second_half_base_command,
            paste0('    --with_lambda=', lambda_opt, ' \\'),
            paste0('    --use_bme=', bme_opt, ' \\'), 
            paste0('    --sample "', sample_opt, '" \\'),
            paste0('    --reference "', ref_opt, '" \\'),
            paste0('    --annotation "', here('ref'), '" \\'),
            paste0('    -w "', opt$base_dir, '/work_second_', index, '" \\'),
            paste0('    --output "', opt$base_dir, '/out_second_', index, '"')
        )
        
        #  Write into shell script and submit as a job
        shell_path = file.path(opt$base_dir, paste0('run_second_', index, '.sh'))
        writeLines(second_half_command, con=shell_path)
        #system(paste('qsub', shell_path))
        
        index = index + 1
    }
}
