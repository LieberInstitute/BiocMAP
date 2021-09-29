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

#  "Major" options (apply to both first and second halves)
sample_opts = c('single', 'paired')
ref_opts = c('hg38', 'hg19')

#  "Minor" options (for first and second halves, respectively)
minor_opts_first = list(
    c('trim_mode'), 
    list(c('skip', 'adaptive', 'force'))
)

minor_opts_second = list(
    c('use_bme', 'with_lambda'),
    list(c('true', 'false'), c('true', 'false'))
)


spec = matrix(
    c("base_dir", "d", 1, "character", "full path to directory to place grid-test files"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

dir.create(opt$base_dir, showWarnings=FALSE)

#  Shuffle minor options randomly
minor_opts_first = list(
    minor_opts_first[[1]],
    lapply(
        minor_opts_first[[2]],
        function(x) sample(x, size=length(x), replace=FALSE)
    )
)

minor_opts_second = list(
    minor_opts_second[[1]],
    lapply(
        minor_opts_second[[2]],
        function(x) sample(x, size=length(x), replace=FALSE)
    )
)


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


index = 1
for (sample_opt in sample_opts) {
    for (ref_opt in ref_opts) {
        for (half in c('first', 'second')) {
            #  Use the appropriate variables for this half
            if (half == 'first') {
                base_command = first_half_base_command
                minor_opts = minor_opts_first
            } else {
                base_command = second_half_base_command
                minor_opts = minor_opts_second
            }
            
            #  Determine command
            log_path = file.path(
                opt$base_dir,
                paste0('grid_test_', half, '_', index, '.log')
            )
        
            command = c(
                paste('#$ -o', log_path),
                paste('#$ -e', log_path),
                base_command,
                paste0('    --sample "', sample_opt, '" \\'),
                paste0('    --reference "', ref_opt, '" \\'),
                paste0('    --annotation "', here('ref'), '" \\'),
                paste0('    -w "', opt$base_dir, '/work_', half, '_', index, '" \\'),
                paste0('    --output "', opt$base_dir, '/out_', half, '_', index, '" \\')
            )
        
            #  Add command options corresponding to the "minor_opts"
            for (i in 1:length(minor_opts[[1]])) {
                #  The possible values this command option takes as an argument,
                #  followed by the particular value to use for this "element" of
                #  the grid test
                opt_values = minor_opts[[2]][[i]]
                this_value = opt_values[(index %% length(opt_values)) + 1]
                
                command_opt = paste0(
                    '    --', minor_opts[[1]][i], ' ', this_value
                )
                
                if (i < length(minor_opts[[1]])) {
                    command_opt = paste0(command_opt, ' \\')
                }
                
                command = c(command, command_opt)
            }
        
            #  Write into shell script and submit as a job
            shell_path = file.path(opt$base_dir, paste0('run_', half, '_', index, '.sh'))
            writeLines(command, con=shell_path)
            #system(paste('qsub', shell_path))
        }
        
        index = index + 1
    }
}
