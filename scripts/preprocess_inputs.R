## Required libraries
library('devtools')
library('jaffelab')

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

get_value = function(rules, key) {
    this_line = rules[grep(key, ss(rules, '='))]
    if (length(this_line) != 1) {
        stop(paste0("Key '", key, "' had ", length(this_line), " values in rules.txt."))
    }
    
    return(gsub(' ', '', ss(this_line, '=', 2)))
}

rules = readLines('rules.txt')

man_path = file.path(get_value(rules, 'base_dir'), get_value(rules, 'manifest'))
manifest = read.table(man_path, header = FALSE, stringsAsFactors = FALSE)
ids = manifest[ncol(manifest)]

#  Get the paths to the Arioc SAMs, substituting in actual ids for each '[id]'
pieces = ss(get_value(rules, 'sam'), '[id]', fixed=TRUE)
arioc_sams = ''
for (i in 1:(length(pieces)-1)) {
    arioc_sams = paste0(pieces[i], arioc_sams)
}
arioc_sams = paste0(arioc_sams, pieces[length(pieces)])
stopifnot(length(arioc_sams) == length(ids))

#  Symbolically link each sam into the current working directory
for (i in 1:length(ids)) {
    command = paste0('ln -s ', arioc_sams[i], ' ', ids[i], '.sam')
    run_command(command)
}
