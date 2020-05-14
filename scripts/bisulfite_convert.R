genome_in = 'lambda.fa'
genome_out = 'lambda_bs_artificial.fa'

fa_lines_in = readLines(genome_in)

bs_convert = function(x) {
    if (length(x) == 0 || substr(x, 1, 1) == ">") return(x)
    
    return(gsub("C", "T", x))
}

fa_lines_out = sapply(fa_lines_in, bs_convert)

writeLines(fa_lines_out, con=genome_out)
