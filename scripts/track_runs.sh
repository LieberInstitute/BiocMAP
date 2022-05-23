#   Script to log metadata about BiocMAP runs submitted at JHPCE. This helps
#   decentralize knowledge about which WGBS datasets have been processed,
#   where they were processed, and other info

log_dir='/dcs04/lieber/lcolladotor/libdData_LIBD001/processed-data/nextflow_logs/BiocMAP'
BiocMAP_log=$1
half=$2

#   First, only continue if the BiocMAP log indicates completion of the
#   pipeline
if [[ $(grep -E "^Completed at: " $BiocMAP_log | wc -l) -ne 1 ]]; then
    exit 0
fi

input_dir=$(
    grep -E "^Input dir *: /.*" $BiocMAP_log |
    tail -n 1 | cut -d ":" -f 2 | tr -d " "
)

#   Test runs should not be logged
is_small_test=$(
    echo $input_dir |
    grep -E ".*/BiocMAP/test/(human|mouse)/(paired|single)$" | wc -l
)
if [[ $is_small_test -eq 1 ]]; then
    exit 0
fi

#   Get the output directory for the run
out_dir=$(
    grep -E "^Output dir *: /.*$" $BiocMAP_log |
    tail -n 1 | cut -d ":" -f 2 | tr -d " "
)

#   Determine the manifest path, then count the number of samples
if [[ $half == "first" ]]; then
    man_path=$input_dir/samples.manifest
elif [[ $half == "second" ]]; then
    man_path=$(
        grep -E "^manifest *= *" $input_dir/rules.txt |
        cut -d "=" -f 2 | tr -d " "
    )
else
    echo "While tracking run: received an invalid 'half' argument '$half'. Should be 'first' or 'second'."
    exit 1
fi

pairedness=$(
    grep -E "^Sample *: (single|paired)$" $BiocMAP_log |
    tail -n 1 | cut -d ":" -f 2 | tr -d " "
)
if [[ "$pairedness" == "single" ]]; then
    num_samples=$(cut -f 3 $man_path | sort -u | wc -l)
elif [[ "$pairedness" == "paired" ]]; then
    num_samples=$(cut -f 5 $man_path | sort -u | wc -l)
else
    echo "While tracking run: could not determine if samples were paired-end."
    exit 1
fi

#   Tab-separated list:
#   [output dir] [date] [number of samples] [user] [module ("first" or "second")]

log_path=$(mktemp -p $log_dir -t run_XXXXXXX.log)
echo -e "${out_dir}\t$(date +%Y-%m-%d,%H:%M)\t${num_samples}\t$(whoami)\t$half" >
    $log_path
chmod 755 $log_path
