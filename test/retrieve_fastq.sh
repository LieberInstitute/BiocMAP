#  This script documents how FASTQ files were retrieved to ultimately form
#  the test data for this repository.

num_reads=100000

module load sratoolkit/2.10.8
module load conda_R/4.1.x

for pairing in (single paired); do
    for species in (human mouse); do
        if [[ $species == "human" ]] && [[ $pairing == "paired" ]]; then
            #   "AgeNeunSortedWGBS" samples are available from synapse
            Rscript retrieve_fastq_human_paired.R
        else
            #  Prefetch
            cat SRR_Acc_List_${pairing}_${species}.txt | xargs -I{} prefetch {}
            
            if [[ $pairing == "paired" ]]; then
                cat SRR_Acc_List.txt | xargs -I{} fastq-dump --split-files --gzip {}
                
                for i in $(seq 1 $(cat SRR_Acc_List_${pairing}_${species}.txt | wc -l)); do
                    sra_id=$(awk "NR==$i" SRR_Acc_List_${pairing}_${species}.txt)
                    
                    #  Subset
                    gunzip -c ${sra_id}_1.fastq.gz | head -n $num_reads | gzip > $species/$pairing/fastq/${sra_id}_1.fastq.gz
                    gunzip -c ${sra_id}_2.fastq.gz | head -n $num_reads | gzip > $species/$pairing/fastq/${sra_id}_2.fastq.gz
                done
            else
                cat SRR_Acc_List.txt | xargs -I{} fastq-dump --gzip {}
                
                for i in $(seq 1 $(cat SRR_Acc_List_${pairing}_${species}.txt | wc -l)); do
                    sra_id=$(awk "NR==$i" SRR_Acc_List_${pairing}_${species}.txt)
                    
                    #  Subset
                    gunzip -c ${sra_id}.fastq.gz | head -n $num_reads | gzip > $species/$pairing/fastq/${sra_id}.fastq.gz
                done
            fi
        fi
    done
done
