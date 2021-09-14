#  Build and push docker images after manually writing dockerfiles for each
#  bundle of software used by the pipeline (note several images have already
#  been created for SPEAQeasy and thus will not be created here).

image_names=("quality_and_trim" "arioc" "filter_alignments" "methyldackel" "bismark" "bioc_kallisto")
image_versions=("0.6.6" "1.43" "1.0" "0.5.2" "0.23.0" "3.13")

num_images=$(($(echo ${image_names[@]} | wc -w) - 1))
for i in $(seq 0 $num_images); do
    docker build \
        -f ${image_names[$i]}_${image_versions[$i]}.dockerfile \
        -t libddocker/${image_names[$i]}:${image_versions[$i]} \
        .
    docker push libddocker/${image_names[$i]}:${image_versions[$i]}
done
