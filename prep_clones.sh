for i in {11..118};
do
python /Users/Nick/conga/scripts/setup_10x_for_conga.py \
--filtered_contig_annotations_csvfile "/Users/Nick/Documents/GitHub/COVID_TCR/data/sequencingRuns/s""$i""/filtered_contig_annotations.csv" \
--output_clones_file "/Users/Nick/Documents/GitHub/COVID_TCR/data/preCoNGA/contigs/s""$i""_clones.tsv" \
--organism human \
--no_kpca
done