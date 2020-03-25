grep "H3K27ac" experiment_report.tsv >  experiment_report_H3K27ac.tsv
# grab celltypes
cut -f7 experiment_report_H3K27ac.tsv > celltypes_H3K27ac.txt
cut -f7 experiment_report_DHS.tsv > celltypes_DHS.txt
# remove whitespace
sed -i '1d' celltypes_H3K27ac.tsv
sed -i '1,2d' celltypes_DHS.tsv
# grab common celltypes:
comm -12 <(sort celltypes_H3K27ac.txt) <(sort celltypes_DHS.txt) | sort -u > common_celltypes.txt


