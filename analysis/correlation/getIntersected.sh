#! bin/bash 

## $1 fed in should contain list of celltypes 
## $2 should be the EnhancerPredictionsAllPutative file of the reference celltype
## $3 should be the outdir for intersected files
while read p;
do
	echo $p
	a=($p)
	bedtools intersect -a $2 -b /mnt/lab_data3/kmualim/PredictionFiles_qnorm/Predictions_${a[0]}_qnorm/EnhancerPredictionsAllPutative.txt.gz -f 0.5  -wa -wb  | gzip -c > $3/${a[0]}_Intersected.tsv.gz
done < $1 
