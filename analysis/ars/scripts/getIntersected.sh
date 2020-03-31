while read p;
do
	echo $p
	a=($p)
	bedtools intersect -a EnhancerPreds_K562_nonthres.tsv -b /mnt/lab_data3/kmualim/PredictionFiles_qnorm/Predictions_${a[0]}_qnorm/EnhancerPredictionsAllPutative.txt.gz -f 0.5  -wa -wb  | gzip -c > non_thres/${a[0]}_Intersected.tsv.gz
done < $1 
