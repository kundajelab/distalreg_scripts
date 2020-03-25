#! bin/bash 

while read p;
do
	a=($p)
	echo ${a[0]} ${a[1]}
	macs2 callpeak -f BAM -g hs -p .05 --call-summits --outdir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_05/Peaks_${a[0]} -t $2/${a[1]}
	macs2 callpeak -f BAM -g hs -p .01 --call-summits --outdir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_01/Peaks_${a[0]} -t $2/${a[1]}
	macs2 callpeak -f BAM -g hs -p .005 --call-summits --outdir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_005/Peaks_${a[0]} -t $2/${a[1]}
	macs2 callpeak -f BAM -g hs -p .1 --call-summits --outdir /mnt/lab_data3/kmualim/PeakAndNeighborhoods/Peaks_${a[0]} -t $2/${a[1]}
done < $1  

conda deactivate
conda activate macs-py2.7

while read p;
do
       	a=($p)	
	python /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_05/Peaks_${a[0]}/NA_peaks.narrowPeak --chrom_sizes /mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes --regions_blacklist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed --regions_whitelist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed --peakExtendFromSummit 250 --nStrongestPeaks 150000 --outDir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_05/Peaks_${a[0]}/ --bam /srv/scratch/kmualim/ABC_data/ENCODEdata/cellline_files/${a[1]} --genome_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed
	
	python /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_01/Peaks_${a[0]}/NA_peaks.narrowPeak --chrom_sizes /mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes --regions_blacklist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed --regions_whitelist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed --peakExtendFromSummit 250 --nStrongestPeaks 150000 --outDir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_01/Peaks_${a[0]}/ --bam /srv/scratch/kmualim/ABC_data/ENCODEdata/cellline_files/${a[1]} --genome_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed
	
	python /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_005/Peaks_${a[0]}/NA_peaks.narrowPeak --chrom_sizes /mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes --regions_blacklist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed --regions_whitelist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed --peakExtendFromSummit 250 --nStrongestPeaks 150000 --outDir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_005/Peaks_${a[0]}/ --bam /srv/scratch/kmualim/ABC_data/ENCODEdata/cellline_files/${a[1]} --genome_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed
	python /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_1/Peaks_${a[0]}/NA_peaks.narrowPeak --chrom_sizes /mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes --regions_blacklist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed --regions_whitelist /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed --peakExtendFromSummit 250 --nStrongestPeaks 150000 --outDir /srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_1/Peaks_${a[0]}/ --bam /srv/scratch/kmualim/ABC_data/ENCODEdata/cellline_files/${a[1]} --genome_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed
done < $1 

