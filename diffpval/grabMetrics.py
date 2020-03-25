#! bin/python3
import sys
from metrics import *

if __name__=="__main__":
    files = sys.argv[1]
    cells = pd.read_csv(files, sep="\t", header=None)
    pvalues = ["01", '005', "05", "1"]
    
    for pval in pvalues:
        for cell in cells[0]:
            macs_peaks = "/srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_{}/Peaks_{}_{}/NA_peaks.narrowPeak".format(str(pval), str(i), str(pval))
            genome_tss = "/users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"
            peak_outdir = "/srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_{}/Peaks_{}_{}".format(str(pval), str(i), str(pval))
            grab_nearest_tss_from_peak(macs_peaks, genome_tss, peak_outdir)
            PeakFileQC(macs_peaks, peak_outdir)
