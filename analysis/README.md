# 03/25/2020


# 2019 Analysis using different peak files and different methods 
Link to the slides that compared barplots across different peak files:
Barplots across different peak files + different methods to calculate “activity”: https://docs.google.com/presentation/d/1x9cZqa5ObNOUM6gy8wZBE_yDNBFie-x1a4FyqKNkHH8/edit?usp=sharing

File specifications + process for making barplots + PR curves: https://docs.google.com/document/d/1NStb6erFkJUQdU3wIaXt0g7Q7CJ_k7zNvFvUuK8VX5M/edit?usp=sharing

Files:
DHS Stam file: https://www.encodeproject.org/experiments/ENCSR000EOT/ (specifically used ENCFF821KDJ (bed narrowPeak)) — lifted these peaks over to hg19 genome build
KundajeLab DHS peak file used:
(ENCODE file: https://www.encodeproject.org/experiments/ENCSR000EOY/) (ENCFF801RMG)
on mitra: http://mitra.stanford.edu/kundaje/projects/atlas/dnase_processed/k562_dnase/hg19/atac/5873d283-2b12-4692-9266-bb28329d83a1/call-macs2/shard-0/execution/ENCFF801RMG.merged.nodup.pval0.01.300K.bfilt.narrowPeak.gz

The above was used for Kundaje Lab peak file (labelled ‘KCandidate’) which, I assume, are the raw peaks. Below should be the overlap peaks located here: http://mitra.stanford.edu/kundaje/projects/atlas/dnase_processed/k562_dnase/hg19/atac/5873d283-2b12-4692-9266-bb28329d83a1/call-overlap/shard-0/execution/rep1_rep2.overlap.bfilt.narrowPeak.gz
