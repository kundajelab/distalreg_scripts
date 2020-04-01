After prior analysis in 2019, we decided to venture into using bed files directly from ENCODE on all tissues/celllines since the performance was comparable across JCandidates (first calling macs peaks using bam) v.s KCandidates (using KundajeLab peak file before running makeCandidateRegions on it)

- Applied to different tissues seen here: https://docs.google.com/presentation/d/1m75d4Kp0vOwGfnXOg1kHn-gAsCRwJ_1UAnu1To_tgZw/edit?usp=sharing
- there were concerns with certain celllines (mainly iPSDF6.9) and iPSDF19.11) which motivated us to just use bam files directly in the macs peak calling algorithm since prior analysis showed comparable results when using Counts. 

In addition, we were also looking into expanding the ABC maps across different tissues/ celllines/ biosamples. 

Code to grab additional celltypes from ENCODE located in experiment
s folder in this repository: https://github.com/kundajelab/distalreg_scripts/tree/master/experiments

Celltypes lookuptable located here: http://mitra.stanford.edu/kundaje/projects/ABC_links/plots/DistributionPlots/newEncodeLookUpTable.txt

Metadata located here: https://docs.google.com/spreadsheets/d/12WsZzuTKpsj0KtxBYEJWWYpSqIbqXMA0d3Sbmu7TZM4/edit?usp=sharing

Upon trying to generate ABC maps across ~70 Celltypes: 
- generated QC Metrics located here : http://mitra.stanford.edu/kundaje/projects/ABC_links/plots/DistributionPlots/SeventyOne_CellTypes_QCMetrics.tsv

Some abnormalities and diagosis can be located here: 

Tissues have different median width of peak and candidate regions as compared to using narrowPeak files diretly from ENCODE: 
Analysis located here: 
https://docs.google.com/document/d/12H4eOWiYSq-88h20Yvr29dUqu2ztZCqc_vV1ZhPXZpU/edit?usp=sharing  
