# enhancer_calling

This series of scripts calls enhancers using TSS- and promoter-distal ATAC-seq peaks that overlap enhancers predicted by CRUP.
The reason for this is that CRUP predictions tend to overlap with annotated TSS / promoters and many predictions do not overlap ATAC-seq peaks, making it necessary to add the mentioned filtering.
