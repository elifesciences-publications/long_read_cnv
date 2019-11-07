### This code is associated with the paper from Bloom et al., "Rare variants contribute disproportionately to quantitative trait variation in yeast". eLife, 2019. http://dx.doi.org/10.7554/eLife.49212
# Long-read hybrid CNV caller (LRHCNV)

LRHCNV utilises long-reads and short-reads to call structural variants. 

The long reads are used to uncover a set of putative CNVs, and then the short reads are used to provide support for the extracted breakpoints.

Simplied version of the pipeline utilised a moleculo genome assembly is shown below.

![alt tag](https://raw.githubusercontent.com/theboocock/long_read_cnv/master/img/LRCNV_pipeline.png)
 
# Input file format for pilon running (tab or comma-delimited).

```
    SAMPLE_NAME SCAFFOLDS PE SE
    SM1 SM.fasta P1*.gz:P2*.gz,P1*.gz:P2*.gz fastq.gz
```

