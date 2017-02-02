# long read pipeline.

Uses a reference genome assembly to call variants and leverages reads where possible.



# Input file format for pilon running (tab or comma-delimited).

```
    SAMPLE_NAME SCAFFOLDS PE SE
    SM1 SM.fasta P1*.gz:P2*.gz,P1*.gz:P2*.gz fastq.gz
```

reference is specified on the command-line. 
