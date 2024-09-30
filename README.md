# Scripts to run certain analyses on McCutcheon et al. 2023 Nat. Genetics

## DESeq2 analysis to detect significant differences in gRNA and gene expression abundances between two conditions

Assuming a `path/to/counts_table.txt` counts table file like this:
```csv
gRNA,gRNA_sequence,Low_R1,Low_R2,Low_R3,High_R1,High_R2,High_R3
IL2RA_gRNA1,CCTTGTTTCAAATGGATTTTC,1723,819,4081,2242,2216,2158
IL2RA_gRNA2,TCGAGAAAATCCATTTGAAAC,1456,851,3703,1552,2262,2197
IL2RA_gRNA3,CGAGAAAATCCATTTGAAACA,902,1226,2840,1475,1859,2178
IL2RA_gRNA4,AGAAGTAGTAATGTTCTAAAA,1387,1283,2903,1484,2856,2149
IL2RA_gRNA5,TTTCTGTAAAGTTGCACTTGT,1328,1241,3314,2207,2940,2143
...
```

And a `path/to/metadata.txt` metadata file like this:
```csv
Sample,id,sorted_bin,replicate
1,Low R1,Low,1
2,Low R2,Low,2
3,Low R3,Low,3
4,High R1,High,1
5,High R2,High,2
6,High R3,High,3
```

You can run DESeq2 like this:
```sh
Rscript run_deseq2.R \
    path/to/counts_table.txt \
    path/to/metadata.txt \
    path/to/output.tsv \
    --skip_columns gRNA_sequence
```
The `output.tsv` file will contain the DESeq2 results.

