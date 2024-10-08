# Scripts to run certain analyses on McCutcheon et al. 2023 Nat. Genetics

## [alignment-22bpSE-no-reverse-comp.sh](alignment-22bpSE-no-reverse-comp.sh)
Alignment shell script to run `Bowtie2` and compute counts with `samtools`.  

## [run_deseq2.R](run_deseq2.R)
DESeq2 analysis to detect significant differences in gRNA and gene expression abundances between two conditions. A combined table of raw counts is expected.

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


## [preprocess_sc_rna_seq_screen.R](preprocess_sc_rna_seq_screen.R)
Preprocess for analsyis scRNA-seq CRISPR screen 10X data computed with `cellranger`. This script used `Seurat::Read10X` function to load data in folder containing 3 files: `features.tsv.gz`, `barcodes.tsv.gz` and `matrix.mtx.gz`. It produces a QC'ed `.rds` file that can be used in [single_job_mast_sc_rna_seq_screen.R](single_job_mast_sc_rna_seq_screen.R) to find perturbed/DE genes.

Run like this:
```sh
Rscript process_scRNAseq.R \
    --input_dir `/path/to/input` \
    --output_dir `/path/to/output` \
    --donor `<DONOR_NAME>` \
    --crispr_app `<CRISPRi/a>`
```


## [single_job_mast_sc_rna_seq_screen.R](single_job_mast_sc_rna_seq_screen.R) 
Script to find perturbed/differentially expressed genes using the `FindMarkers` function in [Seurat](https://satijalab.org/seurat/) and [MAST](https://rglab.github.io/MAST/) to detect significant changes in transcription levels induced by a single perturbation.

For example, to detect significant changes induce by the 3rd gRNA in the set, run this:
```sh
Rscript single_job_mast_sc_rna_seq_screen.R \
    --job 3 \
    --rds-files `/path/to/CRISPRa_D1_qc.rds /path/to/CRISPRa_D2_qc.rds /path/to/CRISPRa_D3_qc.rds` \
    --output-basename /path/to/output_basename \
    --min-umi-thres 4
```