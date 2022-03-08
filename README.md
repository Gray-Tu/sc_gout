# sc_gout
scripts of article

## cellranger

we analysis 10x Genomics data with cellranger v3.1.0 and 6.0.2

| SampleID | NovaSeq-batch | GROUP          | CellRangerVersion |
| -------- | ------------- | -------------- | ----------------- |
| Case01   | 20200609      | Acute          | 3.1.0             |
| Case02   | 20200609      | Acute          | 3.1.0             |
| Case03   | 20200609      | Acute          | 3.1.0             |
| Case04   | 20200609      | Acute          | 3.1.0             |
| Case05   | 20200609      | Acute          | 3.1.0             |
| Case06   | 20200609      | Acute          | 3.1.0             |
| Case07   | 20200609      | Acute          | 3.1.0             |
| YB       | 20200609      | Health         | 3.1.0             |
| N01F     | 20200915      | Health         | 3.1.0             |
| N05S     | 20201112      | Health         | 3.1.0             |
| Gout146  | 20210708      | Synovial fluid | 6.0.2             |
| MH       | 20211103      | Health         | 6.0.2             |
| RG       | 20211103      | Health         | 6.0.2             |
| Gout151  | 20211103      | Synovial fluid | 6.0.2             |
| Gout183  | 20211103      | Synovial fluid | 6.0.2             |

## cell qc

we remain cells follow:

PBMC: (pc30m2h5_cell.qc.R, Gout_10x_detect_doublets.ipynb, and N01F_N05S_10x_detect_doublets.ipynb)

SF: 

1. nFeature_RNA > 500
2. nFeature_RNA < 4500
3. percent.mt < 15 
4. nCount_RNA < 30000
5. not doublets

## SCT and rPCA integration (test PC30 and selecting PC50)

PBMC: pc30m2h5_cell.qc.R

SF: SF_1_cell_qc.R

selecting health controls as reference 

## Annotation main cell type by SingleR

1_PBMC_SYN_re_annot.R

## Abundance test

0_paper_Fig01.R

A0_paper_Fig02-alt.R

## subclustering

3_MainDEG.R

3_add_subclustering_DEG_Func.R

A0_paper_Fig03.R

## IS and ES

A0_paper_Fig05_SCORE.R

## CCI

0_paper_Fig06_CCI_PBMCacute.R

0_paper_FigCCI_combinePlot.R
