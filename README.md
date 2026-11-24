# A3B (APOBEC3B) in TCGA-UCEC: exploratory analysis code

Supports the review "APOBEC3B in Endometrial Carcinoma" (Cancer Pathogenesis and Therapy, manuscript CMACPT-D-26-00113R1).

## Data
- cbio_clinical_data_1.xlsx
  - Sheet "clinical data": cBioPortal export, study ucec_tcga_pan_can_atlas_2018 (500 primary tumours), https://www.cbioportal.org
  - Sheet "protein exp": APOBEC3B mRNA expression, log2(normalised RSEM count + 1), Illumina HiSeq, from the UCSC Xena TCGA Hub (dataset ID TCGA.UCEC.sampleMap/HiSeqV2, version 2017-10-13), https://xenabrowser.net, downloaded on 15 September 2026. Despite the sheet name it contains mRNA, not protein. It carries the high/low label supplied with the download.

## Sample derivation
189 Xena samples (175 primary tumours, 13 solid-tissue normals, 1 other) -> 175 primary tumours -> 145 matched to cBioPortal by sample barcode -> 6 without a molecular subtype excluded -> 139 tumours (POLE 14, MSI 35, CN-low 24, CN-high 66). One sample per patient. MSI subgroup events: OS 5, PFS 8.

## A3B high/low
The label is taken from the downloaded file and is not computed here. High >= 7.882 and low <= 7.861 (log2 units); 74 high, 65 low. It is not the median of the cohort.

## Run order (R; working folder = this folder)
1. step0_excel_to_csv.R - creates the two CSV files
2. A3B_UCEC_analysis_subtype.R - tables, Kaplan-Meier curves, Cox models, forest plot
3. A3B_reviewer_additions.R - Kruskal-Wallis, Dunn-Holm, continuous hazard ratios, interaction tests, adjusted Cox model; prints sessionInfo() for package versions

## Note
All analyses are exploratory. No multiplicity correction was applied to the survival analyses; Holm adjustment was used for the pairwise subtype comparisons.
