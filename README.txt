A3B (APOBEC3B) in TCGA-UCEC: exploratory analysis code
Review: "APOBEC3B in Endometrial Carcinoma" (Cancer Pathogenesis and Therapy, CMACPT-D-26-00113R1)

DATA
- cbio_clinical_data_1.xlsx
  Sheet "clinical data": cBioPortal export, study ucec_tcga_pan_can_atlas_2018 (500 primary tumours).
  Sheet "protein exp": APOBEC3B RNA-seq expression (Illumina HiSeq, log2 scale) with a high/low label and PFI,
  as downloaded (UCSC Xena TCGA Hub) on 15 September 2026. Despite the sheet name it is mRNA, not protein.
- Public sources: https://www.cbioportal.org ; https://xenabrowser.net

SAMPLE DERIVATION
189 expression rows -> joined to clinical data by sample ID (tumour samples only) -> kept samples with an assigned
molecular subtype and A3B value = 139 (POLE 14, MSI 35, CN-low 24, CN-high 66). No duplicate patients.
Events in the MSI subgroup: OS 5, PFS 8.

A3B HIGH/LOW
The high/low label is taken from the downloaded file, not computed here. In the file, high >= 7.882 and
low <= 7.861 (log2 units), i.e. a fixed threshold of about 7.87, not the median of the cohort.

RUN ORDER (R, working folder = this folder)
1. step0_excel_to_csv.R             (creates the two CSV files)
2. A3B_UCEC_analysis_subtype.R      (tables, KM curves, Cox models, forest plot)
3. A3B_reviewer_additions.R         (Kruskal-Wallis, Dunn-Holm, continuous HR, interaction tests, adjusted Cox)
   Prints sessionInfo() for package versions.
All analyses are exploratory; no multiplicity correction except Holm for the pairwise subtype comparisons.
