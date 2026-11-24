## Step 0: convert the cBioPortal Excel file into the two CSV files the analysis scripts expect.
## Put cbio_clinical_data_1.xlsx in your working folder (Session > Set Working Directory) and run this first.
if (!"readxl" %in% installed.packages()[, "Package"]) install.packages("readxl")
library(readxl)
f <- "cbio_clinical_data_1.xlsx"
write.csv(read_excel(f, sheet = "clinical data"), "cbio_clinical_data.csv",   row.names = FALSE)
write.csv(read_excel(f, sheet = "protein exp"),   "cbio_clinical_data-2.csv", row.names = FALSE)
cat("Done: cbio_clinical_data.csv and cbio_clinical_data-2.csv created\n")
