
library(data.table)


#load raw CUD summary statistics
df <- fread(
  "sumStats/dbGAP_EUR_CanMVP_META"
)


#create the standardized file
clean_cud <- data.table(
  SNP = as.character(df$SNP),
  A1  = toupper(df$A1),
  A2  = toupper(df$A2),
  Z   = df$BETA / df$SE,
  P   = df$P.value
)


#remove rows with invalid or missing values
clean_cud <- clean_cud[
  !is.na(SNP) & SNP != "" &
    !is.na(Z) & is.finite(Z) &
    !is.na(P) & is.finite(P) & P >= 0 & P <= 1 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]


#save file to clean folder
fwrite(clean_cud, "clean/cud_cleaned.txt", sep = "\t")


#print results
cat("Rows in cleaned file:", nrow(clean_cud), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_cud$SNP)), "\n")
summary(clean_cud[, .(P, Z)])
