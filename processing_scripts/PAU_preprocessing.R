library(data.table)


df <- fread(
  "sumStats/PAU_EUR_AUG2023.txt",
)


#calculate effective sample size for each variant
df[, N := SampleSize]

#create the standardized file
clean_pau <- data.table(
  SNP  = as.character(df$SNP_ID),
  A1   = toupper(df$EA),
  A2   = toupper(df$Allele2),
  BETA = df$Effect,  
  P    = df$PValue,
  N    = df$SampleSize
)


#remove rows with invalid or missing values
clean_pau <- clean_pau[
  !is.na(SNP) & SNP != "" &
    !is.na(P) & is.finite(P) & P >= 0 & P <= 1 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]




#save file to clean folder
fwrite(clean_pau, "clean/pau_cleaned.txt", sep = "\t")

#print results
cat("Rows in cleaned file:", nrow(clean_pau), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_pau$SNP)), "\n")



