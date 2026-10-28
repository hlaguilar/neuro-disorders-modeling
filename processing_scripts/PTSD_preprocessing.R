library(data.table)


df <- fread(
  "sumStats/ptsd_sumStats.vcf",
  skip = "#CHROM" #skipping to #CHROM header
)

#calculate effective sample size for each variant
df[, N := NEFF]


#create the standardized file
clean_ptsd <- data.table(
  SNP  = as.character(df$ID),
  A1   = toupper(df$A1),
  A2   = toupper(df$A2),
  Z    = df$Z,
  P    = df$P,
  N    = df$N,
  INFO = df$IMPINFO,
  MAF  = pmin(
    df$FREQ, #MAF for control allele freq
    1 - df$FREQ
  )
)


#remove rows with invalid or missing values
clean_ptsd <- clean_ptsd[
  !is.na(SNP) & SNP != "" &
    is.finite(Z) &
    !is.na(P) & is.finite(P) & P >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    is.finite(MAF) & MAF >= 0 & MAF <= 0.5 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]


#save file to clean folder
fwrite(clean_ptsd, "clean/ptsd_cleaned.txt", sep = "\t")


#print results
cat("Rows in cleaned file:", nrow(clean_ptsd), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_ptsd$SNP)), "\n")
cat("Invalid P-vals: ", sum(!is.finite(df$PVAL) | df$PVAL < 0 |
                              df$PVAL > 1), "\n")
cat("Invalid OR:", sum(!is.finite(df$OR) | df$OR <= 0), "\n")
summary(clean_ptsd[, .(Z, P, N, INFO, MAF)])
