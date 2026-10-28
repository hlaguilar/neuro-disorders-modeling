library(data.table)


df <- fread(
  "sumStats/MDD_sumStats.tsv",
  skip = "#CHROM" #skipping to #CHROM header
)

#calculate effective sample size for each variant
df[, N := NEFF]


#create the standardized file
clean_mdd <- data.table(
  SNP  = as.character(df$ID),
  A1   = toupper(df$EA),
  A2   = toupper(df$NEA),
  BETA = df$BETA,
  P    = df$PVAL,
  N    = df$N,
  SE   = df$SE,
  INFO = df$IMPINFO,
  MAF  = pmin(
    df$FCON, #MAF for effect allele freq
    1 - df$FCON
  )
)

#remove rows with invalid or missing values
clean_mdd <- clean_mdd[
  !is.na(SNP) & SNP != "" &
    is.finite(BETA) &
    !is.na(P) & is.finite(P) & P >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    is.finite(BETA) &
    is.finite(SE) & SE > 0 &
    is.finite(INFO) & INFO >= 0 & INFO <= 1 &
    is.finite(MAF) & MAF >= 0 & MAF <= 0.5 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]



#save file to clean folder
fwrite(clean_mdd, "clean/mdd_cleaned.txt", sep = "\t")


#print results
cat("Rows in cleaned file:", nrow(clean_mdd), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_mdd$SNP)), "\n")
cat("Invalid BETA: ", sum(!is.finite(df$BETA)), "\n")
cat("Invalid SE: ", sum(!is.finite(df$SE) | df$SE <= 0), "\n")
cat("Invalid P-vals: ", sum(!is.finite(df$PVAL) | df$PVAL < 0 |
                              df$PVAL > 1), "\n")
cat("Invalid INFO: ", sum(!is.finite(df$IMPINFO) | df$IMPINFO < 0 |
                            df$IMPINFO > 1), "\n")
summary(clean_mdd[, .(BETA, P, N, INFO, MAF)])