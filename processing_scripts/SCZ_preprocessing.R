library(data.table)

df <- fread(
  "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",
  skip = "CHROM" #skipping to #CHROM header
)

#calculate effective sample size for each variant
df[, N := NEFF]


#create the standardized file
clean_scz <- data.table(
  SNP  = as.character(df$ID),
  A1   = toupper(df$A1),
  A2   = toupper(df$A2),
  BETA = df$BETA,
  P    = df$PVAL,
  N    = df$N,
  INFO = df$IMPINFO,
  MAF  = pmin(
    df$FCON, #MAF for control allele freq
    1 - df$FCON
  )
)

#remove rows with invalid or missing values
clean_scz <- clean_scz[
  !is.na(SNP) & SNP != "" &
    is.finite(BETA) &
    !is.na(P) & is.finite(P) & P >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    is.finite(INFO) & INFO >= 0 & INFO <= 1 &
    is.finite(MAF) & MAF >= 0 & MAF <= 0.5 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]



#save file to clean folder
fwrite(clean_scz, "clean/SCZ_cleaned.txt", sep = "\t")


#print results
cat("Rows in cleaned file:", nrow(clean_scz), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_scz$SNP)), "\n")
cat("Invalid BETA: ", sum(!is.finite(df$PVAL) | df$PVAL < 0 | 
                            df$PVAL > 1), "\n")
cat("Invalid P-vals: ", sum(!is.finite(df$PVAL) | df$PVAL < 0 |
                              df$PVAL > 1), "\n")
cat("Invalid INFO: ", sum(!is.finite(df$IMPINFO) | df$IMPINFO < 0 |
                            df$IMPINFO > 1), "\n")
cat("Invalid OR:", sum(!is.finite(df$OR) | df$OR <= 0), "\n")
summary(clean_scz[, .(BETA, P, N, INFO, MAF)])
