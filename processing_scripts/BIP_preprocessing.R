library(data.table)

df <- fread(
  "sumStats/pgc-bip2021-all.vcf.tsv",
)

#calculate effective sample size for each variant
df[, N := 4 / (1 / Nca + 1 / Nco)]

#create the standardized file
clean_bip <- data.table(
  SNP  = as.character(df$SNP),
  A1   = toupper(df$A1),
  A2   = toupper(df$A2),
  BETA = log(df$OR),
  P    = df$P,
  N    = df$N,
  INFO = df$INFO,
  MAF  = pmin(
    df$FRQ_U_729602, #MAF for control allele freq
    1 - df$FRQ_U_729602
  )
)


#remove rows with invalid or missing values
clean_bip <- clean_anx[
  !is.na(SNP) & SNP != "" &
    is.finite(BETA) &
    !is.na(P) & is.finite(P) >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    is.finite(INFO) & INFO >= 0 & INFO <= 1 &
    is.finite(MAF) & MAF >= 0 & MAF <= 0.5 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]



#save file to clean folder
fwrite(clean_bip, "clean/BIP_cleaned.txt", sep = "\t")


#print results
cat("Rows in cleaned file:", nrow(clean_bip), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_bip$SNP)), "\n")
cat("Invalid OR:", sum(!is.finite(df$OR) | df$OR <= 0), "\n")
summary(clean_bip[, .(BETA, P, N, INFO, MAF)])
head(clean_bip)
nrow(clean_bip)