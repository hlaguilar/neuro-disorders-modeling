library(data.table)

df <- fread(
  "sumStats/ANX_2026_daner_fullANX_v12_woUTAH_11022026",
)

#calculate effective sample size for each variant
df[, N := 4 / (1 / Nca + 1 / Nco)]

#create the standardized file
clean_anx <- data.table(
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
clean_anx <- clean_anx[
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
fwrite(clean_anx, "clean/ANX_cleaned.txt", sep = "\t")


#print results
cat("Rows in cleaned file:", nrow(clean_anx), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_anx$SNP)), "\n")
cat("Invalid OR:", sum(!is.finite(df$OR) | df$OR <= 0), "\n")
summary(clean_anx[, .(BETA, P, N, INFO, MAF)])
head(clean_anx)
nrow(clean_anx)
