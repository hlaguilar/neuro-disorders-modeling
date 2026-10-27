library(data.table)

df <- fread(
  "sumStats/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta"
  )

names(df)
head(df)



#calculate effective sample size for each variant
df[, N := 4 / (1 / Nca + 1 / Nco)]

#create the clean file
clean_adhd <- data.table(
  SNP  = as.character(df$SNP),
  A1   = toupper(df$A1),
  A2   = toupper(df$A2),
  BETA = log(df$OR),
  P    = df$P,
  N    = df$N,
  INFO = df$INFO,
  MAF  = pmin(
    #A1 allele freqs
    df$FRQ_U_186843,
    1 - df$FRQ_U_186843
  )
)

#remove rows with invalid or missing values
clean_adhd <- clean_adhd[
  !is.na(SNP) & SNP != "" &
    is.finite(BETA) &
    !is.na(P) & is.finite(P) & P >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]

#create a folder for standardized files
dir.create("clean", showWarnings = FALSE)

#save the standardized ADHD file
fwrite(clean_adhd, "clean/ADHD_cleaned.txt", sep = "\t")

#print results
cat("Rows in cleaned file:", nrow(clean_adhd), "\n")
cat("Duplicate SNP IDs:", sum(duplicated(clean_adhd$SNP)), "\n")
cat("Invalid OR:", sum(!is.finite(df$OR) | df$OR <= 0), "\n")
summary(clean_adhd[, .(BETA, P, N, INFO, MAF)])
head(clean_adhd)
nrow(clean_adhd)
