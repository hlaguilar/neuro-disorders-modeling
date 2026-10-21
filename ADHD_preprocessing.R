library(data.table)

df <- fread(
  "sumStats/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta",
  nrows = 5
)

names(df)
head(df)


#calculate ADHD effective sample size
df[, N_eff := 4 / (1 / Nca + 1 / Nco)]
 
df[, .(SNP, OR, P, Nca, Nco, N_eff)][1:5]


library(data.table)

# 1. Read the original ADHD summary statistics
df <- fread(
  "sumStats/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta"
)

# 2. Calculate effective sample size for each variant
df[, N := 4 / (1 / Nca + 1 / Nco)]

# 3. Create the standardized file
clean_adhd <- data.table(
  SNP  = as.character(df$SNP),
  A1   = toupper(df$A1),
  A2   = toupper(df$A2),
  BETA = log(df$OR),
  P    = df$P,
  N    = df$N,
  INFO = df$INFO,
  
  MAF  = pmin(
    df$FRQ_U_186843,
    1 - df$FRQ_U_186843
  )
)

# 4. Remove rows with invalid or missing values
clean_adhd <- clean_adhd[
  !is.na(SNP) &
    is.finite(BETA) &
    !is.na(P) & P >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]

# 5. Create a folder for standardized files
dir.create("clean", showWarnings = FALSE)

# 6. Save the standardized ADHD file
fwrite(clean_adhd, "clean/ADHD.txt", sep = "\t")

# 7. Check the result
head(clean_adhd)
names(clean_adhd)
nrow(clean_adhd)
