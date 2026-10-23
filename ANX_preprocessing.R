library(data.table)

df <- fread(
  "sumStats/ANX_2026_daner_fullANX_v12_woUTAH_11022026",
  nrows = 5
)


#calculate ANX effective sample size
df[, N_eff := 4 / (1 / Nca + 1 / Nco)]

df[, .(SNP, OR, P, Nca, Nco, N_eff)][1:5]


library(data.table)

df <- fread(
  "sumStats/ANX_2026_daner_fullANX_v12_woUTAH_11022026",
)


# Calculate effective sample size for each variant
df[, N := 4 / (1 / Nca + 1 / Nco)]

# Create the standardized file
clean_anx <- data.table(
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

# Remove rows with invalid or missing values
clean_anx <- clean_anx[
  !is.na(SNP) &
    is.finite(BETA) &
    !is.na(P) & P >= 0 & P <= 1 &
    is.finite(N) & N > 0 &
    A1 %chin% c("A", "C", "G", "T") &
    A2 %chin% c("A", "C", "G", "T") &
    A1 != A2
]


#save file to clean folder
fwrite(clean_anx, "clean/ANX.txt", sep = "\t")
