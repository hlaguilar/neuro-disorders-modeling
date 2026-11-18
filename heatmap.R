install.packages("ggplot2")
library(ggplot2)

S <- as.matrix(LDSCoutput_updated$S)
rownames(S) <- colnames(S)

#convert covariance to genetic correlation
R <- cov2cor(S)

#convert matrix into a format ggplot understands
corr_df <- as.data.frame(as.table(R))

#create heatmap
ggplot(corr_df, aes(x = Var2, y = Var1, fill = Freq)) +
  
  geom_tile(color = "white", linewidth = 0.5) +
  
  #display correlation values
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 3) +
  
  #set color gradient
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Genetic\ncorrelation"
  ) +
  
  #axis formatting
  scale_y_discrete(limits = rev(colnames(R))) +
  scale_x_discrete(position = "top") +
  
  coord_fixed() +
  
  labs(
    title = "Genetic Correlations Across Nine Disorders",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      face = "bold"
    ),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )