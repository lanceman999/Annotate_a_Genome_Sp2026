library(ggplot2)
library(dplyr)

# ============================================
# Load and clean data
# ============================================

# Load basic association results
assoc <- read.table("basic_assoc.assoc", header = TRUE, stringsAsFactors = FALSE)
head(assoc)

# Load logistic regression results
logistic <- read.table("logistic_assoc.assoc.logistic", header = TRUE, stringsAsFactors = FALSE)
head(logistic)

# Filter logistic to SNP results only (TEST == "ADD")
logistic_clean <- logistic %>%
  filter(TEST == "ADD") %>%
  filter(!is.na(P))

# ============================================
# Prepare data for Manhattan plot
# ============================================

# For basic association
assoc_clean <- assoc %>%
  filter(!is.na(P)) %>%
  arrange(CHR, BP)

# Calculate cumulative position for x-axis
assoc_manhattan <- assoc_clean %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(assoc_clean, ., by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BP_cum = BP + tot)

# Get chromosome center positions for x-axis labels
axis_df <- assoc_manhattan %>%
  group_by(CHR) %>%
  summarise(center = (max(BP_cum) + min(BP_cum)) / 2)

# ============================================
# Manhattan Plot - Basic Association
# ============================================

manhattan_basic <- ggplot(assoc_manhattan, aes(x = BP_cum, y = -log10(P), color = as.factor(CHR))) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = rep(c("#1B9E77", "#666666"), 11)) +
  scale_x_continuous(labels = axis_df$CHR, breaks = axis_df$center, expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(assoc_manhattan$P)) + 1)) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P)),
    title = "Manhattan Plot - Basic Association (Unadjusted)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 16, face = "bold")
  )

manhattan_basic

# ============================================
# Manhattan Plot - Logistic Regression (Adjusted)
# ============================================
# Prepare logistic data
logistic_manhattan <- logistic_clean %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(logistic_clean, ., by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BP_cum = BP + tot)

axis_df_log <- logistic_manhattan %>%
  group_by(CHR) %>%
  summarise(center = (max(BP_cum) + min(BP_cum)) / 2)

manhattan_logistic <- ggplot(logistic_manhattan, aes(x = BP_cum, y = -log10(P), color = as.factor(CHR))) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = rep(c("#D95F02", "#666666"), 11)) +
  scale_x_continuous(labels = axis_df_log$CHR, breaks = axis_df_log$center, expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(logistic_manhattan$P)) + 1)) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P)),
    title = "Manhattan Plot - Logistic Regression (MDS Adjusted)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 16, face = "bold")
  )

manhattan_logistic

# ============================================
# QQ Plot - Basic Association
# ============================================
# Calculate expected vs observed p-values
qq_basic <- assoc_clean %>%
  arrange(P) %>%
  mutate(
    expected = -log10(ppoints(n())),
    observed = -log10(P)
  )

# Calculate genomic inflation factor (lambda)
chisq <- qchisq(1 - assoc_clean$P, 1)
lambda_basic <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

qq_plot_basic <- ggplot(qq_basic, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.5, size = 1, color = "#1B9E77") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P)),
    title = "QQ Plot - Basic Association (Unadjusted)",
    subtitle = paste0("λ = ", round(lambda_basic, 3))
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  coord_fixed(ratio = 1)

qq_plot_basic

# ============================================
# QQ Plot - Logistic Regression (Adjusted)
# ============================================
qq_logistic <- logistic_clean %>%
  arrange(P) %>%
  mutate(
    expected = -log10(ppoints(n())),
    observed = -log10(P)
  )

# Calculate lambda for logistic
chisq_log <- qchisq(1 - logistic_clean$P, 1)
lambda_logistic <- median(chisq_log, na.rm = TRUE) / qchisq(0.5, 1)

qq_plot_logistic <- ggplot(qq_logistic, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.5, size = 1, color = "#D95F02") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = 1, linewidth = 1) +
  labs(
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P)),
    title = "QQ Plot - Logistic Regression (MDS Adjusted)",
    subtitle = paste0("λ = ", round(lambda_logistic, 3))
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  coord_fixed(ratio = 1)

qq_plot_logistic

# ============================================
# Print lambda values
# ============================================

cat("Genomic inflation factor (lambda):\n")
cat("Basic association:", round(lambda_basic, 3), "\n")
cat("Logistic regression:", round(lambda_logistic, 3), "\n")