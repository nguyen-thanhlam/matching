library(ggplot2)
library(gridExtra)
library(dplyr)

col = c('strategy',
        'algo',
        'approach',
        'n',
        'n.cova',
        'event.prob',
        'ttm.prob',
        'seed',
        'mh_pval',
        'clog_coef',
        'clog_pval',
        'clog_lower',
        'clog_upper')
# Ensure data is available
beta = 1.5
resa1 = read.csv("match_resA1.csv")
resa2 = read.csv("match_resA2.csv")
# Convert to data frame if needed
res <- rbind(resa1[,c(col)],resa2[,c(col)])

# Treat variables as categorical
res$strategy <- factor(res$strategy, levels = c("matching", "counter-matching"), labels = c("m", "cm"))
res$algo <- factor(res$algo, levels = c("A1", "A2"), labels = c("A1", "A2"))
res$approach <- factor(res$approach)
res$n.cova <- factor(res$n.cova)
res$event.prob <- factor(res$event.prob)
res$ttm.prob <- factor(res$ttm.prob)


# Combine strategy and approach for plotting, and order for bars

# Order for bars: m-drs, cm-drs, m-ps, cm-ps, m-mhl, cm-mhl
approach_levels <- c("drs", "ps", "mhl")
strat_app_levels <- c("m-drs", "cm-drs", "m-ps", "cm-ps", "m-mhl", "cm-mhl")
res$strat_app <- factor(paste0(res$strategy, "-", res$approach), levels = strat_app_levels)

# Function to calculate type 1 error with confidence interval
calc_type1_error <- function(data, alpha = 0.05) {
  n_total <- sum(!is.na(data))
  n_reject <- sum(data < alpha, na.rm = TRUE)
  prop <- n_reject / n_total
  
  # Wilson score confidence interval
  z <- qnorm(1 - 0.05/2)
  denom <- 1 + z^2/n_total
  center <- (prop + z^2/(2*n_total)) / denom
  margin <- z * sqrt((prop*(1-prop) + z^2/(4*n_total))/n_total) / denom
  
  return(data.frame(
    prop = prop,
    lower = max(0, center - margin),
    upper = min(1, center + margin),
    n = n_total
  ))
}


# Function to calculate coverage (proportion of true value 0 between lower and upper)
calc_coverage <- function(data, beta) {
  n_total <- sum(!is.na(data$clog_lower) & !is.na(data$clog_upper))
  n_cover <- sum(data$clog_lower <= beta & data$clog_upper >= beta, na.rm = TRUE)
  prop <- n_cover / n_total
  # Wilson score confidence interval for binomial proportion
  z <- qnorm(1 - 0.05/2)
  denom <- 1 + z^2/n_total
  center <- (prop + z^2/(2*n_total)) / denom
  margin <- z * sqrt((prop*(1-prop) + z^2/(4*n_total))/n_total) / denom
  return(data.frame(
    prop = prop,
    lower = max(0, center - margin),
    upper = min(1, center + margin),
    n = n_total
  ))
}



# Summarize for type 1 error (MH)
type1_mh_df <- res %>%
  dplyr::group_by(strat_app, algo, n.cova, event.prob, ttm.prob) %>%
  dplyr::summarise(calc_type1_error(mh_pval), .groups = 'drop')

# Summarize for type 1 error (clogit)
type1_clog_df <- res %>%
  dplyr::group_by(strat_app, algo, n.cova, event.prob, ttm.prob) %>%
  dplyr::summarise(calc_type1_error(clog_pval), .groups = 'drop')

# Summarize for coverage
coverage_df <- res %>%
  dplyr::group_by(strat_app, algo, n.cova, event.prob, ttm.prob) %>%
  dplyr::summarise(calc_coverage(cur_data(),beta=beta), .groups = 'drop')

# Prepare for coefficient distribution
coef_all <- res

# Custom facet labels
algo_labeller <- function(x) {
  c("A1" = "Algorithm 1", "A2" = "Algorithm 2")[x]
}

ttm_labeller <- function(x) {
  c("0.1" = "Treatment probability 0.1",
    "0.5" = "Treatment probability 0.5")[x]
}

# Plot 1: Type 1 Error (MH p-value) with fixed 95% reference band and 5% line
plot_type1_mh <- ggplot(type1_mh_df, aes(x = event.prob, y = prop, fill = strat_app)) +
  #annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05 - 1.96*sqrt(0.05*0.95/unique(type1_mh_df$n)[1]), ymax = 0.05 + 1.96*sqrt(0.05*0.95/unique(type1_mh_df$n)[1]), alpha = 0.15, fill = "grey") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, width = 0.8) +
  geom_hline(yintercept = c(0.05 - 1.96*sqrt(0.05*0.95/unique(type1_mh_df$n)[1]), 0.05 + 1.96*sqrt(0.05*0.95/unique(type1_mh_df$n)[1])), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 0.5) +
  facet_grid(n.cova ~ algo + ttm.prob, labeller = labeller(algo = algo_labeller, ttm.prob = ttm_labeller, .default = label_value), switch = "y") +
  scale_fill_manual(values = c("#1B5E20", "#66BB6A", "#B71C1C", "#EF5350", "#0D47A1", "#42A5F5")) +
  labs(title = "Type 1 Error (MH p-value)", x = "Event Probability", y = "Proportion p < 0.05", fill = "Approach") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), strip.background = element_rect(fill = "grey95"))




# Plot 2: Type 1 Error (Clogit p-value) with fixed 95% reference band and 5% line
plot_type1_clog <- ggplot(type1_clog_df, aes(x = event.prob, y = prop, fill = strat_app)) +
  #annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05 - 1.96*sqrt(0.05*0.95/unique(type1_clog_df$n)[1]), ymax = 0.05 + 1.96*sqrt(0.05*0.95/unique(type1_clog_df$n)[1]), alpha = 0.15, fill = "grey") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, width = 0.8) +
  geom_hline(yintercept = c(0.05 - 1.96*sqrt(0.05*0.95/unique(type1_mh_df$n)[1]), 0.05 + 1.96*sqrt(0.05*0.95/unique(type1_mh_df$n)[1])), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 0.5) +
  facet_grid(n.cova ~ algo + ttm.prob, labeller = labeller(algo = algo_labeller, ttm.prob = ttm_labeller, .default = label_value), switch = "y") +
  scale_fill_manual(values = c("#1B5E20", "#66BB6A", "#B71C1C", "#EF5350", "#0D47A1", "#42A5F5")) +
  labs(title = "Type 1 Error (Clogit p-value)", x = "Event Probability", y = "Proportion p < 0.05", fill = "Approach") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), strip.background = element_rect(fill = "grey95"))



# Plot 3: Distribution of clog_coef (boxplot)
plot_coef_dist <- ggplot(coef_all, aes(x = event.prob, y = clog_coef, fill = strat_app)) +
  geom_boxplot(position = position_dodge(width = 0.9), outlier.size = 0.8, width = 0.7) +
  geom_hline(yintercept = beta, linetype = "dashed", color = "red", size = 0.5) +
  facet_grid(n.cova ~ algo + ttm.prob, labeller = labeller(algo = algo_labeller, ttm.prob = ttm_labeller, .default = label_value), switch = "y") +
  scale_fill_manual(values = c("#1B5E20", "#66BB6A", "#B71C1C", "#EF5350", "#0D47A1", "#42A5F5")) +
  labs(title = "Distribution of Clogit Coefficients", x = "Event Probability", y = "Coefficient Value", fill = "Approach") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), strip.background = element_rect(fill = "grey95")) +
  coord_cartesian(ylim = c(-3, 3))


# Plot 4: Coverage (bar)
plot_coverage <- ggplot(coverage_df, aes(x = event.prob, y = prop, fill = strat_app)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, width = 0.8) +
  #geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.9), width = 0.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 0.5) +
  facet_grid(n.cova ~ algo + ttm.prob, labeller = labeller(algo = algo_labeller, ttm.prob = ttm_labeller, .default = label_value), switch = "y") +
  scale_fill_manual(values = c("#1B5E20", "#66BB6A", "#B71C1C", "#EF5350", "#0D47A1", "#42A5F5")) +
  labs(title = "Coverage (95% CI includes true parameter)", x = "Event Probability", y = "Coverage Proportion", fill = "Approach") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), strip.background = element_rect(fill = "grey95")) +
  coord_cartesian(ylim = c(0.8, 1))
  #scale_y_continuous(limits = c(0.5, 1.1), breaks = seq(0.5, 1.1, 0.05))


# ============================================================================
# Save plots
# ============================================================================
# Save as PDF with multiple pages
pdf("des/plot_results.pdf", width = 12, height = 8)

print(plot_type1_mh)
print(plot_type1_clog)
print(plot_coef_dist)
print(plot_coverage)

dev.off()

# Also save as individual PNG files
ggsave("des/plot_type1_mh.png", plot_type1_mh, width = 12, height = 8, dpi = 300)
ggsave("des/plot_type1_clog.png", plot_type1_clog, width = 12, height = 8, dpi = 300)
ggsave("des/plot_coef_dist.png", plot_coef_dist, width = 14, height = 10, dpi = 300)
ggsave("des/plot_coverage.png", plot_coverage, width = 12, height = 8, dpi = 300)

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("\nType 1 Error (MH p-value):\n")
print(type1_mh_df)

cat("\nType 1 Error (Clogit p-value):\n")
print(type1_clog_df)

cat("\nCoverage:\n")
print(coverage_df)

cat("\nPlots saved to:\n")
cat("- plot_results.pdf (all plots)\n")
cat("- plot_type1_mh.png\n")
cat("- plot_type1_clog.png\n")
cat("- plot_coef_dist.png\n")
cat("- plot_coverage.png\n")
