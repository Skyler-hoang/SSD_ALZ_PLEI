#!/bin/bash
#SBATCH --job-name=alzheimer_pgs_analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=alzheimer_pgs_%j.out
#SBATCH --error=alzheimer_pgs_%j.err

# Load R module
module load R/4.4.2-gfbf-2024a

# Run R script
R --vanilla << 'EOF'

# Load required libraries
library(data.table)

# Set up graphics for headless environment
options(bitmapType = "cairo")

# Set directories
results_dir <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/03_B_Omni2013_Michigan_1KG30X_PGS"
clinical_file <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/A_Phonology_GWAS_435_samples.csv"
scores_file <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/03_A_Omni2013_Michigan_1KG30X_PGS/scores.txt"
output_results_dir <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/06_Alzheimer_results"
output_plots_dir <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/06_Alzheimer_plots"

# Create output directories
dir.create(output_results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_plots_dir, showWarnings = FALSE, recursive = TRUE)

print("=== ALZHEIMER'S PGS ANALYSIS ===")

# Define Alzheimer's disease PGS IDs
alzheimer_pgs <- c(
  "PGS002249", "PGS002753", "PGS003953", "PGS003957", "PGS003958", 
  "PGS003992", "PGS004008", "PGS004034", "PGS004062", "PGS004092", 
  "PGS004228", "PGS004229", "PGS004146"
)

print(paste("Target Alzheimer's PGS count:", length(alzheimer_pgs)))

# Define clinical measures
clinical_measures <- c(
  "NTSM_MSW_PPC_Percentage_age_adjusted_t",
  "NTSM_NSW_PPC_Percentage_age_adjusted_t", 
  "NTSM_NSW_Percentage_age_adjusted",
  "NTSM_MSW_Percentage_age_adjusted",
  "CELF_TOLD_Receptive",
  "CELF_TOLD_Expressive",
  "CTPP_PA_Elision_Z_t",
  "CTPP_RAN_Colors_Z_t",
  "WRMT_Ident_Standard_t",
  "WRMT_Attack_Standard_t",
  "TWS_Total_Standard_t",
  "EOWPVT_Standard",
  "PPVT_Standard",
  "FL_OSMCP_Function_comb_z_t",
  "FL_OSMCP_Function_comb_z",
  "Age_FL_OSMCP_Function_comb_z",
  "FL_OSMCP_Function_comb_z_qn_t"
)

# Load data for stratification analysis
print("Loading data for stratification analysis...")
clinical_data <- fread(clinical_file)
scores_data <- fread(scores_file)

# Create combined IDs for merging
clinical_data[, combined_id := paste(FID, IID, sep = "_")]
scores_data[, combined_id := sample]

print(paste("Clinical data samples:", nrow(clinical_data)))
print(paste("Scores data samples:", nrow(scores_data)))

# Function to read linear regression results from TSV files
read_linear_regression_results <- function(clinical_measure, alzheimer_pgs) {
  tsv_file <- file.path(results_dir, paste0(clinical_measure, "_PGS_results.tsv"))
  
  if (!file.exists(tsv_file)) {
    print(paste("Warning: TSV file not found for", clinical_measure))
    return(NULL)
  }
  
  # Read TSV file
  results <- fread(tsv_file)
  
  # Filter for Alzheimer's PGS only
  alzheimer_results <- results[PGS %in% alzheimer_pgs]
  
  if (nrow(alzheimer_results) == 0) {
    print(paste("Warning: No Alzheimer's PGS found in", clinical_measure, "results"))
    return(NULL)
  }
  
  # Add clinical measure column
  alzheimer_results[, clinical_measure := clinical_measure]
  
  return(alzheimer_results[, .(PGS, clinical_measure, p_value, beta, r_squared)])
}

# Function to perform stratified analysis for one PGS-clinical pair
perform_stratified_analysis <- function(pgs_name, clinical_measure, clinical_data, scores_data, output_plots_dir) {
  
  # Check if clinical measure exists
  if (!clinical_measure %in% names(clinical_data)) {
    print(paste("Clinical measure", clinical_measure, "not found - skipping"))
    return(NULL)
  }
  
  # Check if PGS exists
  if (!pgs_name %in% names(scores_data)) {
    print(paste("PGS", pgs_name, "not found - skipping"))
    return(NULL)
  }
  
  # Merge clinical and PGS data
  analysis_data <- merge(clinical_data[, .(combined_id, get(clinical_measure))], 
                        scores_data[, .(combined_id, get(pgs_name))], 
                        by = "combined_id")
  
  # Rename columns for clarity
  setnames(analysis_data, c("combined_id", "clinical_score", "pgs_score"))
  
  # Remove missing values
  analysis_data <- analysis_data[!is.na(clinical_score) & !is.na(pgs_score)]
  
  if (nrow(analysis_data) < 20) {
    print(paste("Insufficient data for", pgs_name, "vs", clinical_measure, "- skipping"))
    return(NULL)
  }
  
  # Stratify into bottom 50% and top 50% based on clinical measure
  n_samples <- nrow(analysis_data)
  n_per_group <- floor(n_samples / 2)
  
  # Sort by clinical measure and select top and bottom 50%
  analysis_data <- analysis_data[order(clinical_score)]
  
  bottom_50 <- analysis_data[1:n_per_group]
  top_50 <- analysis_data[(n_samples - n_per_group + 1):n_samples]
  
  # Test normality using Shapiro-Wilk test
  shapiro_bottom <- tryCatch({
    shapiro.test(bottom_50$pgs_score)
  }, error = function(e) {
    list(p.value = NA)
  })
  
  shapiro_top <- tryCatch({
    shapiro.test(top_50$pgs_score)
  }, error = function(e) {
    list(p.value = NA)
  })
  
  # Determine if data is normally distributed (p > 0.05 in Shapiro-Wilk)
  normal_bottom <- !is.na(shapiro_bottom$p.value) && shapiro_bottom$p.value > 0.05
  normal_top <- !is.na(shapiro_top$p.value) && shapiro_top$p.value > 0.05
  both_normal <- normal_bottom & normal_top
  
  # Perform statistical tests
  # Welch t-test
  t_test_result <- tryCatch({
    t.test(top_50$pgs_score, bottom_50$pgs_score, var.equal = FALSE)
  }, error = function(e) {
    list(p.value = NA)
  })
  
  # Wilcoxon test
  wilcox_result <- tryCatch({
    wilcox.test(top_50$pgs_score, bottom_50$pgs_score)
  }, error = function(e) {
    list(p.value = NA)
  })
  
  # Determine recommendation
  recommendation <- ifelse(both_normal, "Welch", "Wilcoxon")
  
  # Create density plot
  plot_filename <- paste0(pgs_name, "_vs_", clinical_measure, "_density_plot.png")
  plot_path <- file.path(output_plots_dir, plot_filename)
  
  tryCatch({
    png(plot_path, width = 3600, height = 2400, res = 300, type = "cairo")
    
    # Calculate density curves
    bottom_density <- density(bottom_50$pgs_score)
    top_density <- density(top_50$pgs_score)
    
    # Set up plot limits
    x_range <- range(c(bottom_density$x, top_density$x))
    y_range <- range(c(bottom_density$y, top_density$y))
    
    # Create plot
    plot(bottom_density, 
         col = "blue", 
         lwd = 2,
         xlim = x_range,
         ylim = y_range,
         main = paste("Alzheimer's PGS:", pgs_name),
         sub = paste("vs", gsub("_", " ", clinical_measure), "| Welch p =", 
                    ifelse(is.na(t_test_result$p.value), "NA", round(t_test_result$p.value, 4)), 
                    "| Wilcoxon p =", 
                    ifelse(is.na(wilcox_result$p.value), "NA", round(wilcox_result$p.value, 4)),
                    "| Recommendation:", recommendation),
         xlab = "PGS Score",
         ylab = "Density")
    
    # Add top 50% density
    lines(top_density, col = "red", lwd = 2)
    
    # Fill areas under curves
    polygon(bottom_density, col = rgb(0, 0, 1, 0.3), border = NA)
    polygon(top_density, col = rgb(1, 0, 0, 0.3), border = NA)
    
    # Add legend
    legend("topright", 
           legend = c(paste("Bottom 50% (n =", nrow(bottom_50), ")"),
                     paste("Top 50% (n =", nrow(top_50), ")")),
           col = c("blue", "red"),
           lwd = 2,
           cex = 0.8)
    
    # Add grid
    grid(col = "lightgray", lty = 1)
    
    dev.off()
    
  }, error = function(e) {
    dev.off()
    print(paste("Error creating plot:", e$message))
  })
  
  # Return results
  return(data.table(
    PGS_ID = pgs_name,
    Clinical_Measure = clinical_measure,
    Welch_p_value = ifelse(is.na(t_test_result$p.value), NA, t_test_result$p.value),
    Wilcoxon_p_value = ifelse(is.na(wilcox_result$p.value), NA, wilcox_result$p.value),
    Recommendation = recommendation,
    Sample_size = nrow(analysis_data)
  ))
}

# Step 1: Read linear regression results for all clinical measures
print("Step 1: Reading linear regression results...")

all_linear_results <- list()

for (clinical_measure in clinical_measures) {
  linear_result <- read_linear_regression_results(clinical_measure, alzheimer_pgs)
  
  if (!is.null(linear_result)) {
    all_linear_results[[length(all_linear_results) + 1]] <- linear_result
  }
}

if (length(all_linear_results) == 0) {
  stop("No linear regression results found for any clinical measures")
}

# Combine all linear regression results
combined_linear_results <- rbindlist(all_linear_results)

print(paste("Total linear regression results:", nrow(combined_linear_results)))
print(paste("Unique PGS found:", length(unique(combined_linear_results$PGS))))
print(paste("Unique clinical measures found:", length(unique(combined_linear_results$clinical_measure))))

# Step 2: Perform stratified analysis for all PGS-clinical combinations
print("Step 2: Performing stratified analyses...")

all_stratified_results <- list()
total_combinations <- length(alzheimer_pgs) * length(clinical_measures)
current_combination <- 0

for (pgs_name in alzheimer_pgs) {
  for (clinical_measure in clinical_measures) {
    current_combination <- current_combination + 1
    
    if (current_combination %% 20 == 0) {
      print(paste("Processing combination", current_combination, "of", total_combinations))
    }
    
    stratified_result <- perform_stratified_analysis(pgs_name, clinical_measure, clinical_data, scores_data, output_plots_dir)
    
    if (!is.null(stratified_result)) {
      all_stratified_results[[length(all_stratified_results) + 1]] <- stratified_result
    }
  }
}

if (length(all_stratified_results) == 0) {
  stop("No stratified analyses could be performed")
}

# Combine all stratified results
combined_stratified_results <- rbindlist(all_stratified_results)

print(paste("Total stratified analyses completed:", nrow(combined_stratified_results)))

# Step 3: Merge linear regression and stratified results
print("Step 3: Merging results...")

# Merge on PGS and clinical_measure
final_results <- merge(
  combined_linear_results,
  combined_stratified_results,
  by.x = c("PGS", "clinical_measure"),
  by.y = c("PGS_ID", "Clinical_Measure"),
  all.x = TRUE
)

# Reorder columns as requested
final_table <- final_results[, .(
  PGS_ID = PGS,
  Clinical_Measure = clinical_measure,
  p_value = p_value,
  beta = beta,
  r_squared = r_squared,
  Welch_p_value = Welch_p_value,
  Wilcoxon_p_value = Wilcoxon_p_value,
  Recommendation = Recommendation,
  Sample_size = Sample_size
)]

# Sort by PGS_ID and Clinical_Measure
final_table <- final_table[order(PGS_ID, Clinical_Measure)]

# Step 4: Save results
print("Step 4: Saving results...")

output_file <- file.path(output_results_dir, "alzheimer_pgs_complete_analysis.tsv")
fwrite(final_table, output_file, sep = "\t")

print("\n=== ANALYSIS COMPLETE ===")
print(paste("Final results saved to:", output_file))
print(paste("Density plots saved in:", output_plots_dir))
print(paste("Total analyses in final table:", nrow(final_table)))
print(paste("Alzheimer's PGS analyzed:", length(unique(final_table$PGS_ID))))
print(paste("Clinical measures analyzed:", length(unique(final_table$Clinical_Measure))))

# Display summary
print("\n=== SUMMARY STATISTICS ===")
print("PGS coverage:")
print(table(final_table$PGS_ID))

print("\nClinical measure coverage:")
print(table(final_table$Clinical_Measure))

print("\nFirst 10 rows of final results:")
print(head(final_table, 10))

print("\nTest recommendation distribution:")
print(table(final_table$Recommendation, useNA = "ifany"))

EOF
