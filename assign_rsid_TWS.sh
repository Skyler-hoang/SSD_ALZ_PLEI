#!/bin/bash
#SBATCH --job-name=assign_rsid_TWS
#SBATCH --cpus-per-task=3
#SBATCH --mem=128G
#SBATCH --time=84:00:00
#SBATCH --output=assign_rsid_TWS_%j.out
#SBATCH --error=assign_rsid_TWS_%j.err

# Load R module
module load R/4.4.2-gfbf-2024a

# Run R script
R --vanilla << 'EOF'

# Load required libraries
library(data.table)
library(Rsamtools)  # For tabix functionality

# Set directory - now using only one directory for both input and output
work_dir <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/06_GWAS_Assigning_rsid"
dbsnp_file <- file.path(work_dir, "00-common_all.vcf.gz")

# Process only TWS phenotype
phenotypes <- c("TWS")

print("=== FORMATTING GWAS RESULTS FOR GENOMICSEM ===")
print(paste("Working directory:", work_dir))
print(paste("dbSNP reference file:", dbsnp_file))

# Check if dbSNP file exists
if (!file.exists(dbsnp_file)) {
  stop("dbSNP file not found. Please make sure it's in the correct location.")
}
if (!file.exists(paste0(dbsnp_file, ".tbi"))) {
  stop("dbSNP index file (.tbi) not found. Please make sure it's in the correct location.")
}

# List the .Robj files in the directory to confirm they're accessible
robj_files <- list.files(work_dir, pattern = "\\.Robj$", full.names = TRUE)
print(paste("Found", length(robj_files), "Robj files in the working directory:"))
for (file in robj_files) {
  print(paste("  ", basename(file)))
}

# Function to extract rsID from VCF based on chromosome and position
get_rsids <- function(chrom, pos, vcf_file) {
  # Ensure consistent data types
  chrom <- as.character(chrom)
  pos <- as.integer(pos)
  
  # Remove any 'chr' prefix for consistency
  chrom <- gsub("^chr", "", chrom)
  
  # Create result vector
  rsids <- rep(NA_character_, length(chrom))
  
  # Process in smaller chunks to avoid memory issues
  chunk_size <- 500
  n_chunks <- ceiling(length(chrom) / chunk_size)
  
  for (i in 1:n_chunks) {
    start_idx <- (i-1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(chrom))
    
    chunk_chrom <- chrom[start_idx:end_idx]
    chunk_pos <- pos[start_idx:end_idx]
    
    # Process each SNP in the chunk
    for (j in 1:length(chunk_chrom)) {
      idx <- start_idx + j - 1
      
      # Skip if chromosome or position is invalid
      if (is.na(chunk_chrom[j]) || is.na(chunk_pos[j]) || chunk_chrom[j] == "" || chunk_pos[j] <= 0) {
        next
      }
      
      # Try to fetch from tabix - ensure chromosome format matches VCF
      tryCatch({
        # Try without 'chr' prefix first
        chr_query <- chunk_chrom[j]
        param <- GRanges(chr_query, IRanges(chunk_pos[j], chunk_pos[j]))
        vcf_records <- scanTabix(vcf_file, param=param)
        
        # If no match, try with 'chr' prefix
        if (length(vcf_records[[1]]) == 0) {
          chr_query <- paste0("chr", chunk_chrom[j])
          param <- GRanges(chr_query, IRanges(chunk_pos[j], chunk_pos[j]))
          vcf_records <- scanTabix(vcf_file, param=param)
        }
        
        # If we found a match
        if (length(vcf_records[[1]]) > 0) {
          # Parse the VCF line to extract rsID (3rd column)
          vcf_line <- unlist(strsplit(vcf_records[[1]][1], "\t"))
          if (length(vcf_line) >= 3 && vcf_line[3] != ".") {
            rsids[idx] <- vcf_line[3]
          }
        }
      }, error = function(e) {
        # Just continue with NA if there's an error
      })
    }
    
    # Print progress
    if (i %% 50 == 0) {
      print(paste("    Processed chunk", i, "of", n_chunks))
    }
  }
  
  return(rsids)
}

# Alternative function using fread - with fixed data type handling
get_rsids_fread <- function(chrom, pos, vcf_file) {
  print("  Using fread method to read VCF - this may take a while...")
  
  # Ensure consistent data types
  chrom <- as.character(chrom)
  pos <- as.integer(pos)
  
  # Remove any 'chr' prefix for consistency
  chrom <- gsub("^chr", "", chrom)
  
  tryCatch({
    # Read VCF file with only the columns we need
    print("  Reading VCF file...")
    vcf_data <- fread(cmd = paste("zcat", vcf_file, "| grep -v '^##'"), 
                     select = c(1,2,3), 
                     col.names = c("CHROM", "POS", "rsID"),
                     colClasses = c("character", "integer", "character"))
    
    # Ensure chromosome format consistency
    vcf_data[, CHROM := gsub("^chr", "", CHROM)]
    
    # Create query data.table with consistent types
    query <- data.table(CHROM = chrom, POS = pos)
    
    # Set keys for efficient joining
    setkey(vcf_data, CHROM, POS)
    setkey(query, CHROM, POS)
    
    print("  Matching positions...")
    # Merge to get rsIDs
    result <- vcf_data[query, nomatch=NA]
    
    # Replace "." with NA for missing rsIDs
    result[rsID == ".", rsID := NA]
    
    return(result$rsID)
    
  }, error = function(e) {
    print(paste("  Error in fread method:", e$message))
    return(rep(NA_character_, length(chrom)))
  })
}

# Function to check if chromosome format is numeric or "chr" prefix
detect_chr_format <- function(vcf_file) {
  tryCatch({
    # Read the first few data lines to check chromosome format
    vcf_sample <- fread(cmd = paste("zcat", vcf_file, "| grep -v '^##' | head -n 5"), 
                        header = TRUE, sep = "\t", nrows = 5)
    
    if (ncol(vcf_sample) > 0 && nrow(vcf_sample) > 0) {
      first_chrom <- as.character(vcf_sample[[1]][1])
      if (grepl("^chr", first_chrom)) {
        return("chr")
      } else {
        return("numeric")
      }
    }
  }, error = function(e) {
    print(paste("Warning: Could not detect chromosome format:", e$message))
  })
  
  # Default to numeric if we can't determine
  return("numeric")
}

# Determine chromosome format in VCF
chr_format <- detect_chr_format(dbsnp_file)
print(paste("Detected chromosome format in VCF:", chr_format))

# Process each phenotype
for (pheno in phenotypes) {
  # Construct file path to the copied files
  file_path <- file.path(work_dir, paste0(pheno, "_withParents_allCHR_MAF05.Robj"))
  
  # Skip if file doesn't exist
  if (!file.exists(file_path)) {
    print(paste("Skipping", pheno, "- file not found:", file_path))
    next
  }
  
  print(paste("\nProcessing phenotype:", pheno))
  print(paste("  Loading file:", file_path))
  
  # Load the Robj file
  tryCatch({
    load(file_path)
    
    # Get all objects in the environment
    all_objs <- ls()
    
    # Remove known objects from this script
    script_objs <- c("work_dir", "dbsnp_file", "phenotypes", "pheno", 
                    "file_path", "chr_format", "get_rsids", "get_rsids_fread", 
                    "detect_chr_format", "all_objs", "script_objs", "robj_files")
    
    gwas_objs <- setdiff(all_objs, script_objs)
    
    if (length(gwas_objs) == 0) {
      print(paste("  Error: No objects found in", file_path))
      next
    }
    
    # Try to identify the GWAS results object
    gwas_obj_name <- NULL
    if ("data05" %in% gwas_objs) {
      gwas_obj_name <- "data05"
    } else if ("my_r" %in% gwas_objs) {
      gwas_obj_name <- "my_r"
    } else {
      # Otherwise take the first object
      gwas_obj_name <- gwas_objs[1]
    }
    
    print(paste("  Found data object:", gwas_obj_name))
    
    # Get the GWAS data
    gwas_data <- get(gwas_obj_name)
    
    # Convert to data.table if not already
    if (!is.data.table(gwas_data)) {
      gwas_data <- as.data.table(gwas_data)
    }
    
    # Print original column names to understand the data
    print("  Original column names:")
    print(colnames(gwas_data))
    
    # Verify we have the required columns
    if (!all(c("CHROM", "POS") %in% colnames(gwas_data))) {
      # Try to accommodate variations in column naming
      if ("CHR" %in% colnames(gwas_data) && "BP" %in% colnames(gwas_data)) {
        setnames(gwas_data, c("CHR", "BP"), c("CHROM", "POS"))
        print("  Renamed CHR/BP to CHROM/POS")
      } else {
        print("  Column names in GWAS data:")
        print(colnames(gwas_data))
        stop(paste("  Required columns not found in", pheno, "file"))
      }
    }
    
    # Ensure chromosome and position are correct data types
    gwas_data[, CHROM := as.character(CHROM)]
    gwas_data[, POS := as.integer(POS)]
    
    # Identify columns
    ref_col <- "REF"
    alt_col <- "ALT"
    effect_col <- "Beta"
    pval_col <- "Pvalue"
    
    # Verify these columns exist
    if (!all(c(ref_col, alt_col, effect_col, pval_col) %in% colnames(gwas_data))) {
      print("  Warning: Some expected columns not found. Checking alternatives...")
      
      # Look for alternative column names
      if (!"REF" %in% colnames(gwas_data)) {
        if ("A1" %in% colnames(gwas_data)) {
          ref_col <- "A1"
        } else if ("A2" %in% colnames(gwas_data)) {
          ref_col <- "A2"
        }
      }
      
      if (!"ALT" %in% colnames(gwas_data)) {
        if ("A2" %in% colnames(gwas_data)) {
          alt_col <- "A2"
        } else if ("A1" %in% colnames(gwas_data)) {
          alt_col <- "A1"
        }
      }
      
      if (!"Beta" %in% colnames(gwas_data)) {
        if ("beta" %in% colnames(gwas_data)) {
          effect_col <- "beta"
        } else if ("BETA" %in% colnames(gwas_data)) {
          effect_col <- "BETA"
        }
      }
      
      if (!"Pvalue" %in% colnames(gwas_data)) {
        if ("P" %in% colnames(gwas_data)) {
          pval_col <- "P"
        } else if ("pvalue" %in% colnames(gwas_data)) {
          pval_col <- "pvalue"
        }
      }
    }
    
    print(paste("  GWAS data has", nrow(gwas_data), "SNPs"))
    print(paste("  REF column:", ref_col))
    print(paste("  ALT column:", alt_col))
    print(paste("  Effect column:", effect_col))
    print(paste("  P-value column:", pval_col))
    
    # Get rsIDs for all SNPs
    print("  Assigning rsIDs to SNPs...")
    
    # Try tabix method first
    rsids <- get_rsids(gwas_data$CHROM, gwas_data$POS, dbsnp_file)
    
    # Count successful assignments
    n_with_rsid <- sum(!is.na(rsids))
    success_rate <- n_with_rsid / nrow(gwas_data)
    
    print(paste("  Tabix method assigned rsIDs to", n_with_rsid, "out of", nrow(gwas_data), 
               "SNPs (", round(success_rate * 100, 1), "%)"))
    
    # If tabix method didn't find many rsIDs, try fread method
    if (success_rate < 0.10) {
      print("  Low success rate with tabix method, trying fread method...")
      rsids_fread <- get_rsids_fread(gwas_data$CHROM, gwas_data$POS, dbsnp_file)
      
      # Use fread results if they're better
      n_with_rsid_fread <- sum(!is.na(rsids_fread))
      if (n_with_rsid_fread > n_with_rsid) {
        rsids <- rsids_fread
        n_with_rsid <- n_with_rsid_fread
        print(paste("  Fread method assigned rsIDs to", n_with_rsid, "out of", nrow(gwas_data), 
                   "SNPs (", round(n_with_rsid/nrow(gwas_data)*100, 1), "%)"))
      }
    }
    
    # Add rsIDs to GWAS data
    gwas_data[, SNP := rsids]
    
    # Create a unique identifier for SNPs without rsID
    gwas_data[is.na(SNP), SNP := paste0("chr", CHROM, ":", POS, "_", get(ref_col), "_", get(alt_col))]
    
    # Create new data.table with GenomicSEM format
    genomicsem_data <- data.table(
      SNP = gwas_data$SNP,
      A1 = gwas_data[[alt_col]],  # Effect allele
      A2 = gwas_data[[ref_col]],  # Non-effect allele
      BETA = gwas_data[[effect_col]],
      P = gwas_data[[pval_col]]
    )
    
    # Add additional columns if available
    if ("N_INFORMATIVE" %in% colnames(gwas_data)) {
      genomicsem_data[, N := gwas_data$N_INFORMATIVE]
    } else if ("N" %in% colnames(gwas_data)) {
      genomicsem_data[, N := gwas_data$N]
    }
    
    if ("AF" %in% colnames(gwas_data)) {
      genomicsem_data[, AF := gwas_data$AF]
    }
    
    if ("BetaVar" %in% colnames(gwas_data)) {
      genomicsem_data[, SE := sqrt(gwas_data$BetaVar)]
    } else if ("SE" %in% colnames(gwas_data)) {
      genomicsem_data[, SE := gwas_data$SE]
    }
    
    # Remove any rows with missing essential data
    initial_rows <- nrow(genomicsem_data)
    genomicsem_data <- genomicsem_data[!is.na(SNP) & !is.na(BETA) & !is.na(P)]
    final_rows <- nrow(genomicsem_data)
    
    if (initial_rows != final_rows) {
      print(paste("  Removed", initial_rows - final_rows, "rows with missing essential data"))
    }
    
    # Save to output file
    output_file <- file.path(work_dir, paste0(pheno, "_GenomicSEM_format.txt"))
    fwrite(genomicsem_data, output_file, sep="\t", na="NA")
    
    print(paste("  Saved", nrow(genomicsem_data), "SNPs to:", output_file))
    
    # Also save the full data with rsIDs for reference
    full_output_file <- file.path(work_dir, paste0(pheno, "_withRSIDs_full.txt"))
    fwrite(gwas_data, full_output_file, sep="\t", na="NA")
    
    print(paste("  Saved full data with rsIDs to:", full_output_file))
    
    # Clean up to avoid conflicts with the next file
    rm(list = gwas_objs)
    
  }, error = function(e) {
    print(paste("  Error processing", pheno, ":", e$message))
    print("  Skipping this phenotype")
  })
}

# Create a summary of the process
print("\n=== SUMMARY ===")
output_files <- list.files(work_dir, pattern = "_GenomicSEM_format\\.txt$", full.names = TRUE)
print(paste("Total GenomicSEM format files created:", length(output_files)))

# Check file sizes to make sure they're not empty
if (length(output_files) > 0) {
  file_sizes <- file.info(output_files)$size
  print("Output file sizes (bytes):")
  for (i in seq_along(output_files)) {
    print(paste("  ", basename(output_files[i]), ":", file_sizes[i]))
  }
}

print("\n=== GENOMICSEM FORMAT CONVERSION COMPLETE ===")

EOF
