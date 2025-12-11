

#devtools::install_local("/Users/skodikara/Documents/GitHub/tens-priv")
library(sva)
library(mixOmics)
library(tensorOmics)
library(dplyr)
library(tidyr)
library(plotly)
library(gridExtra)
library(patchwork)
library(tibble)
library(ggrepel)
library(PLSDAbatch)
library(Harman)
library(scales)

set.seed(9991)
load("data_filtered/fmtStudy.RData") 

custom_colors <- c("FMT" = "orange", "Placebo" = "blue")
weeks <- sort(unique(metadata_common$Week_mod))

# -------------------------------------------------------------------
# Helper function: merge metadata & feature matrix
# -------------------------------------------------------------------
merge_meta <- function(meta, dat) {
  merge(meta, dat, by = "row.names", all = TRUE)
}

merged_OTU_df     <- merge_meta(metadata_common, mOTU_clr_final)
merged_KO_df      <- merge_meta(metadata_common, ko_log_final)
merged_Protein_df <- merge_meta(metadata_common, prot_scale_final)

# -------------------------------------------------------------------
# Extract n, p, t
# -------------------------------------------------------------------
n <- length(unique(metadata_common$SampleID))
p_motu <- ncol(mOTU_clr_final)
p_ko   <- ncol(ko_log_final)
p_prot <- ncol(prot_scale_final)
t      <- length(weeks)

# -------------------------------------------------------------------
# Build batch variable
# (If possible, replace this with actual batch info)
# -------------------------------------------------------------------
subject_info <- metadata_common %>%
  select(SampleID, Group) %>%
  distinct() %>%
  mutate(batch_group = c("Batch2","Batch1","Batch2","Batch1","Batch2","Batch1","Batch1","Batch1","Batch2",
                         "Batch1","Batch1","Batch2","Batch1","Batch2","Batch2","Batch1","Batch1","Batch2",
                         "Batch2","Batch1","Batch2"))

group <- subject_info$Group
batch <- subject_info$batch_group

# -------------------------------------------------------------------
# Allocate arrays
# -------------------------------------------------------------------
mic_array        <- array(NA, dim = c(n, p_motu, t))
mic_nobatch_array<- array(NA, dim = c(n, p_motu, t))
colnames(mic_array)<-colnames(mic_nobatch_array)<-colnames(mOTU_clr_final)

ko_array         <- array(NA, dim = c(n, p_ko, t))
ko_nobatch_array <- array(NA, dim = c(n, p_ko, t))
colnames(ko_array)<-colnames(ko_nobatch_array)<-colnames(ko_log_final)

prot_array       <- array(NA, dim = c(n, p_prot, t))
colnames(prot_array)<-colnames(prot_scale_final)

# -------------------------------------------------------------------
# Helper: apply PLSDA batch correction to a matrix X
# -------------------------------------------------------------------
run_batch_correction <- function(X, group, batch) {
  PLSDA_batch(X, Y.trt = group, Y.bat = batch,
              ncomp.trt = 1, ncomp.bat = 1)$X.nobatch
}

# -------------------------------------------------------------------
# Fill arrays for each week
# -------------------------------------------------------------------
for (i in seq_along(weeks)) {
  
  wk <- weeks[i]
  
  # Sample order for this week
  week_ids <- paste0(subject_info$SampleID, "_W", wk)
  idx      <- match(week_ids, rownames(metadata_common))
  
  # Raw arrays
  mic_array[,,i]  <- as.matrix(mOTU_clr_final[idx, ])
  ko_array[,,i]   <- as.matrix(ko_log_final[idx, ])
  prot_array[,,i] <- as.matrix(prot_scale_final[idx, ])
  
  # Batch-corrected arrays
  mic_nobatch_array[,,i] <- run_batch_correction(mOTU_clr_final[idx, ], group, batch)
  ko_nobatch_array[,,i]  <- run_batch_correction(ko_log_final[idx, ], group, batch)
}

# -------------------------------------------------------------------
# Helper: plot PCA for any matrix
# -------------------------------------------------------------------
plot_pca <- function(mat, group, batch, week_label, title_prefix) {
  
  p <- pca(mat, ncomp = 2)
  df <- data.frame(
    Group = group,
    Batch = batch,
    PC1   = p$variates$X[,1],
    PC2   = p$variates$X[,2],
    scale =TRUE
  )
  
  ggplot(df, aes(PC1, PC2, color = Group, shape = Batch)) +
    geom_point(size = 3) +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = c(1, 17)) +
    labs(title = paste0(title_prefix, " — Week ", week_label),
         x = "PC1", y = "PC2") +
    theme_bw() +
    theme(legend.position = "none")
}

# -------------------------------------------------------------------
# Generate PCA plots for all omics × weeks × batch/non-batch
# -------------------------------------------------------------------
plots <- list()

for (i in seq_along(weeks)) {
  wk <- weeks[i]
  
  plots[[paste0("mOTU_batch_w", wk)]]     <- plot_pca(mic_array[,,i], group, batch, wk, "mOTU (batch)")
  plots[[paste0("mOTU_nobatch_w", wk)]]   <- plot_pca(mic_nobatch_array[,,i], group, batch, wk, "mOTU (no batch)")
  
  plots[[paste0("KO_batch_w", wk)]]       <- plot_pca(ko_array[,,i], group, batch, wk, "KO (batch)")
  plots[[paste0("KO_nobatch_w", wk)]]     <- plot_pca(ko_nobatch_array[,,i], group, batch, wk, "KO (no batch)")
  
  plots[[paste0("Protein_w", wk)]]        <- plot_pca(prot_array[,,i], group, batch, wk, "Protein")
}

# -------------------------------------------------------------------
#Week 0
plots$mOTU_batch_w0
plots$mOTU_nobatch_w0

plots$KO_batch_w0
plots$KO_nobatch_w0

plots$Protein_w0

#Week 8
plots$mOTU_batch_w8
plots$mOTU_nobatch_w8

plots$KO_batch_w8
plots$KO_nobatch_w8

plots$Protein_w8

#Week 24
plots$mOTU_batch_w24
plots$mOTU_nobatch_w24

plots$KO_batch_w24
plots$KO_nobatch_w24

plots$Protein_w24

# -------------------------------------------------------------------
# Tensor approach
# -------------------------------------------------------------------

plot_tpca <- function(tensor, group, batch_group, title = "") {
  
  # run tPCA
  tpca_res <- tpca(tensor, ncomp = 2)
  
  # create dataframe for plotting
  df <- data.frame(
    Group = group,
    Batch = batch_group,
    PC1 = tpca_res$variates[, 1],
    PC2 = tpca_res$variates[, 2]
  )
  
  # ggplot
  ggplot(df, aes(PC1, PC2, color = Group, shape = Batch)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("orange", "blue")) +
    scale_shape_manual(values = c(1, 17)) +
    labs(
      title = title,
      x = "Comp1", y = "Comp2"
    ) +
    theme_bw() +
    theme(legend.position = "none")
}

# microbiome
plot_tpca(mic_array, group, subject_info$batch_group, "tPCA - microbiome")

# microbiome (no batch)
plot_tpca(mic_nobatch_array, group, subject_info$batch_group, "tPCA - microbiome (no batch)")

# KO
plot_tpca(ko_array, group, subject_info$batch_group, "tPCA - KO")

# KO (no batch)
plot_tpca(ko_nobatch_array, group, subject_info$batch_group, "tPCA - KO (no batch)")

# Protein
plot_tpca(prot_array, group, subject_info$batch_group, "tPCA - Protein")




# -------------------------------------------------------------------
# Save final objects
# -------------------------------------------------------------------
save(
  taxonomy_data, ko_info,
  mic_nobatch_array, ko_nobatch_array,
  prot_array, metadata_common, subject_info,
  file = "data_filtered/fmtStudy_noBatch.RData"
)

