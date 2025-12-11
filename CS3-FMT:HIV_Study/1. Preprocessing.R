library(mixOmics)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(readxl)

# Function to perform low count removal
low.removal <- function(
    data, # df of size n (sample) x p (features)
    percent=0.01 # cutoff chosen
) 
{
  keep.features = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.features]
  return(list(data.filter = data.filter, keep.features = keep.features))
}

# ============================================================================
# LOAD CLINICAL METADATA
# ============================================================================

metadata_all <- read_xlsx("data/Info patients REFRESH.xlsx", n_max = 30)


# ============================================================================
# PROCESS MICROBIOME DATA (mOTUs)
# ============================================================================

# ---- Load raw mOTU file ----
mOTU_raw <- read_tsv(
  "data/phyloseq_profile_all_motus_profile_counts_g1l70_TT.tsv",
  show_col_types = FALSE
)

# ---- Extract taxonomy (first 7 columns) ----
taxonomy_data <- mOTU_raw %>%
  select(1:7) %>%
  mutate(Taxa = paste0("Taxa", row_number()))   # create unique feature names

# ---- Transpose abundance table (samples x features) ----
mOTU_data <- mOTU_raw %>%
  select(-1:-7) %>%   # remove taxonomy columns
  t()

colnames(mOTU_data) <- taxonomy_data$Taxa       # assign labelled feature names

# ---- Create metadata for samples extracted from rownames ----
metadata_mOTU <- tibble(id = rownames(mOTU_data))


# Split IDs like "R24W8" → SampleID = "R24", Week = "8"
# And create a new week variable with numeric values only 
# (e.g., "24X" and "24Y" both become "24")
# And filter 24Y 
split_ids <- strsplit(metadata_mOTU$id, "W", fixed = TRUE)

metadata_mOTU <- metadata_mOTU %>%
  mutate(
    SampleID = sapply(split_ids, `[`, 1),
    Week = sapply(split_ids, function(x) ifelse(length(x) >= 2, x[2], NA)),
    Week_clean = case_when(
      Week %in% c("24X", "24Y") ~ "24",
      TRUE ~ Week
    ),
    Week_mod = as.numeric(Week_clean)
  )%>%select( -Week_clean)

#---- Checking metadata_OTU with actual metadata ----
sum(!(unique(metadata_mOTU$SampleID) %in% unique(metadata_all$id)))

#There are 3 SampleIDs in metadata_mOTU that are not in metadata_all: 
#"0505_0101", "543_46G1", "585_57G1"

#---- Filter out these 3 SampleIDs from metadata_mOTU ----
metadata_mOTU <- metadata_mOTU %>%
  filter(!(SampleID %in% c("0505_0101", "543_46G1", "585_57G1")))

# ---- Remove duplicated second visit 24Y ----
metadata_mOTU <- filter(metadata_mOTU, !(Week == "24Y"))

# ---- Check number of time points per individual ----
table(metadata_mOTU$SampleID,
      metadata_mOTU$Week_mod)%>%rowSums()

metadata_mOTU %>%
  count(SampleID) %>%
  filter(n < 3)

# ---- List of samples with incomplete time series ----
incomplete_samples <- c("R5", "R27", "R26", "R25", "R22", "R13")

metadata_mOTU_filtered <- metadata_mOTU %>%
  filter(!SampleID %in% incomplete_samples)

# ---- Identify samples missing required weeks 0, 8, 24 ----
samples_without_weeks <- metadata_mOTU_filtered %>%
  group_by(SampleID) %>%
  filter(!all(c(0, 8, 24) %in% Week_mod)) %>%
  pull(SampleID) %>%
  unique()

# ---- Fix inconsistent week codes for specific participants ----
metadata_mOTU_filtered %>%
  filter(SampleID %in% samples_without_weeks) %>%
  { table(.$SampleID, .$Week_mod) }

# R10, R16, R19, R21, R3, R6, R8, R9 => Missing 8 has 7 instead, thus recode 7 → 8
# "R24"  => Missing 0 has 1 instead, thus recode 1 → 0

metadata_mOTU_filtered <- metadata_mOTU_filtered %>%
  mutate(
    Week_mod = case_when(
      SampleID == "R24" & Week_mod == 1 ~ 0,        # manual correction
      SampleID %in% setdiff(samples_without_weeks, "R24") & Week_mod == 7 ~ 8,
      TRUE ~ Week_mod
    )
  ) %>%
  filter(Week_mod %in% c(0, 8, 24))

# ---- Merge with clinical metadata ----
metadata_mOTU_filtered <- metadata_mOTU_filtered %>%
  left_join(metadata_all, by = c("SampleID" = "id")) %>%
  arrange(SampleID, Week_mod) %>%
  # Generate consistent unique ID "SampleID_Week"
  mutate(newid = paste0(SampleID, "_W", Week_mod), .before = 1) %>%
  column_to_rownames("newid")

# ---- Subset mOTU matrix to matched metadata order ----
mOTU_data_filtered <- mOTU_data[
  match(metadata_mOTU_filtered$id, rownames(mOTU_data)),
  , drop = FALSE
]

# Ensure perfect alignment
sum(rownames(mOTU_data_filtered) != metadata_mOTU_filtered$id) == 0

# Update rownames to newid (This is needed as we approximated week for some samples above)
rownames(mOTU_data_filtered) <- rownames(metadata_mOTU_filtered)

# ---- Low count filtering ----
mOTU_filtered <- low.removal(mOTU_data_filtered, percent = 0.01)$data.filter

# ---- CLR transformation (compositional normalisation) ----
mOTU_filtered.clr <- logratio.transfo(
  as.matrix(mOTU_filtered),
  logratio = "CLR",
  offset = 1               # pseudocount to avoid log(0)
)


# ============================================================================
# PROCESS KEGG KO DATA
# ============================================================================


ko_raw <- read_tsv("data/ko_relabundance.tsv", show_col_types = FALSE)

ko_info <- select(ko_raw, 1:3)    # KO ID + functional annotation

ko_data <- ko_raw %>%
  select(-1:-3) %>%
  t()

colnames(ko_data) <- ko_info$KEGG_ko

# ---- Match sample IDs to mOTU metadata ----
ko_data_filtered <- ko_data[
  match(metadata_mOTU_filtered$id, rownames(ko_data)),
  , drop = FALSE
]

sum(rownames(ko_data_filtered) != metadata_mOTU_filtered$id) == 0

# change rownames to newid
rownames(ko_data_filtered) <- rownames(metadata_mOTU_filtered)


# ---- Low count filtering ----
ko_filtered <- low.removal(ko_data_filtered, percent = 0.01)$data.filter

# ---- Log transform ----
ko_filtered.log <- log(as.matrix(ko_filtered) + 1)


# ============================================================================
# PROCESS PROTEOMIC DATA
# ============================================================================

prot_raw <- read_csv("data/protdata_processed.csv", show_col_types = FALSE)

# ---- Harmonise sample ID formats ----
protdata <- prot_raw %>%
  rename(id = SampleID) %>%
  mutate(
    id       = str_replace(id, " w", "W"),  # "R24 w8" → "R24W8"
    SampleID = metadata_id
  )

# ---- Keep unique metadata rows only ----
metadata_prot <- protdata %>%
  distinct(SampleID, patientid, week, week_cat, id)

table(metadata_prot$SampleID,
      metadata_prot$week)

# ---- Remove samples flagged as incomplete ----
# R18 has no observations after 8
metadata_prot_filtered <- metadata_prot %>%
  filter(!SampleID %in% c("R18")) %>%
  mutate(Week_mod = week_cat) %>%   # week_cat has already done some approximation 
  # Merge with clinical data
  left_join(metadata_all, by = c("SampleID" = "id")) %>%
  
  arrange(SampleID, Week_mod) %>%
  mutate(newid = paste0(SampleID, "_W", Week_mod), .before = 1) %>%
  column_to_rownames("newid")

# ---- Convert long-format NPX into matrix form ----
prot_mat <- xtabs(NPX ~ id + Assay, data = protdata) %>%
  as.matrix()

# ---- Match NPX matrix to metadata ----
protdata_filtered <- prot_mat[
  match(metadata_prot_filtered$id, rownames(prot_mat)),
  , drop = FALSE
]

sum(rownames(protdata_filtered) != metadata_prot_filtered$id) == 0

rownames(protdata_filtered) <- rownames(metadata_prot_filtered)

# ---- Low count removal ----
prot_filtered <- low.removal(protdata_filtered, percent = 0.01)$data.filter

# ---- Scale features ----
prot_scaled <- scale(prot_filtered)%>%as.data.frame.matrix()


# ============================================================================
# KEEP ONLY INDIVIDUALS PRESENT IN BOTH PROTEOMICS & mOTUs
# ============================================================================
# Ensures all blocks have identical rows (samples),
# which is required for multi-block methods (e.g., DIABLO/tensor methods)
# ---------------------------------------------------------------------------

individuals_with_prot_mOTU <- intersect(
  metadata_prot_filtered$SampleID,
  metadata_mOTU_filtered$SampleID
)

metadata_common <- metadata_mOTU_filtered %>%
  filter(SampleID %in% individuals_with_prot_mOTU) %>%
  arrange(SampleID)

# ---- Final aligned data blocks ----
mOTU_clr_final <- mOTU_filtered.clr[
  match(rownames(metadata_common), rownames(mOTU_filtered.clr)),
]

ko_log_final <- ko_filtered.log[
  match(rownames(metadata_common), rownames(ko_filtered.log)),
]

prot_scale_final <- prot_scaled[
  match(rownames(metadata_common), rownames(prot_scaled)),
]

dim(metadata_common)
dim(mOTU_clr_final)
dim(ko_log_final)
dim(prot_scale_final)



# ---- Save preprocessed data ----

save(taxonomy_data, ko_info,
  mOTU_clr_final, ko_log_final,
  prot_scale_final,metadata_common,
  file = "Data_filtered/fmtStudy.RData")
