### setup ----------------------------------------------------------------------
# library ######################################################################
lib_lst <- c("dplyr",
             "purrr",
             "stringr",
             "reticulate")

invisible(purrr::walk(lib_lst, library, character.only = TRUE))
rm(lib_lst)

# functions ####################################################################
# variable #####################################################################
source("config.R")
remotes::install_github("frvirard/signatR", auth_token = token)
remotes::install_github("frvirard/datascrapR", auth_token = token)

maf_dr <- file.path(data_dr, "maf")

### Download files from GDC ----------------------------------------------------
#' I downloaded expression and mutations files from GDC
dataset_lst <- c("TCGA-LIHC")
dataset <- "TCGA-LIHC"
manifest <- datascrapR::gdc_download(project = dataset,
                                     file_cat = "Simple Nucleotide Variation",
                                     sample_type = "Primary Tumor",
                                     cache = file.path(file_dr, "tmp"),
                                     download = TRUE)

#' create maf ------------------------------------------------------------------
maf_file <- paste0(dataset, "_", manifest|>
                     pull(type)|>
                     unique(), ".maf")

if(!file.exists(file.path(maf_dr, maf_file))){
  message("-> Building maf file ...")
  maf <- pmap_dfr(manifest, function(id, file_name, ...){
    file.path(file_dr, "tmp", dataset, id, file_name)|>
      vroom::vroom(comment = "#", 
                   col_types = readr::cols(.default = "c"),
                   show_col_types = FALSE)
  }, .progress = TRUE)|>
    mutate(project = dataset, .before = 1)
  
  message("-> Saving maf file ...")
  vroom::vroom_write(maf, file.path(maf_dr, maf_file))
}
 
# create matrix ----------------------------------------------------------------
col_select <- c(Project = "project",
                Sample = "Tumor_Sample_Barcode",
                Genome = "NCBI_Build",
                chrom = "Chromosome",
                pos_start = "Start_Position",
                pos_end = "End_Position",
                ref = "Reference_Allele",
                alt = "Tumor_Seq_Allele2")

maf_lst <- list.files(file.path(maf_dr), 
                      pattern = ".maf", 
                      full.names = TRUE)

walk(maf_lst, signatR::sigpro_matrix, 
     matrix_path = file.path(data_dr, "sigprofiler/matrix"),
     col = col_select, clean = TRUE)

# extract signatures -----------------------------------------------------------
context_lst <- c("SBS96", "SBS288", "SBS384", "SBS1536")
walk(context_lst, function(c){
  set_lst <- dir(file.path(data_dr, "sigprofiler/matrix"))|>
    map(signatR::sigpro_extract, context = c, 
        solution_fit = TRUE, 
        compatibility = TRUE,
        sigpro_path = file.path(data_dr, "sigprofiler"))
})
