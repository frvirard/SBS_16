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
expr_dr <- file.path(data_dr, "expr")

#' download mutation files -----------------------------------------------------
#' mutations files were downloaded from GDC 20230127
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

#' create matrix ---------------------------------------------------------------
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

#' extract signatures ----------------------------------------------------------
context_lst <- c("SBS96", "SBS288", "SBS384", "SBS1536")
walk(context_lst, function(c){
  set_lst <- dir(file.path(data_dr, "sigprofiler/matrix"))|>
    map(signatR::sigpro_extract, context = c,
        solution_fit = TRUE,
        compatibility = TRUE,
        sigpro_path = file.path(data_dr, "sigprofiler"))
})

#' download expression files ---------------------------------------------------
#' expression files were downloaded from GDC 20230128
#' need to aggregate tumoral and tumoral expr for differential exp
dataset_lst <- c("TCGA-LIHC")
dataset <- "TCGA-LIHC"
manifest <- datascrapR::gdc_download(project = dataset,
                                     file_cat = "Transcriptome Profiling", 
                                     feature_type = "Gene Expression Quantification",
                                     sample_type = "Primary Tumor",
                                     cache = file.path(file_dr, "tmp"),
                                     download = TRUE)

#' build expression file -------------------------------------------------------
#' need to add sample id and sample type
expr_file <- paste0(dataset, "_", manifest|>
                     pull(type)|>
                     unique(), ".tsv")

if(!file.exists(file.path(expr_dr, expr_file))){
  message("-> Building expression file ...")
  expr <- pmap_dfr(manifest, function(id, file_name, ...){
    file.path(file_dr, "tmp", dataset, id, file_name)|>
      vroom::vroom(comment = "#", 
                   col_types = readr::cols(.default = "c"),
                   show_col_types = FALSE)
  }, .progress = TRUE)|>
    mutate(project = dataset, .before = 1)
  
  message("-> Saving tsv file ...")
  vroom::vroom_write(expr, file.path(expr_dr, expr_file))
}
