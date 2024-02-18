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

### Download TCGA-LIHC from GDC ------------------------------------------------
#' I downloaded mutations files from GDC (20240128)
dataset <- "TCGA-LIHC"
manifest <- datascrapR::gdc_download(project = dataset,
                                     file_cat = "Simple Nucleotide Variation",
                                     sample_type = "Primary Tumor",
                                     cache = file.path(file_dr, "tmp"),
                                     download = TRUE)

#' create maf ------------------------------------------------------------------
maf_file <- paste0(dataset, "_", manifest|>
                     pull(type)|>
                     unique(), ".maf.gz")

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

list.files(file.path(maf_dr),
                      pattern = ".maf.gz",
                      full.names = TRUE)|>
  str_subset("TCGA-LIHC")|>
  walk(signatR::sigpro_matrix,
     matrix_path = file.path(data_dr, "sigprofiler/matrix"),
     col = col_select, clean = TRUE)

### LICA-FR and LIRI-JP from ICGC ----------------------------------------------
#' release 28 (downloaded 20240202)
#' warning a correction in November 26, 2019 was uploaded on ICGC for
#' LICA-FR expression files
dataset_lst <- c("LICA-FR", "LIRI-JP")

walk(dataset_lst, function(dataset){
  message("oo Processing ", dataset)
  vroom::vroom(file.path(file_dr, dataset,
                         paste0("simple_somatic_mutation.open.", dataset, ".tsv.gz")),
               show_col_types = FALSE)|>
    vroom::vroom_write(file.path(data_dr, "maf",
                                 paste0("simple_somatic_mutation.open.", dataset, ".maf.gz")))
  message("Done")
})

# All set ----------------------------------------------------------------------
dataset_lst <- c("LICA-FR", "LIRI-JP", "TCGA-LIHC")
col_names <- c("Project", "Sample", "Genome", "chrom", "pos_start", "pos_end", "ref", "alt")
tcga_col <- c("project",
              "Tumor_Sample_Barcode",
              "NCBI_Build",
              "Chromosome",
              "Start_Position",
              "End_Position",
              "Reference_Allele",
              "Tumor_Seq_Allele2")

icgc_col <- c("project_code",
              "icgc_sample_id",
              "assembly_version",
              "chromosome",
              "chromosome_start",
              "chromosome_end",
              "mutated_from_allele",
              "mutated_to_allele")

rename_lst <- list('TCGA-LIHC' = tcga_col,
                   'LICA-FR' = icgc_col,
                   'LIRI-JP' = icgc_col)

rename_lst <- map(rename_lst, function(x){
  names(x) <- col_names
  return(x)
})

output <- file.path(data_dr, "maf",
                    paste0("simple_somatic_mutation.",
                           paste0(dataset_lst, collapse = "_"),
                           ".maf.gz"))
if(!file.exists(output)){
  map_dfr(dataset_lst, function(dataset){
    message("oo Processing ", dataset)
    dt <- list.files(file.path(maf_dr),
               pattern = ".maf.gz",
               full.names = TRUE)|>
      str_subset(dataset)|>
      vroom::vroom(show_col_types = FALSE,
                   col_types = readr::cols(.default = "c"),
                   col_select = rename_lst[[dataset]])

    if(dataset == "TCGA-LIHC"){
      # liftover ---------------------------------------------------------------------
      message("converting to GRCh37 genome ...")
      chain <- rtracklayer::import.chain(file.path(file_dr,"hg38ToHg19.over.chain"))

      dt <- dt |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                                 seqnames.field = "chrom",
                                                 start.field="pos_start",
                                                 end.field = "pos_end",
                                                 strand.field = "Strand")|>
        rtracklayer::liftOver(chain)|>
        as_tibble()|>
        mutate(across(everything(), as.character))|>
        mutate(Genome = "GRCh37")|>
        select(-c(group, group_name, width, strand))|>
        rename(chrom = seqnames, pos_start = start, pos_end = end)

      message("done")
    }
    return(dt)
  })|>
    readr::type_convert()|>
    vroom::vroom_write(output)
}else{
  message(output, " already exists, skipping")
}

# create matrix ----------------------------------------------------------------
names(col_names) <- col_names
signatR::sigpro_matrix(output,
                       matrix_path = file.path(data_dr, "sigprofiler/matrix"),
                       col = col_names, clean = TRUE)

### extract signatures ---------------------------------------------------------
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
dataset_lst <- c("TCGA-LIHC", "TCGA-HNSC")|>
  set_names()

walk(dataset_lst, function(dataset){
  manifest <- datascrapR::gdc_download(project = dataset,
                                       file_cat = "Transcriptome Profiling",
                                       feature_type = "Gene Expression Quantification",
                                       sample_type = "Primary Tumor",
                                       cache = file.path(file_dr, "tmp"),
                                       download = TRUE)

  #' build expression file -----------------------------------------------------
  #' need to add sample id and sample type
  expr_file <- paste0(dataset, "_", manifest|>
                        pull(type)|>
                        unique(), ".tsv.gz")

  if(!file.exists(file.path(expr_dr, expr_file))){
    message("-> Building expression file ...")
    expr <- pmap_dfr(manifest, function(id, file_name, submitter_id, ...){
      file.path(file_dr, "tmp", dataset, id, file_name)|>
        vroom::vroom(comment = "#",
                     col_types = readr::cols(.default = "c"),
                     show_col_types = FALSE)|>
        mutate(sample = submitter_id)
    }, .progress = TRUE)|>
      mutate(project = dataset, .before = 1)

    message("-> Saving tsv file ...")
    vroom::vroom_write(expr, file.path(expr_dr, expr_file))
  }
})

#' get cohort sample ID --------------------------------------------------------
#' patient>sample>aliquot

## TCGA
TCGA_sample <- map_dfr(dataset_lst, function(dataset){
  datascrapR::gdc_clinical(dataset)|>
    select(project = project_id,
           patient = submitter_id,
           sample = samples_submitter_id,
           aliquots = aliquots_submitter_id,
           sample_type = samples_tissue_type)
})

## LIRI
LIRI_specimen <- vroom::vroom(file.path(file_dr, "LIRI-JP", "specimen.LIRI-JP.tsv.gz"),
                              show_col_types = FALSE,
                              col_select = c(aliquots = icgc_specimen_id,
                                             sample_type = specimen_type))|>
  mutate(sample_type = case_when(str_detect(sample_type, "Normal") ~ "Normal",
                                 str_detect(sample_type, "tumour") ~ "Tumor"))

LIRI_sample <- vroom::vroom(file.path(file_dr, "LIRI-JP", "sample.LIRI-JP.tsv.gz"),
                            show_col_types = FALSE)|>
  select(project = project_code,
         patient = icgc_donor_id,
         sample = icgc_sample_id,
         aliquots = icgc_specimen_id)|>
  left_join(LIRI_specimen, join_by(aliquots))

## LICA
LICA_specimen <- vroom::vroom(file.path(file_dr, "LICA-FR", "specimen.LICA-FR.tsv.gz"),
                              show_col_types = FALSE,
                              col_select = c(aliquots = icgc_specimen_id,
                                             sample_type = specimen_type))|>
  mutate(sample_type = case_when(str_detect(sample_type, "Normal") ~ "Normal",
                                 str_detect(sample_type, "tumour") ~ "Tumor"))

LICA_sample <- vroom::vroom(file.path(file_dr, "LICA-FR", "sample.LICA-FR.tsv.gz"),
                            col_select = c(project = project_code,
                            patient = icgc_donor_id,
                            sample = icgc_sample_id,
                            aliquots = icgc_specimen_id),
                            show_col_types = FALSE)|>
  left_join(LICA_specimen, join_by(aliquots))

biospecimen_id <- rbind(TCGA_sample, LIRI_sample, LICA_sample)

vroom::vroom_write(biospecimen_id, file.path(data_dr, "biospecimen_id.txt"))

#' get cohort clinical data ----------------------------------------------------
patient_data <- vroom::vroom("data/files/LICA-FR/donor.LICA-FR.tsv.gz")|>
  select(project = project_code,
         patient = icgc_donor_id,
         genre = donor_sex,
         age = donor_age_at_enrollment)

clinical <- vroom::vroom("data/files/LICA-FR/specimen.LICA-FR.tsv.gz")

patient_data$patient |> unique()
