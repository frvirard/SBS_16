### setup ----------------------------------------------------------------------
# library ######################################################################
lib_lst <- c("dplyr",
             "purrr",
             "readr",
             "stringr",
             "reticulate")

invisible(purrr::walk(lib_lst, library, character.only = TRUE))
rm(lib_lst)

# functions ####################################################################
# variable #####################################################################
source("config.R")
remotes::install_github("frvirard/signatR", auth_token = token, upgrade = "never")
remotes::install_github("frvirard/datascrapR", auth_token = token, upgrade = "never")
remotes::install_github("frvirard/phenotypeR", auth_token = token, upgrade = "never")

maf_dr <- file.path(data_dr, "maf")
expr_dr <- file.path(data_dr, "expr")

### TCGA-LIHC ------------------------------------------------------------------
#" download files
#' I downloaded mutations files from GDC (20240128)
dataset <- "TCGA-LIHC"
manifest <- datascrapR::gdc_download(project = dataset,
                                     file_cat = "Simple Nucleotide Variation",
                                     sample_type = "Primary Tumor",
                                     cache = file.path(file_dr, "tmp"),
                                     download = TRUE)

#' create maf ------------------------------------------------------------------
# maf_file <- paste0(dataset, "_", manifest|>
#                      pull(type)|>
#                      unique(), ".maf.gz")
maf_file <- paste0("simple_somatic_mutation.open.", dataset, ".maf.gz")

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
  str_subset("LICA-FR_LIRI-JP_TCGA-LIHC", negate = TRUE)|>
  walk(signatR::sigpro_matrix,
     matrix_path = file.path(data_dr, "sigprofiler/matrix"),
     col = col_select, clean = FALSE)

### LICA-FR and LIRI-JP --------------------------------------------------------
#' release 28 (downloaded 20240202) from ICGC
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

### CLCA_HCC -------------------------------------------------------------------
#' data from CLCA paper (PMID: 38355797) were provided by the authors
#' downloaded 20240402
dataset <- "CLCA_HCC"
vroom::vroom(file.path(file_dr, "CLCA_HCC", "CLCA_HCC_494.vcf"),
                   show_col_types = FALSE)|>
  mutate(project = dataset,
         assembly_version = "GRCh37", .after = 1)|>
  vroom::vroom_write(file.path(data_dr, "maf",
                               paste0("simple_somatic_mutation.open.",
                                      dataset, ".maf.gz")))

### All set --------------------------------------------------------------------
dataset_lst <- c("LICA-FR", "LIRI-JP", "TCGA-LIHC", "CLCA_HCC")
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

clca_col <- c("project",
              "Sample",
              "assembly_version",
              "Chr",
              "Start",
              "End",
              "Ref",
              "Alt")

rename_lst <- list('TCGA-LIHC' = tcga_col,
                   'LICA-FR' = icgc_col,
                   'LIRI-JP' = icgc_col,
                   'CLCA_HCC' = clca_col)

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
      str_subset(paste0("open.", dataset, ".maf.gz"))|>
      vroom::vroom(show_col_types = FALSE,
                   col_types = readr::cols(.default = "c"),
                   col_select = rename_lst[[dataset]])

    convert <- c("TCGA-LIHC")
    if(dataset %in% convert){
      # liftover ---------------------------------------------------------------
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
    dt <- dt|>
      dplyr::distinct()

    return(dt)
  })|>
    readr::type_convert()|>
    vroom::vroom_write(output)
}else{
  message(output, " already exists, skipping")
}

# Add annotation to all set ----------------------------------------------------
# Note that insertion have been removed during the process
# doesn't understand why there is duplicated lines
file_lst <- list.files(file.path(data_dr, "maf"), pattern = "maf")|>
  str_subset("simple_somatic_mutation.LICA-FR_LIRI-JP_TCGA-LIHC_CLCA_HCC.maf.gz")|>
  str_subset("anno", negate = TRUE)

walk(file_lst, function(file){
  prf <- str_remove(file, ".maf")|>
    str_remove(".gz")
  output <- file.path(data_dr, "maf", paste0(prf, "_anno.maf.gz"))

  if(file.exists(output)){
    message(basename(output), " file already present, skipping ...")
  }else{
    message("- creating annotation for ",basename(output))
    tmp_file <- paste0(prf, ".tmp")
    multianno_file <- paste0(prf, ".hg19_multianno.txt")

    # create a temp file
    vroom::vroom(file.path(data_dr, "maf", file),
                 col_select = c(chrom, pos_start, pos_end, ref, alt),
                 show_col_types = FALSE)|>
      vroom::vroom_write(file.path(data_dr, "maf", tmp_file), col_names = FALSE)

    # run annovar
    phenotypeR::run_annovar(input = tmp_file,
                genome = "hg19",
                db = "ensGene",
                op = "g",
                prefix = prf,
                annovar_tool = "/Users/francoisvirard/Documents/Code/annovar",
                path = file.path(data_dr, "maf"))

    # merge result
    init <- vroom::vroom(file.path(data_dr, "maf", file), show_col_types = FALSE,
                         col_types = cols(.default = col_character()))|>
      rename(Chr = chrom, Start = pos_start, End = pos_end, Ref = ref, Alt = alt)

    multi_annot <- vroom::vroom(file.path(data_dr, "maf", multianno_file),
                                show_col_types = FALSE,
                                col_types = cols(.default = col_character()))

    init |>
      left_join(multi_annot, relationship = "many-to-many",
                by = join_by(Chr, Start, End, Ref, Alt))|>
      vroom::vroom_write(output)

    # clean
    unlink(file.path(data_dr, "maf", multianno_file))
    unlink(file.path(data_dr, "maf", tmp_file))
    unlink(file.path(data_dr, "maf", paste0(prf, ".log")))
  }
})

# create matrix ----------------------------------------------------------------
names(col_names) <- col_names
signatR::sigpro_matrix(output,
                       matrix_path = file.path(data_dr, "sigprofiler/matrix"),
                       col = col_names, clean = FALSE)

### extract signatures ---------------------------------------------------------
context_lst <- c("SBS96", "SBS288", "SBS384", "SBS1536")
walk(context_lst, function(c){
  set_lst <- dir(file.path(data_dr, "sigprofiler/matrix"))|>
    map(signatR::sigpro_extract, context = c,
        n_max = 20,
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

  if(!file.exists(file.path(file_dr, expr_file))){
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
    vroom::vroom_write(expr, file.path(file_dr, expr_file))
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
           aliquot = aliquots_submitter_id,
           sample_type = samples_tissue_type)
})

## LIRI
LIRI_specimen <- vroom::vroom(file.path(file_dr, "LIRI-JP", "specimen.LIRI-JP.tsv.gz"),
                              show_col_types = FALSE,
                              col_select = c(aliquot = icgc_specimen_id,
                                             sample_type = specimen_type))|>
  mutate(sample_type = case_when(str_detect(sample_type, "Normal") ~ "Normal",
                                 str_detect(sample_type, "tumour") ~ "Tumor"))

LIRI_sample <- vroom::vroom(file.path(file_dr, "LIRI-JP", "sample.LIRI-JP.tsv.gz"),
                            show_col_types = FALSE)|>
  select(project = project_code,
         patient = icgc_donor_id,
         sample = icgc_sample_id,
         aliquot = icgc_specimen_id)|>
  left_join(LIRI_specimen, join_by(aliquot))

## LICA
LICA_specimen <- vroom::vroom(file.path(file_dr, "LICA-FR", "specimen.LICA-FR.tsv.gz"),
                              show_col_types = FALSE,
                              col_select = c(aliquot = icgc_specimen_id,
                                             sample_type = specimen_type))|>
  mutate(sample_type = case_when(str_detect(sample_type, "Normal") ~ "Normal",
                                 str_detect(sample_type, "tumour") ~ "Tumor"))

LICA_sample <- vroom::vroom(file.path(file_dr, "LICA-FR", "sample.LICA-FR.tsv.gz"),
                            col_select = c(project = project_code,
                            patient = icgc_donor_id,
                            sample = icgc_sample_id,
                            aliquot = icgc_specimen_id),
                            show_col_types = FALSE)|>
  left_join(LICA_specimen, join_by(aliquot))

biospecimen_id <- rbind(TCGA_sample, LIRI_sample, LICA_sample)

vroom::vroom_write(biospecimen_id, file.path(data_dr, "biospecimen_id.txt"))

#' get cohort clinical data ----------------------------------------------------

#' LICA-FR
#' supplementary from
#'  DNA Methylation Signatures Reveal the Diversity of Processes Remodeling
#'  Hepatocellular Carcinoma Methylomes
#' Léa Meunier, Théo Z. Hirsch, Stefano Caruso, Sandrine Imbeaud, Quentin Bayard,
#' Amélie Roehrig, Gabrielle Couchy, Jean-Charles Nault, Josep M. Llovet,
#' Jean-Frédéric Blanc
#' https://doi.org/10.1002/hep.31796
#'
# 13 samples were absent from ICGC cohort
LICA_specimen_id <- vroom::vroom(file.path(file_dr, "LICA-FR", "specimen.LICA-FR.tsv.gz"),
                              show_col_types = FALSE,
                              col_select = c(aliquot = icgc_specimen_id,
                                             submitter_id = submitted_specimen_id))

lica_clin <- readxl::read_xlsx(file.path(file_dr, "LICA-FR",
                                         "hep31796-sup-0001-tables1.xlsx"), skip=1)|>
  select(submitter_id = "Sample",
         gender = "Gender",
         age = "Age at sampling",
         alcohol = "Alcohol intake",
         HBV = "Hepatitis B",
         HCV = "Hepatitis C",
         tobacco = "Tobacco")|>
  mutate(HBV = if_else(HBV == "yes", "HBV", "NB"),
         HCV = if_else(HCV == "yes", "HCV", "NC"))|>
  unite(virus, c(HBV, HCV))|>
  mutate(virus = if_else(virus == "NB_NC", "NBNC", virus))|>
  mutate(virus = str_remove(virus, "NB_|_NC"))|>
  left_join(LICA_specimen_id)|>
  select(-submitter_id)|>
  relocate(aliquot, .before=1)|>
  filter(!is.na(aliquot))

#' LIRI-JP
#' supplementary from :
#' Whole-genome mutational landscape and characterization of noncoding and
#' structural mutations in liver cancer
#' Akihiro Fujimoto, Mayuko Furuta, Yasushi Totoki, Tatsuhiko Tsunoda,
#' Mamoru Kato, Yuichi Shiraishi, Hiroko Tanaka, Hiroaki Taniguchi,
#' Yoshiiku Kawakami, Masaki Ueno, Kunihito Gotoh, Shun-ichi Ariizumi,
#' Christopher P Wardell, Shinya Hayami, Toru Nakamura, Hiroshi Aikata,
#' Koji Arihiro, Keith A Boroevich, Tetsuo Abe, Kaoru Nakano,
#' Kazuhiro Maejima, Aya Sasaki-Oku, Ayako Ohsawa, Tetsuo Shibuya, Hidewaki Nakagawa
#'
#' Nature Genetics volume 48, pages 500–509 (2016)
LIRI_specimen_id <- vroom::vroom(file.path(file_dr, "LIRI-JP", "specimen.LIRI-JP.tsv.gz"),
                                 show_col_types = FALSE,
                                 col_select = c(aliquot = icgc_specimen_id,
                                                submitter_id = submitted_specimen_id,
                                                type = specimen_type))|>
  filter(type == "Primary tumour - solid tissue")

liri_clin <- readxl::read_xlsx(file.path(file_dr, "LIRI-JP",
"03_differential_expression.xlsx"), sheet = 2, skip = 1)|>
  select(submitter_id = "ID",
         gender = "Gender",
         age = "Age",
         virus = "Virus infection",
         alcohol = "Alcohol intakee",
         tobacco = "Smoking")|>
  mutate(alcohol = case_match(alcohol,
                              "0" ~ "no",
                             c("1", "2") ~ "yes",
                             .default = NA))|>
  mutate(tobacco = case_match(tobacco,
                              "0" ~ "no",
                              c("1", "2") ~ "yes",
                              .default = NA))|>
  left_join(LIRI_specimen_id)
