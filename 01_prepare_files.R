### setup ----------------------------------------------------------------------
source("config.R")
remotes::install_github("frvirard/datascrapR", auth_token = token)

# library ######################################################################
lib_lst <- c("dplyr",
             "purrr",
             "stringr",
             "reticulate")

invisible(purrr::walk(lib_lst, library, character.only = TRUE))
rm(lib_lst)

# functions ####################################################################
# variable #####################################################################

reticulate::use_python("~/Documents/Code/swansea_twistrand/renv/python/virtualenvs/renv-python-3.8/bin/python", required = TRUE)

### Download files from GDC ----------------------------------------------------
# I downloaded expression and mutations files from GDC
manifest <- datascrapR::gdc_download(project = "TCGA-LIHC",
                                     file_cat = "Simple Nucleotide Variation",
                                     sample_type = "Primary Tumor",
                                     download = FALSE)


