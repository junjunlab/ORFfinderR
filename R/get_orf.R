globalVariables(c("end", "frame", "sense", "sp_start", "st_start", "st_stop", "start", "t_len",
                  "orf", "sp_stop"))

#' Get Open Reading Frames (ORFs) from a DNA sequence
#'
#' This function identifies potential open reading frames (ORFs) within a given DNA sequence
#' based on specified start and stop codons. Optionally, it can filter out nested ORFs and/or
#' trim trailing sequences.
#'
#' @param sequence `character`. The input DNA sequence as a string. Default is `NULL`.
#' @param translation_table `integer`. Genetic code translation table ID. Default is `1` (Standard genetic code).
#' @param minimum_length `integer`. Minimum ORF length in base pairs (bp). Default is `75`.
#' @param start_codons `character vector`. Start codons used to identify ORF initiation points. Default is `c("ATG")`.
#' @param stop_codons `character vector`. Stop codons used to identify ORF termination points. Default is `c("TAA", "TAG", "TGA")`.
#' @param remove_nested `logical`. If `TRUE`, nested ORFs (those fully contained within a larger ORF) will be removed. Default is `FALSE`.
#' @param trim_trailing `logical`. If `TRUE`, trims trailing sequences at the ORF boundaries for incomplete codons. Default is `FALSE`.
#' @param output_file The output file name.
#' @param pythonPath Python path.
#' @return A data frame.
#'
#' @details
#' The function uses the Python Biopython package (via the `reticulate` package) to carry out the ORF finding
#' and protein translation. If Python or the required Biopython package is not available on the system, the function
#' attempts to install Biopython automatically. The function works by running a Python script (`orf_finder.py`)
#' which looks for potential ORFs in the given input sequence.
#'
#'
#'
#' @importFrom reticulate py_available py_install use_python py_module_available py_config
#' @importFrom utils read.delim
#'
#' @export
get_orf <- function(sequence = NULL,
                    translation_table = 1,
                    minimum_length = 75,
                    start_codons = c("ATG"),
                    stop_codons = c("TAA", "TAG", "TGA"),
                    remove_nested = FALSE,
                    trim_trailing = FALSE,
                    output_file = NULL,
                    pythonPath = NULL){
  reticulate::py_config()
  if(reticulate::py_available() == FALSE){
    message("Please install python first!")
  }else{
    if(!is.null(pythonPath)){
      reticulate::use_python(pythonPath)
    }

    # check modual
    if (!reticulate::py_module_available("biopython")) {
      reticulate::py_install("biopython")
    }

    # run code
    pyscript.path = system.file("extdata", "orf_finder.py", package = "ORFfinderR")
    reticulate::source_python(pyscript.path)


    reticulate::py$batch_get_orf(sequence = sequence,
                                 minimum_length = reticulate::r_to_py(minimum_length),
                                 start_codons = reticulate::r_to_py(start_codons),
                                 stop_codons = reticulate::r_to_py(stop_codons),
                                 remove_nested = reticulate::r_to_py(remove_nested),
                                 trim_trailing = reticulate::r_to_py(trim_trailing),
                                 translation_table = reticulate::r_to_py(translation_table),
                                 output_file = output_file)

  }

  # load into R
  df <- read.delim(output_file,header = FALSE)
  colnames(df) <- c("id","t_len","start","end","frame","sense","length","AA_seq")

  df <- df %>%
    dplyr::mutate(AA_len = length/3 - 1,
                  start = dplyr::if_else(sense == "-",start - 1,start),
                  end = dplyr::if_else(sense == "+",end - 1,end))

  df$orf <- paste("orf",1:nrow(df),"|",df$frame,"|",df$sense,sep = " ")

  return(df)
}
