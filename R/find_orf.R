#' Find ORFs (Open Reading Frames) in DNA Sequences
#'
#' This function identifies Open Reading Frames (ORFs) from a given DNA sequence or a set of sequences.
#' It uses specified start and stop codons to locate possible ORFs.
#'
#' @param seq_obj A DNA sequence or a DNAStringSet object. If passed as a file path, the file format should be FASTA.
#' @param start_codons A character vector of the start codons used to identify ORFs. Default is "ATG".
#' @param stop_codons A character vector of the stop codons used to terminate ORFs. Default is c("TAA", "TAG", "TGA").
#' @param minim_length Minimum length for an ORF (in nucleotides). Default is 75 nucleotides.
#' @param return_AA Logical. Whether to return the translated amino acid sequence of the ORF. Default is FALSE.
#' @param filter_orf Logical. If TRUE, nested ORFs are filtered to retain only the largest possible ORF between start and stop codon pairs. Default is TRUE.
#'
#' @return A data frame with information about the ORFs found, including start and stop positions, reading frame, strand orientation, ORF length, and optionally the amino acid sequence.
#'
#' The returned data frame contains the following columns:
#' \item{id}{The sequence identifier.}
#' \item{st_start}{The start position of the ORF based on start codon. (1-based offset)}
#' \item{st_stop}{The stop position of the ORF based on stop codon. (1-based offset)}
#' \item{st_seq}{The start codon sequence.}
#' \item{sp_start}{The start position of the stop codon.}
#' \item{sp_stop}{The stop position of the stop codon.}
#' \item{sp_seq}{The stop codon sequence.}
#' \item{seq_len}{The length of the ORF in nucleotides.}
#' \item{aa_len}{The length of the translated protein in amino acids.}
#' \item{frame}{The reading frame (1, 2, or 3).}
#' \item{sense}{The strand orientation, '+' for forward, '-' for reverse strand.}
#' \item{t_len}{The total length of the input sequence.}
#' \item{aa_seq}{(Optional) The amino acid sequence of the found ORF (if return_AA is TRUE).}
#' \item{orf}{The ORF identifier including its frame and strand orientation.}
#'
#' @details
#' This function looks for ORFs in both the forward and reverse strands of the input DNA sequence(s).
#' It searches for start codons on both strands provided (along with their reverse complements) and matches them
#' with possible stop codons to define ORFs. The user can specify a minimum length requirement and whether nested ORFs
#' should be filtered. Optionally, the corresponding amino acid sequences can be returned.
#'
#' @seealso \code{\link[Biostrings]{DNAStringSet}}, \code{\link[Biostrings]{translate}}, \code{\link[Biostrings]{reverseComplement}}
#'
#' @examples
#' \dontrun{# Example with a sequence string
#' seq <- "ATGAGAGTAGTAA"
#' find_orf(seq_obj = seq, return_AA = TRUE)
#'
#' # Example with fasta file input
#' # find_orf(seq_obj = "example.fasta", return_AA = TRUE)}
#'
#' @import dplyr
#' @importFrom Biostrings readDNAStringSet DNAStringSet reverseComplement translate
#' @importFrom stringr str_locate_all
#' @importFrom purrr map_df map2_chr
#'
#' @export
find_orf <- function(seq_obj = NULL,
                     start_codons = c("ATG"),
                     stop_codons = c("TAA","TAG","TGA"),
                     minim_length = 75,
                     return_AA = FALSE,
                     filter_orf = TRUE){
  # read seqeunce
  if(inherits(seq_obj,"DNAStringSet")){
    seq <- toupper(as.character(seq_obj))
  }else{
    seq <- toupper(as.character(Biostrings::readDNAStringSet(filepath = seq_obj,format = "fasta")))
  }


  # reverse_strand <- as.character(reverseComplement(DNAStringSet(x = seq)))
  # reverse_strand <- as.character(Biostrings::complement(DNAStringSet(x = seq)))

  strand <- c("+","-")

  start_codons <- c(start_codons,
                    as.character(Biostrings::reverseComplement(DNAStringSet(x = start_codons))))
  start_codons <- c(paste(start_codons[1:(length(start_codons)/2)],sep = "|",collapse = "|"),
                    paste(start_codons[(length(start_codons)/2 + 1):length(start_codons)],sep = "|",collapse = "|"))

  stop_codons <- c(stop_codons,
                   as.character(Biostrings::reverseComplement(DNAStringSet(x = stop_codons))))
  stop_codons <- c(paste(stop_codons[1:(length(stop_codons)/2)],sep = "|",collapse = "|"),
                   paste(stop_codons[(length(stop_codons)/2 + 1):length(stop_codons)],sep = "|",collapse = "|"))


  # loop for each sequencev
  # sq = 1
  purrr::map_df(seq_along(seq),function(sq){
    tmp_seq <- seq[sq]

    # lt = 1
    purrr::map_df(1:2,function(lt){
      seq_it <- tmp_seq

      start_codon_pos <- data.frame(stringr::str_locate_all(string = seq_it,pattern = start_codons[lt])[[1]])
      stop_codon_pos <- data.frame(stringr::str_locate_all(string = seq_it,pattern = stop_codons[lt])[[1]])

      stop_codon_pos_seq <- purrr::map2_chr(.x = stop_codon_pos$start,
                                            .y = stop_codon_pos$end,
                                            .f = substr,
                                            x = seq_it)

      start_codon_pos_seq <- purrr::map2_chr(.x = start_codon_pos$start,
                                             .y = start_codon_pos$end,
                                             .f = substr,
                                             x = seq_it)

      if(strand[lt] == "-"){
        tmp <- Biostrings::DNAStringSet(x = stop_codon_pos_seq)
        stop_codon_pos_seq <- as.character(Biostrings::reverseComplement(tmp))

        tmp2 <- Biostrings::DNAStringSet(x = start_codon_pos_seq)
        start_codon_pos_seq <- as.character(Biostrings::reverseComplement(tmp2))
      }

      stop_codon_pos$sp_seq <- stop_codon_pos_seq
      start_codon_pos$st_seq <- start_codon_pos_seq


      # x = 1
      purrr::map_df(1:nrow(start_codon_pos),function(x){
        df <- data.frame(st_start = start_codon_pos$start[x],
                         st_stop = start_codon_pos$end[x],
                         st_seq = as.vector(start_codon_pos[x,3]),
                         sp_start = as.vector(stop_codon_pos[,1]),
                         sp_stop = as.vector(stop_codon_pos[,2]),
                         sp_seq = as.vector(stop_codon_pos[,3]))

        if(strand[lt] == "-"){
          df <- df %>%
            dplyr::filter(sp_start < st_stop)
        }else{
          df <- df %>%
            dplyr::filter(sp_start > st_stop)
        }

        return(df)
      }) -> pos_info


      pos_info_ft <- pos_info %>%
        dplyr::mutate(seq_len = abs(sp_start - st_start)) %>%
        dplyr::mutate(frame = seq_len %% 3,aa_len = seq_len/3) %>%
        dplyr::filter(frame == 0) %>%
        dplyr::select(-frame)

      if(filter_orf){
        tmp1 <-
          pos_info_ft %>%
          dplyr::group_by(sp_start) %>%
          dplyr::filter(seq_len == max(seq_len)) %>%
          dplyr::group_by(st_start) %>%
          dplyr::filter(seq_len == min(seq_len))

        tmp2 <-
          pos_info_ft %>%
          dplyr::group_by(st_start) %>%
          dplyr::filter(seq_len == min(seq_len)) %>%
          dplyr::group_by(sp_start) %>%
          dplyr::filter(seq_len == max(seq_len))

        res <- unique(rbind(tmp1,tmp2)) %>%
          dplyr::ungroup()
      }else{
        res <- pos_info_ft
      }

      # add frame
      seq_id <- sapply(strsplit(x = names(seq_it),split = " "),"[",1)
      seq_length <- nchar(seq_it)

      if(strand[lt] == "-"){
        res <- res %>%
          dplyr::mutate(id = seq_id,t_len = seq_length) %>%
          dplyr::mutate(frame = (abs(st_stop - t_len) %% 3) + 1,
                        sense = strand[lt])
      }else{
        res <- res %>%
          dplyr::mutate(id = seq_id,t_len = seq_length) %>%
          dplyr::mutate(frame = ((st_start - 1) %% 3) + 1,
                        sense = strand[lt])
      }

      return(res)
    }) %>% dplyr::filter(seq_len + 3 >= minim_length)-> pos_info_all

    # ==========================================================================
    # getting amino acid
    # ==========================================================================
    if(return_AA){
      pos <- subset(pos_info_all,sense == "+")
      pos_seq_aa_ip <- Biostrings::translate(Biostrings::DNAStringSet(x = tmp_seq,start = pos$st_start,end = pos$sp_stop))
      pos$aa_seq <- as.character(pos_seq_aa_ip)

      neg <- subset(pos_info_all,sense == "-")
      neg_seq_aa_ip <- Biostrings::translate(Biostrings::reverseComplement(Biostrings::DNAStringSet(x = tmp_seq,start = neg$sp_start,end = neg$st_stop)))
      neg$aa_seq <- as.character(neg_seq_aa_ip)

      cmb <- rbind(pos,neg)
    }else{
      cmb <- pos_info_all
    }

    cmb$orf <- paste("orf",1:nrow(cmb),"|",cmb$frame,"|",cmb$sense,sep = " ")

    return(cmb)
  }) -> final_res


}
