#' Plot ORFs (Open Reading Frames) for Each Transcript
#'
#' This function generates visualizations of Open Reading Frames (ORFs) across transcript sequences using `ggplot2`.
#' For each transcript (specified by the "id" column), ORFs are shown as colored rectangles based on their strand direction ('+' or '-').
#' A separate plot is created for each transcript, and all plots are combined together.
#'
#' @param data A data frame containing information about the ORFs. It must contain the columns `start`, `end`, `sense`, `orf`, and `id` or alternatively the columns `st_start` and `sp_stop`.
#' @param trans_col A color string representing the color for the transcript line in the plot. Default is `#996600`.
#'
#'
#' @return A combined ggplot object using `cowplot::plot_grid` to visualize the ORFs for all transcripts.
#' The resulting plot will display individual ORFs as colored rectangles along the transcript axis, with different colors for `+` and `-` strands.
#'
#'
#' @import ggplot2
#' @import ggside
#' @importFrom cowplot plot_grid
#' @importFrom dplyr filter select mutate group_by
#' @export
orf_plot <- function(data = NULL,trans_col = "#996600"){
  if("start" %in% colnames(data) & "end" %in% colnames(data)){
    layer <- geom_rect(aes(xmin = start,xmax = end,ymin = -0.5,ymax = 0.5,
                           fill = sense),
                       color = "black")
  }else{
    layer <- geom_rect(aes(xmin = st_start,xmax = sp_stop,ymin = -0.5,ymax = 0.5,
                           fill = sense),
                       color = "black")
  }

  # plot
  lapply(unique(data$id),function(x){
    tmp_data <- subset(data,id == x)

    struc <- tmp_data %>%
      dplyr::select(id,orf,t_len) %>%
      unique()

    ggplot(tmp_data) +
      geom_point(aes(x = 0,y = 0),color = NA) +
      layer +
      ggside::geom_xsidesegment(data = struc,
                                mapping = aes(x = 1,y = 1,xend = t_len,yend = 1),
                                linewidth = 3,color = trans_col) +
      # facet_wrap(~orf,ncol = 1,switch = "y") +
      facet_grid(orf~id,switch = "y",scales = "free_x") +
      # ggh4x::facet_grid2(orf~id,switch = "y",scales = "free_x",render_empty = FALSE) +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold",size = rel(1)),
            strip.text.y.left = element_text(angle = 0),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            # panel.spacing.y = unit(0,"mm"),
            ggside.panel.background = element_blank(),
            ggside.panel.border = element_blank(),
            axis.text = element_text(colour = "black"),
            axis.ticks.y = element_blank()) +
      ylab("") + xlab("Transcript position") +
      scale_fill_manual(values = c("+" = "#CC99CC","-" = "#99CCCC"))
  }) -> plist

  # do.call("+",plist)
  cowplot::plot_grid(plotlist = plist)
}
