#' Correlation Hex Plot Function
#'
#' This function generates correlation hexagon plots for a given data frame.
#'
#' @description
#' The function calculates pairwise correlations between numeric columns in the input data frame `dt`,
#' and generates hexagon plots showing these correlations. It also includes facet-wrapped plots for each
#' unique pair of variables, with correlation coefficients annotated on the plots.
#'
#' @param dt A data frame containing the data to be analyzed.
#' @param id.col The name of the ID column in the data frame. Default is "protein".
#' @param cor.method The correlation method to use. Options are "pearson" (default), "kendall", or "spearman".
#' @param savefile Directory to save the generated plots. Default is "outputfile".
#' @param singleplotsize Size of each individual plot in inches (width, height). Default is c(3, 2.5).
#' @param facetplotsize Size of facet plot in inches (width, height). Default is c(9, 7.5).
#' @param bin Number of bins for hex plot. Default is 50.
#'
#' @return A list containing two elements:
#'   - facet_plt: Facet-wrapped correlation hexagon plot.
#'   - singleplot: List of individual correlation hexagon plots.
#'
#' @examples
#' dt <- data.frame(
#'   protein = c("A", "B", "C"),
#'   var1 = c(1, 2, 3),
#'   var2 = c(4, 5, 6),
#'   var3 = c(7, 8, 9)
#' )
#' cor_hex(dt)
#'
#' @importFrom dplyr select everything mutate filter distinct group_split left_join any_of
#' @importFrom tidyr pivot_longer
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_hex scale_fill_viridis_c theme_bw labs facet_wrap geom_text
#' @importFrom ggplot2 ggsave element_text element_rect element_line margin theme
#' @importFrom ggtext element_textbox
#' @importFrom purrr map map_chr walk2
#' @importFrom export graph2ppt
#' @importFrom grid unit
#' @importFrom stats cor
#'
#' @export
cor_hex=function(dt=dtplot,
                 id.col="protein",
                 cor.method=c("pearson", "kendall", "spearman")[1],
                 savefile="outputfile",
                 singleplotsize=c(3,2.5),#width height
                 facetplotsize=c(3*3,2.5*3),#width height
                 bin=50
){



  dt=dt |> dplyr::select(any_of(id.col),dplyr::everything())
  if (!is.data.frame(dt)) {
    dt <- as.data.frame(dt)
  }

  dt[-1] <- lapply(dt[-1], function(x) {
    if (is.numeric(x)) {
      x[is.infinite(x) | x == 0 | x == "" | is.nan(x)] <- NA
    }
    return(x)
  })


  dt[-1] <- lapply(dt[-1], function(x) as.numeric(as.character(x)))

  corvalue<- round(cor(as.matrix(dt[-1]), method = cor.method,use = "pairwise.complete.obs"), 3)



  corInofrs=reshape2::melt(corvalue, na.rm = TRUE) |> #dplyr::rename("xlab"="Var1","ylab"="Var2") |>
    dplyr::mutate(compars=paste0(Var1,"_vs_",Var2))



  dt_filtered=dt |>
    tidyr::pivot_longer(!id.col, names_to = "xlab", values_to = "x") |>
    dplyr::left_join(dt |> tidyr::pivot_longer(!id.col, names_to = "ylab", values_to = "y"),by = id.col)  |>
    dplyr::filter(xlab > ylab) |>
    dplyr::mutate(compars=paste0(xlab,"_vs_",ylab),
                  x=log10(x+1),
                  y=log10(y+1))  |>
    dplyr::left_join(corInofrs, by = "compars")


  #  dt_filtered <- dt.plot1  |> dplyr::filter(xlab > ylab)

  dt_label <- dt_filtered  |>
    dplyr::distinct(compars, .keep_all = TRUE)

  plt=ggplot2::ggplot(dt_filtered, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hex(bins = bin) +
    ggplot2::scale_fill_viridis_c(option = "viridis", guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = margin(t = 0.5)),
      axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = margin(r = 0.5)),
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      axis.line.x = ggplot2::element_line(color = "black"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = "italic"),
      panel.grid.major.y =ggplot2:: element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "right",
      strip.text = ggtext::element_textbox(
        size = 8, face = 'italic', color = "grey20",
        hjust = 0.5, halign = 0.5, r = grid::unit(5, "pt"),
        width = grid::unit(5.5, "npc"),
        padding = margin(3, 0, 3, 0),
        margin = margin(1, 1, 1, 1)
      ),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
      panel.spacing = grid::unit(1, 'lines')
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    ggplot2::facet_wrap(~ compars, scales = "free", drop = TRUE) +
    ggplot2::geom_text(data = dt_label,  ggplot2::aes(label = paste0("R^2 = ", round(value, 3))),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 2, color = "#9A2631")




  plots <- dt_filtered |>
    dplyr::group_split(compars) |>
    purrr::map(~ ggplot2::ggplot(.x, ggplot2::aes(x = x, y = y)) +
                 ggplot2::geom_hex(bins = bin) +
                 ggplot2::scale_fill_viridis_c(option = "viridis", guide = "none") +
                 ggplot2::theme_bw() +
                 ggplot2::theme(
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
                   axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = margin(t = 0.5)),
                   axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = margin(r = 0.5)),
                   axis.text = ggplot2::element_text(size = 8, color = "black"),
                   axis.line.x = ggplot2::element_line(color = "black"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = "italic"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   legend.position = "right",
                   strip.text = ggtext::element_textbox(
                     size = 8, face = 'bold', color = "grey20",
                     hjust = 0.5, halign = 0.5, r = grid::unit(5, "pt"),
                     width = grid::unit(5.5, "npc"),
                     padding = margin(3, 0, 3, 0),
                     margin = margin(1, 1, 1, 1)
                   ),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
                   panel.spacing = grid::unit(1, 'lines')
                 ) +
                 ggplot2::labs(
                   title = NULL,
                   x = .x$xlab[1],
                   y = .x$ylab[1]
                 ) +
                 ggplot2::geom_text(
                   ggplot2::aes(label = paste0("R^2 = ", round(.x$value[1], 3))),
                   x = -Inf, y = Inf,
                   hjust = -0.1, vjust = 1.1,
                   size = 2, color = "#9A2631"
                 )
    )



  if (savefile == "" || is.null(savefile)) {
    savepath <- getwd()
  } else {
    if (dir.exists(savefile)) {
      savepath <- savefile
    } else {
      dir.create(savefile, recursive = TRUE)
      savepath <- savefile
    }
  }



  plot_names <- dt_filtered |>
    dplyr::group_split(compars) |>
    purrr::map_chr(~ unique(.x$compars))

  #  savefile
  # singleplotsize
  #  savepath


  purrr::walk2(plots, plot_names, ~ {
    ggplot2::ggsave(filename = paste0(savepath,"/",.y, ".png"), plot = .x, width = singleplotsize[1], height = singleplotsize[2], dpi = 300)
    ggplot2::ggsave(filename = paste0(savepath,"/",.y, ".pdf"), plot = .x, width = singleplotsize[1], height = singleplotsize[2],device = "pdf")
    export::graph2ppt(x =.x, file =paste0(savepath,"/",.y, ".ppt"),
                      vector.graphic = TRUE,width = singleplotsize[1], height = singleplotsize[2], aspectr = sqrt(2), append = FALSE)
  })

  .save_zcp <- function(Fig,FigName,outputfile,widths,heights){
    Filepaths=paste0(outputfile,"/",FigName,c(".pdf",".png",".ppt"))
    ggplot2::ggsave(Filepaths[1], width =widths, plot=Fig,height = heights,device = "pdf")
    ggplot2::ggsave(Filepaths[2], width =widths,  plot=Fig,height = heights,device = "png")
    export::graph2ppt(x = Fig, file =Filepaths[3],
                      vector.graphic = TRUE, width =widths, height =heights, aspectr = sqrt(2), append = FALSE)
  }

  .save_zcp(Fig =plt,FigName = "facetwrapPlot_cor",outputfile =savepath,widths = facetplotsize[1],heights = facetplotsize[2])

  result=list(facet_plt=plt,singleplot=plots)
  return(result)

}
