#' Correlation Hex Plot Function
#'
#' This function generates correlation hexagon or density point plots for a given data frame,
#' including pairwise correlation calculations and visualization.
#'
#' @description
#' The function calculates pairwise correlations between numeric columns in the input data frame `dt`,
#' and generates either hexagon plots or density-based point plots to visualize these correlations.
#' Facet-wrapped plots are created for each unique pair of variables, with correlation coefficients annotated.
#' It also supports kernel density estimation for creating density-based point plots.
#'
#' @param dt A data frame containing the data to be analyzed.
#' @param id.col The name of the ID column in the data frame. Default is "protein".
#' @param cor.method The correlation method to use. Options are "pearson" (default), "kendall", or "spearman".
#' @param savefile Directory to save the generated plots. Default is "outputfile".
#' @param singleplotsize Size of each individual plot in inches (width, height). Default is c(3, 2.5).
#' @param facetplotsize Size of facet plot in inches (width, height). Default is c(9, 7.5).
#' @param bin Number of bins for hex plot. Default is 50.
#' @param logL Logical. Whether to log-transform the data (log10). Default is TRUE.
#' @param plottype Character. Type of plot to generate: "hex" for hexagon plots, "point" for density point plots. Default is "point".
#' @param pointcolor Character. Color scheme for density point plots. Options include "A", "B", "C", "D", "E". Default is "C".
#' @param pointsize Numeric. Size of points in density point plots. Default is 0.5.
#' @param kde2d.n Numeric. Grid size for kernel density estimation using `MASS::kde2d`. Default is 50.
#' @param smooth.color Color of the smoothing line in scatter and hexbin plots. Default is "blue".
#' @param smooth.width Width of the smoothing line in scatter and hexbin plots. Default is `0.8`.
#' @param savePPT Logical, whether to save the plot as a PowerPoint file. Default is `FALSE`.
#'
#' @return A list containing two elements:
#'   - facet_plt: Facet-wrapped plot (hexagon or point density).
#'   - singleplot: List of individual plots (hexagon or point density).
#'
#' @details
#' The function uses the following key dependencies:
#' - **MASS::kde2d**: Performs 2D kernel density estimation for density-based point plots.
#' - **spatstat.geom::as.im** and **spatstat.geom::interp.im**: Converts the output of `MASS::kde2d` into an interpolated spatial image for density calculation.
#' - **ggplot2**: Creates the hexagon or point density plots and handles all visualization components.
#' - **dplyr** and **tidyr**: Perform data wrangling, including pivoting and filtering.
#' - **reshape2::melt**: Converts correlation matrices into long format for easier processing.
#' - **purrr**: Used to map over combinations of data for generating individual plots.
#'
#' @examples
#' # Load the demo data included in the package
#' data("demo", package = "corhex")
#'
#' # Run the cor_hex function on the demo data
#' cor_hex(dt=demo,
#' id.col="Index",
#' cor.method=c("pearson", "kendall", "spearman")[3],
#' savefile="outputfile",
#' singleplotsize=c(3,2.5),
#' facetplotsize=c(3*3,2.5*3),
#' bin=50,
#' logL=FALSE,
#' pointcolor=c("A","B","C","D","E")[3],
#' pointsize=0.5,
#' kde2d.n=50,
#' plottype=c("hex","point")[2])
#'
#' @importFrom dplyr select everything mutate filter distinct group_split left_join any_of
#' @importFrom tidyr pivot_longer
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_hex scale_fill_viridis_c theme_bw labs facet_wrap geom_text
#' @importFrom ggplot2 ggsave element_text element_rect element_line margin theme
#' @importFrom ggtext element_textbox
#' @importFrom purrr map map_chr walk2
#' @importFrom export graph2ppt
#' @importFrom MASS kde2d
#' @importFrom spatstat.geom as.im interp.im
#' @importFrom grid unit
#' @importFrom stats cor
#'
#' @export
cor_hex=function(dt=dtplot,
                 id.col="Index",
                 cor.method=c("pearson", "kendall", "spearman")[1],
                 savefile="outputfile",
                 singleplotsize=c(3,2.5),#width height
                 facetplotsize=c(3*3,2.5*3),#width height
                 bin=50,
                 logL=TRUE,
                 plottype=c("hex","point")[2],
                 pointcolor=c("A","B","C","D","E")[3],
                 pointsize=0.5,
                 kde2d.n=50,
                 smooth.color="red",
                 smooth.width=0.8,
                 savePPT=FALSE
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


  column_combinations <- combn(names(dt)[-1], 2, simplify = FALSE)


  dt_filtered <- purrr::map_dfr(column_combinations, function(cols) {
    dt |>
      dplyr::select(any_of(id.col), all_of(cols)) |>
      dplyr::rename(x = !!cols[1], y = !!cols[2]) |>
      dplyr::mutate(
        compars = paste0(cols[1], "_vs_", cols[2]),
        x = if (!logL) log10(x + 1) else x,
        y = if (!logL) log10(y + 1) else y
      )
  }) |>
    dplyr::filter(!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)) |>
    dplyr::left_join(corInofrs |> dplyr::select(value, compars), by = "compars")



  # geompoint data


  kde_result <- MASS::kde2d(dt_filtered$x, dt_filtered$y, n = kde2d.n)
  dt_point_plot <- dt_filtered |>
  dplyr::mutate(
    density = spatstat.geom::interp.im(
      spatstat.geom::as.im(kde_result$z, xcol = kde_result$x, yrow = kde_result$y),
      x, y
    )
  )



  dt_label <- dt_filtered  |>
    dplyr::distinct(compars, .keep_all = TRUE)

switch (plottype,
  "hex" = {

    plt=ggplot2::ggplot(dt_filtered, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_hex(bins = bin) +
      ggplot2::geom_smooth(color=smooth.color,linewidth=smooth.width,method="lm",formula = y ~ x,se = TRUE, level = 0.95)+  #     smooth.color="blue"  smooth.width=0.8 =========
      #ggpubr::stat_cor(method = cor.method, ggplot2::aes(x =x, y =y),size=2.5,color="#9A2631")+
      ggplot2::scale_fill_viridis_c(option = "viridis", guide = "none") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
        axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(t = 0.5)),
        axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(r = 0.5)),
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
          padding = ggplot2::margin(3, 0, 3, 0),
          margin = ggplot2::margin(1, 1, 1, 1)
        ),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
        panel.spacing = grid::unit(1, 'lines')
      ) +
      ggplot2::labs(
        title = NULL,
        x = NULL,
        y = NULL
      ) +
      ggplot2::facet_wrap(~ compars, scales = "free", drop = TRUE) +
      ggplot2::geom_text(data = dt_label,  ggplot2::aes(label = paste0("R = ", round(value, 3))),
                         x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 2.5, color = "#9A2631")




    plots <- dt_filtered |>
      dplyr::group_split(compars) |>
      purrr::map(~ ggplot2::ggplot(.x, ggplot2::aes(x = x, y = y)) +
                   ggplot2::geom_hex(bins = bin) +
                   ggplot2::geom_smooth(color=smooth.color,linewidth=smooth.width,method="lm",formula = y ~ x,se = TRUE, level = 0.95)+
                   ggplot2::scale_fill_viridis_c(option = "viridis", guide = "none") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(
                     plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                     panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                     plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                     axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(t = 0.5)),
                     axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(r = 0.5)),
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
                       padding = ggplot2::margin(3, 0, 3, 0),
                       margin = ggplot2::margin(1, 1, 1, 1)
                     ),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
                     panel.spacing = grid::unit(1, 'lines')
                   ) +
                   ggplot2::labs(
                     title = NULL,
                     x = .x$xlab[1],
                     y = .x$ylab[1]
                   ) +
                   ggplot2::geom_text( data = dt_label |>
                                         dplyr::filter(compars == .x$compars[1]),
                                       ggplot2::aes(label = paste0("R = ", round(value, 3))),
                                       x = -Inf, y = Inf,
                                       hjust = -0.1, vjust = 1.1,
                                       size = 2.5, color = "#9A2631"
                   )
      )



  },
  "point"={

    plt=ggplot2::ggplot(dt_point_plot, ggplot2::aes(x = x, y = y,color=density)) +
      ggplot2::geom_point(alpha = 0.8, size = pointsize)+
      ggplot2::geom_smooth(color=smooth.color,linewidth=smooth.width,method="lm",formula = y ~ x,se = TRUE, level = 0.95)+
      ggplot2::scale_color_viridis_c(option = pointcolor, guide = "none") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
        axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(t = 0.5)),
        axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(r = 0.5)),
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
          padding = ggplot2::margin(3, 0, 3, 0),
          margin = ggplot2::margin(1, 1, 1, 1)
        ),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
        panel.spacing = grid::unit(1, 'lines')
      ) +
      ggplot2::labs(
        title = NULL,
        x = NULL,
        y = NULL
      ) +
      ggplot2::facet_wrap(~ compars, scales = "free", drop = TRUE) +
      ggplot2::geom_text(data = dt_label,  ggplot2::aes(label = paste0("R = ", round(value, 3))),
                         x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 2.5, color = "#9A2631")

    #single fig
    plots <- dt_point_plot |>
      dplyr::group_split(compars) |>
      purrr::map(~ ggplot2::ggplot(.x, ggplot2::aes(x = x, y = y,color=density)) +
                   ggplot2::geom_point(alpha = 0.8, size = pointsize) +
                   ggplot2::geom_smooth(color=smooth.color,linewidth=smooth.width,method="lm",formula = y ~ x,se = TRUE, level = 0.95)+
                   ggplot2::scale_color_viridis_c(option = pointcolor, guide = "none") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(
                     plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                     panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                     plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                     axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(t = 0.5)),
                     axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(r = 0.5)),
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
                       padding = ggplot2::margin(3, 0, 3, 0),
                       margin = ggplot2::margin(1, 1, 1, 1)
                     ),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
                     panel.spacing = grid::unit(1, 'lines')
                   ) +
                   ggplot2::labs(
                     title = NULL,
                     x = .x$xlab[1],
                     y = .x$ylab[1]
                   ) +
                   ggplot2::geom_text( data = dt_label |>
                                         dplyr::filter(compars == .x$compars[1]),
                                       ggplot2::aes(label = paste0("R = ", round(value, 3))),
                                       x = -Inf, y = Inf,
                                       hjust = -0.1, vjust = 1.1,
                                       size = 2.5, color = "#9A2631"
                   )
      )

  }
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
    if (savePPT) {
      export::graph2ppt(x =.x, file =paste0(savepath,"/",.y, ".ppt"),
                        vector.graphic = TRUE, width = singleplotsize[1], height = singleplotsize[2], aspectr = sqrt(2), append = FALSE)
    }
  })


  .save_zcp <- function(Fig,FigName,outputfile,widths,heights,ppt=FALSE){
    Filepaths=paste0(outputfile,"/",FigName,c(".pdf",".png",".ppt"))
    ggplot2::ggsave(Filepaths[1], width =widths, plot=Fig,height = heights,device = "pdf")
    ggplot2::ggsave(Filepaths[2], width =widths,  plot=Fig,height = heights,device = "png")

    if (ppt) {
      export::graph2ppt(x = Fig, file =Filepaths[3],
                        vector.graphic = TRUE, width =widths, height =heights, aspectr = sqrt(2), append = FALSE)
    }
  }

  .save_zcp(Fig =plt,FigName = "facetwrapPlot_cor",outputfile =savepath,widths = facetplotsize[1],heights = facetplotsize[2],ppt=savePPT)

  result=list(facet_plt=plt,singleplot=plots)
  return(result)

}
