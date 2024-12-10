#' @title corpairs: Correlation Pairs Plotting Function
#' @description This function provides correlation matrix visualizations using `GGally::ggpairs`.
#' It supports plotting scatter plots with density-based color mappings or hexbin plots.
#' @param dt A data frame to be used for plotting. The first column should be an identifier column.
#' @param id.col The column name of the identifier column in `dt`. Default is "Index".
#' @param cor.method The correlation method to use: "pearson", "kendall", or "spearman". Default is "spearman".
#' @param savefile The directory where output plots will be saved. Default is "outputfile".
#' @param corpairplotsize A numeric vector of length 2 specifying the width and height of the plot in inches. Default is `c(6, 5)`.
#' @param bin Number of bins for hexbin plots. Default is `50`.
#' @param logL Logical, whether to log-transform the data. Default is `FALSE`.
#' @param plottype The type of plot for lower panels: "hex" or "point". Default is "point".
#' @param pointcolor A vector of color palette names for scatter plots, with a default of "C".
#' @param pointsize Size of points in scatter plots. Default is `1`.
#' @param kde2d.n Grid size for `MASS::kde2d` density calculation. Default is `50`.
#' @return A `GGally::ggpairs` plot object.
#' @examples
#' # Load the demo data included in the package
#' data("demo", package = "corhex")
#' corpairs(dt=demo,
#' id.col="Index",
#' cor.method=c("pearson", "kendall", "spearman")[3],
#' savefile="outputfile",
#' corpairplotsize=c(3*2,2.5*2),#width height
#' bin=50,
#' logL=FALSE,
#' plottype=c("hex","point")[2],
#' pointcolor=c("A","B","C","D","E")[3],
#' pointsize=1,
#' kde2d.n=50)
#' @importFrom dplyr select any_of everything
#' @importFrom ggplot2 ggplot geom_hex scale_fill_viridis_c theme_bw theme aes geom_point scale_color_viridis_c annotate theme_void theme_minimal labs ggsave
#' @importFrom GGally ggpairs wrap
#' @importFrom MASS kde2d
#' @importFrom pracma interp2
#' @importFrom tibble column_to_rownames
#' @importFrom stats cor complete.cases
#' @importFrom scales col_numeric
#' @importFrom export graph2ppt
#' @export
corpairs=function(dt=demo,
                  id.col="Index",
                  cor.method=c("pearson", "kendall", "spearman")[3],
                  savefile="outputfile",
                  corpairplotsize=c(3*3,2.5*3),#width height
                  bin=50,
                  logL=FALSE,
                  plottype=c("hex","point")[2],
                  pointcolor=c("A","B","C","D","E")[3],
                  pointsize=1,
                  kde2d.n=50){


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

  #hex function
  geomhex <- function(data, mapping,binss, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_hex(bins = binss) +
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
      )

  }


 #geompoint-related function
  compute_density <- function(x, y) {
    dens <- MASS::kde2d(x, y, n = kde2d.n)
    interp_dens <- pracma::interp2(dens$x, dens$y, dens$z, x, y)
    return(interp_dens)
  }

  geompoint_density <- function(data, mapping, ...) {
    x_var_name <- rlang::as_label(mapping$x)
    y_var_name <- rlang::as_label(mapping$y)

    x_var <- data[[x_var_name]]
    y_var <- data[[y_var_name]]

    valid_indices <- stats::complete.cases(x_var, y_var)
    x_var <- x_var[valid_indices]
    y_var <- y_var[valid_indices]

    density <- compute_density(x_var, y_var)

    data <- data[valid_indices, , drop = FALSE]  # 保留有效行
    data$density <- density

    ggplot2::ggplot(data, mapping) +
      ggplot2::geom_point(ggplot2::aes(color = density), alpha = 0.8, size = pointsize) +  # 动态大小
      ggplot2::scale_color_viridis_c(option = pointcolor, guide = "none") +               # 动态调色板
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
        axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(t = 0.5)),
        axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = ggplot2::margin(r = 0.5)),
        axis.text = ggplot2::element_text(size = 8, color = "black"),
        axis.line.x = ggplot2::element_line(color = "black"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = "italic"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "right",
        strip.text = ggplot2::element_text(size = 8, face = 'italic', color = "grey20")
      ) +
      ggplot2::labs(
        title = NULL,
        x = NULL,
        y = NULL
      )
  }


  #diag function
  custom_densityDiag <- function(data, mapping, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_density(alpha = 0.8, fill = "#26828E", ...) +
      ggplot2::theme_minimal() +
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
          size = 8, face = 'italic', color = "grey20",
          hjust = 0.5, halign = 0.5, r = grid::unit(5, "pt"),
          width = grid::unit(5.5, "npc"),
          padding = ggplot2::margin(3, 0, 3, 0),
          margin = ggplot2::margin(1, 1, 1, 1)
        ),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold'),
        panel.spacing = grid::unit(1, 'lines')
      )
  }


  cor_data <- as.data.frame(as.table(corvalue))
  cor_data <- cor_data[as.numeric(cor_data$Var1) < as.numeric(cor_data$Var2), ]
  cor_data$Var1 <- as.character(cor_data$Var1)
  cor_data$Var2 <- as.character(cor_data$Var2)

  color_mapping <- function(value) {
    if (value > 0) {
      scales::col_numeric(palette = c("black", "#CD4071"), domain = c(0, 1))(value)
    } else {
      scales::col_numeric(palette = c("black", "#31688E"), domain = c(-1, 0))(value)
    }
  }
  custom_upper_plots <- lapply(seq_len(nrow(cor_data)), function(idx) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = sprintf("%.2f", cor_data$Freq[idx]),
                        #size = abs(cor_data$Freq[idx]) *4,
                        size = max(min(abs(cor_data$Freq[idx]) * 4, 10), 3),
                        color =color_mapping(cor_data$Freq[idx])) +
      ggplot2::theme_void() +
      ggplot2::theme(aspect.ratio = 1)
  })



  #plotting data

  dt.ggpairs=dt |> tibble::column_to_rownames(var = id.col)


  if (!logL){
    dt.ggpairs=log10(dt.ggpairs+1)
  }else{
    dt.ggpairs
  }


  switch (plottype,
    "hex" = {

      ggpairs_plot <- GGally::ggpairs(
        dt.ggpairs,
        upper = list(continuous = function(data, mapping, ...) {
          x_var <- as.character(rlang::get_expr(mapping$x))
          y_var <- as.character(rlang::get_expr(mapping$y))
          # print(paste("x_var:", x_var, "y_var:", y_var))
          idx <- which((cor_data$Var1 == x_var & cor_data$Var2 == y_var) |
                         (cor_data$Var1 == y_var & cor_data$Var2 == x_var))
          # print(paste("idx:", idx))
          if (length(idx) == 1) {
            custom_upper_plots[[idx]]
          } else {
            ggplot2::ggplot() + ggplot2::theme_void()
          }
        }),
        lower = list(continuous = GGally::wrap(geomhex, binss = bin)),
        diag = list(continuous = custom_densityDiag)
      )

    },
    "point"={

      ggpairs_plot <- GGally::ggpairs(
        dt.ggpairs,
        upper = list(continuous = function(data, mapping, ...) {
          x_var <- as.character(rlang::get_expr(mapping$x))
          y_var <- as.character(rlang::get_expr(mapping$y))
          # print(paste("x_var:", x_var, "y_var:", y_var))
          idx <- which((cor_data$Var1 == x_var & cor_data$Var2 == y_var) |
                         (cor_data$Var1 == y_var & cor_data$Var2 == x_var))
          # print(paste("idx:", idx))
          if (length(idx) == 1) {
            custom_upper_plots[[idx]]
          } else {
            ggplot2::ggplot() + ggplot2::theme_void()
          }
        }),
        lower = list(continuous = GGally::wrap(geompoint_density)),
        diag = list(continuous = custom_densityDiag)
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

  .save_zcp <- function(Fig,FigName,outputfile,widths,heights){
    Filepaths=paste0(outputfile,"/",FigName,c(".pdf",".png",".ppt"))
    ggplot2::ggsave(Filepaths[1], width =widths, plot=Fig,height = heights,device = "pdf")
    ggplot2::ggsave(Filepaths[2], width =widths,  plot=Fig,height = heights,device = "png")
    export::graph2ppt(x = Fig, file =Filepaths[3],
                      vector.graphic = TRUE, width =widths, height =heights, aspectr = sqrt(2), append = FALSE)
  }

  .save_zcp(Fig =ggpairs_plot,FigName = "pairsCorplot",outputfile =savepath,widths = corpairplotsize[1],heights = corpairplotsize[2])

  return(ggpairs_plot)

}
