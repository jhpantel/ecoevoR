#' Plot gen3sis results, compare across multiple simulation runs
#'
#' @param output a sgen3sis output object resulting from a gen3sis simulation (i.e. run_simulation; see gen3sis::run_simulation)
#' @param summary_title summary plot title as character. If NULL, title is computed from input name (see gen3sis::run_simulation)
#' @param summary_legend either a staring with _\_n for new lines or NULL. If NULL, provides default summary and simulation information (see gen3sis::run_simulation)
#' @param max_gam A maximum gamma diversity value across multiple simulation runs. Sets the right-y-axis for the gamma diversity plot.
#' @param max_alpha A maximum alpha diversity value across multiple simulation runs. Sets the color scale for the alpha diversity map plot.
#'
#' @return no values returned, only for plotting
#' @export
#'
#' @examples
#' # load existing summary example
#' datapath <- system.file(file.path("extdata", "WorldCenter"), package = "gen3sis")
#' output <- readRDS(file.path(datapath, "output/config_worldcenter/sgen3sis.rds"))
#' # plot output summary
#' max_gam <- max(output$summary$phylo_summary[,"alive"])
#' max_alpha <- max(output$summary$`richness-final`[,3],na.rm=T)
#' plot_gen3sis(output,max_gam=max_gam, max_alpha=max_alpha)
plot_gen3sis <- function (output, summary_title = NULL, summary_legend = NULL, max_gam=500, max_alpha=100)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (class(output) != "gen3sis_output") {
    stop("this is not  a gen3sis_output object")
  }
  {
    layout(matrix(c(1, 3, 3, 1, 3, 3, 1, 3, 3, 2, 2, 2, 2, 2, 2), ncol = 3, byrow = TRUE))
    par(mar = c(4, 3, 3, 7), oma = c(0.1, 0.8, 0.3, 0.8))
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(0,10), ylim = c(0, 10), axes = FALSE, ann = FALSE)
    if (is.null(summary_title)) {
      summary_title <- str_split(output$parameters$directories$input,"/")[[1]]
      summary_title <- paste(tail(summary_title, 2)[1],
                             collapse = ">")
    }
    if (is.null(summary_legend)) {
      sum_names <- names(output$parameters$gen3sis$general)
      sumss <- output$parameters$gen3sis$general
      sumar <- output$summary
      phylo <- sumar$phylo_summary[c(1, nrow(sumar$phylo_summary)),
      ]
      col_sum <- colSums(sumar$phylo_summary)
      summary_legend = paste(paste(sum_names[2], sumss[2],sep = ": "), paste("end_time;", tail(names(sumar$occupancy),1)), paste("traits", paste0(sumss[7][[1]], collapse = ","),sep = ": "), paste("world_habited_present", paste0(round(sumar$occupancy[length(sumar$occupancy)],1) * 100, "%"), sep = ": "), paste("initial_richness", phylo[1, "total"], sep = ": "), paste("cumulative_richness", phylo[2, "total"], sep = ": "), paste("richness_present", phylo[2, "alive"], sep = ": "), paste("speciation", paste0(round((col_sum["speciations"] - phylo["initial", "total"])/phylo[1, "total"] * 100, 0), "%"), sep = ": "), paste("extinction", paste0(round((col_sum["extinctions"])/phylo[1, "total"] * 100, 0), "%"), sep = ": "), paste(names(output$system)[1], round(output$system$`runtime-hours`, 2), sep = ": "), sep = "\n")
    }
    par(xpd = TRUE)
    # legend("topleft", inset = c(0, -0.4), title = paste0("Summary [",
    #                                                      summary_title, "]"), legend = summary_legend, bty = "n",
    #        title.adj = 0)
    d <- output$summary$phylo_summary[-1, -1]
    plot(d[, "alive"], xlab = "", ylab = "", type = "l", col = "black", lwd = 4, frame.plot = FALSE, xaxt = "n", yaxt = "n",ylim=c(0,round_up(max_gam)))
    axis(4, line = -1, cex = 1, cex.axis = 1, col = "black")
    mtext(side = 4, text = "γ richness", col = "black",
          line = 2, cex = 1.1)
    par(new = TRUE)
    plot(d[, "speciations"], pch = 3, col = rgb(0, 0, 1, 0.5), xlab = "", ylab = "", type = "b", frame.plot = FALSE, xaxt = "n", yaxt = "n", ylim = c(0,round_up(max_gam)))
    points(d[, "extinctions"], pch = 4, col = rgb(1, 0, 0, 0.5), type = "b")
    axis(2, line = -1, cex = 1, cex.axis = 1, col = "black")
    mtext(side = 2, text = "Evolutionary events", col = "black", line = 1.5, cex = 1.1)
    legend(x = 1, y = 80, legend = c("Richness", "Speciation", "Extinction"), col = c("black", rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), pch = c(NA, 3, 4), lty = c(1, NA, NA), lwd = c(4, NA, NA), bty = "n")
    axis_lab <- seq(from = 1, to = nrow(d), length.out = max((nrow(d)/20), 2))
    axis(1, at = axis_lab, labels = rownames(d)[axis_lab])
    mtext(side = 1, text = "Time steps", line = 2.5, cex = 1.1)
    ras <- rasterFromXYZ(output$summary$`richness-final`)
    max_ras <- max_alpha # Max value from my 4 simulations
    #max_ras <- max(ras@data@values, na.rm = TRUE)
    min_ras <- min(ras@data@values, na.rm = TRUE)
    zerorichness_col <- "navajowhite3"
    if (max_ras == 0) {
      rc <- zerorichness_col
    } else {
      rc <- color_richness(max_ras)
    }
    if (min_ras == 0) {
      rc <- c(zerorichness_col, rc)
    }
    #image(ras, col = rc, bty = "o", xlab = "", ylab = "",
    #      las = 1, asp = 1, breaks=breaks,col=cols)
    breaks <- seq(0, max_alpha, by = 0.01)
    cols <- colorRampPalette(rc)(length(breaks) - 1)
    image(rasterFromXYZ(output$summary$`richness-final`),breaks=breaks, col = cols)
    mtext(4, text = "Final α richness", line = 1, cex = 1.2)
    # Plot scale and colors
    rng <- seq(0,max_alpha,10)
    arg <- list(at=rng, labels=round(rng, 4))

    raster::plot(rasterFromXYZ(output$summary$`richness-final`),breaks=breaks,legend.only = TRUE, add = TRUE, col = cols, axis.args=arg)
  }
}
