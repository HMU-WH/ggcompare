`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  } else {
    return(x)
  }
}


#' @export
#'
#' @title Add Brackets with Labels to a ggplot
#'
#' @description
#' Add brackets with labeled annotations to a ggplot.
#'
#' @inheritParams ggplot2::layer
#' @param ... additional arguments passed on to [layer()][ggplot2::layer()].
#' @param arrow `arrow`, the arrow appear at both ends of the brackets, created by [grid::arrow()].
#' @param parse `logical`, whether to parse the labels as expressions and displayed as described in [plotmath][grDevices::plotmath].
#' @param bracket `logical`, whether to display the bracket. If `FALSE`, only the label will be displayed.
#'
#' @returns `LayerInstance`, a layer that can be added to a ggplot.
#'
#' @section Aesthetics:
#' - required: `xmin`, `xmax`, `ymin`, `ymax`, `label`
#' - optional: `angle`, `alpha`, `hjust`, `vjust`, `colour`, `family`, `fontsize`, `fontface`, `linetype`, `linewidth`, `lineheight`
#'
#' @author HMU-WH
#'
#' @examplesIf interactive()
#' ggplot(mpg, aes(class, displ)) +
#'   geom_boxplot() +
#'   annotate("bracket", xmin = 2, xmax = 4, ymin = 4.5, ymax = 5, label = "label",
#'            arrow = grid::arrow(type = "closed", length = unit(0.1, "inches")))
geom_bracket <- function(mapping = NULL, data = NULL,
                         stat = "identity", position = "identity", ...,
                         arrow = NULL, parse = FALSE, bracket = TRUE, inherit.aes = TRUE) {
  ggplot2::layer(stat = stat,
                 data = data,
                 geom = "bracket",
                 mapping = mapping,
                 position = position,
                 show.legend = FALSE,
                 inherit.aes = inherit.aes,
                 params = list(arrow = arrow,
                               parse = parse,
                               bracket = bracket, ...))
}


#' @export
#'
#' @title Add Mean Comparison for Groups to a ggplot
#'
#' @description
#' Add group mean comparisons to a ggplot. The comparisons can be performed using the t-test, Wilcoxon rank-sum test, one-way ANOVA, or Kruskal-Wallis test.
#'
#' @inheritParams ggplot2::layer
#' @param ... additional arguments passed on to [geom_bracket()].
#' @param nudge `numeric`, the nudge of start position in fraction of scale range.
#' @param start `numeric`, the bracket start position. Defaults to the maximum value of `y`.
#' @param breaks `numeric`, the breaks for p-value labels, like c(0, 0.001, 0.01, 0.05, 1).
#' @param labels `character`, the labels for p-value breaks, like c("***", "**", "*", "ns").
#' @param cutoff `numeric`, the cutoff for p-value, labels above this value will be removed.
#' @param method `function`, the method for the test; it should support formula interface and return a list with components `p.value` and `method` (name).
#' @param ref_group `character`, the reference group for comparison. other groups will be compared to this group.
#' @param tip_length `numeric`, the length of the bracket tips in fraction of scale range.
#' @param parametric `logical`, whether to use parametric test (t-test, One-way ANOVA) or non-parametric test (Wilcoxon rank sum test, Kruskal-Wallis test). Applicable only when `method` is NULL.
#' @param correction `character`, the method for p-value adjustment; options include [p.adjust.methods][stats::p.adjust.methods] with "`none`" as the default.
#' @param panel_indep `logical`, whether to correct the p-value only at the panel level. If `FALSE`, the p-value will be corrected at the layer level.
#' @param method_args `list`, additional arguments to be passed to the test method.
#' @param comparisons `list`, a list of comparisons to be made. Each element should contain two groups to be compared.
#' @param step_increase `numeric`, the step increase in fraction of scale range for every additional comparison, in order to avoid overlapping brackets.
#'
#' @returns `LayerInstance`, a layer that can be added to a ggplot.
#'
#' @details
#' Usually you do not need to specify the test method, you only need to tell `stat_compare()` whether you want to perform a parametric test or a nonparametric test, and `stat_compare()` will automatically choose the appropriate test method based on your data.
#' For comparisons between two groups, the p-value is calculated by t-test (parametric) or Wilcoxon rank sum test (nonparametric). For comparisons among more than two groups, the p-value is calculated by One-way ANOVA (parametric) or Kruskal-Wallis test (nonparametric).
#'
#' @section Aesthetics:
#' - required: `x`, `y`
#'
#' @section Computed variables:
#' - **`p`**: p-value of the test.
#' - **`q`**: adjusted p-value of the test.
#' - **`label`**: the label of the bracket.
#' - **`method`**: the method name of the test.
#' - **`xmin`**, **`xmax`**, **`ymin`**, **`ymax`**: the position of the bracket.
#'
#' @author HMU-WH
#'
#' @examplesIf interactive()
#' p <- ggplot(mpg, aes(class, displ, color = class)) +
#'   geom_boxplot(show.legend = FALSE) +
#'   theme_test()
#'
#' # Global comparison: Each x has only one group.
#' p + stat_compare()
#' # If you just want to display text, you can set parameters "bracket" to FALSE.
#' p + stat_compare(bracket = FALSE)
#' # If you want to display the test method, you can do this.
#' p + stat_compare(aes(label = after_stat(sprintf("%s: %s", method, label))))
#'
#' # Comparison between two groups: specify a reference group.
#' p + stat_compare(ref_group = "minivan")
#' # If you only want to display the p-value less or equal to 0.01, you can do this.
#' p + stat_compare(ref_group = "minivan", cutoff = 0.01)
#' # if you want to display the significance level, you can do this.
#' p + stat_compare(ref_group = "minivan", breaks = c(0, 0.001, 0.01, 0.05, 1))
#'
#' # Comparison between two groups: specify the comparison group.
#' p + stat_compare(tip_length = 0.05,
#'                  step_increase = 0,
#'                  comparisons = list(c("compact", "midsize"), c("pickup", "suv")),
#'                  arrow = grid::arrow(type = "closed", length = unit(0.1, "inches")))
#'                  # Yeah, this supports adding arrows.
#'
#' # Within-group (grouped by the x-axis) population comparison.
#' ggplot(mpg, aes(drv, displ, fill = class)) +
#'   geom_boxplot() +
#'   stat_compare() +
#'   stat_compare(aes(group = drv), nudge = 0.1, color = "gray") + # add global comparison
#'   theme_test()
#'
#' # Better adaptation to faceting.
#' ggplot(mpg, aes(drv, displ)) +
#'   geom_boxplot() +
#'   stat_compare(comparisons = combn(unique(mpg$drv), 2, simplify = FALSE)) +
#'   facet_grid(cols = vars(class), scales = "free") +
#'   theme_test()
#'
#' # P-value correction
#' p <- ggplot(mpg, aes(class, displ)) +
#'   geom_boxplot() +
#'   facet_grid(cols = vars(cyl), scales = "free") +
#'   theme_test()
#' # Layer-level P-value correction
#' p + stat_compare(ref_group = 1, correction = "fdr")
#' # Panel-level P-value correction
#' p + stat_compare(ref_group = 1, correction = "fdr", panel_indep = TRUE)
stat_compare <- function(mapping = NULL, data = NULL, position = "identity", ...,
                         nudge = 0, start = NULL, breaks = NULL, labels = NULL, cutoff = NULL, method = NULL, ref_group = NULL, tip_length = 0.02,
                         parametric = FALSE, correction = "none", panel_indep = FALSE, method_args = NULL, comparisons = NULL, step_increase = 0.1, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    geom = "bracket",
    stat = "compare",
    mapping = mapping,
    position = position,
    show.legend = FALSE,
    inherit.aes = inherit.aes,
    params = list(nudge = nudge,
                  start = start,
                  breaks = breaks,
                  labels = labels,
                  cutoff = cutoff,
                  method = method,
                  ref_group = ref_group,
                  tip_length = tip_length,
                  parametric = parametric,
                  correction = correction,
                  panel_indep = panel_indep,
                  method_args = method_args,
                  comparisons = comparisons,
                  step_increase = step_increase, ...))
}


GeomBracket <- ggplot2::ggproto(
  "GeomBracket",
  ggplot2::Geom,
  extra_params = c("na.rm", "orientation"),
  required_aes = c("xmin", "xmax", "ymin", "ymax", "label"),
  default_aes = ggplot2::aes(angle = 0,
                             alpha = NA,
                             hjust = 0.5,
                             vjust = 0.5,
                             family = "",
                             fontface = 1,
                             linetype = 1,
                             fontsize = 3.88,
                             linewidth = 0.5,
                             lineheight = 1.2,
                             colour = "#000000"),
  setup_params = function(data, params) {
    if ("arrow" %in% names(params)) {
      arrow <- params[["arrow"]]
      stopifnot(is.null(arrow) || inherits(arrow, "arrow"))
    }
    if ("parse" %in% names(params)) {
      parse <- params[["parse"]]
      stopifnot(is.logical(parse) && length(parse) == 1 && ! is.na(parse))
    }
    if ("bracket" %in% names(params)) {
      bracket <- params[["bracket"]]
      stopifnot(is.logical(bracket) && length(bracket) == 1 && ! is.na(bracket))
    }
    params[["flipped"]] <- ggplot2::has_flipped_aes(data, params)
    return(params)
  },
  draw_key = function(data, params, size) {
    grid::nullGrob()
  },
  draw_group = function(data, panel_params, coord, parse = FALSE, arrow = NULL, bracket = TRUE, flipped = FALSE) {
    data <- coord[["transform"]](data, panel_params)
    if (flipped == inherits(coord, "CoordFlip")) {
      x <- mean(c(data[["xmin"]], data[["xmax"]]))
      y <- ifelse(bracket, data[["ymax"]], data[["ymin"]])
      x0 <- c(data[["xmin"]], data[["xmax"]], data[["xmin"]])
      y0 <- c(data[["ymax"]], data[["ymax"]], data[["ymax"]])
      x1 <- c(data[["xmin"]], data[["xmax"]], data[["xmax"]])
      y1 <- c(data[["ymin"]], data[["ymin"]], data[["ymax"]])
    } else {
      data[["angle"]] <- data[["angle"]] - 90
      y <- mean(c(data[["ymin"]], data[["ymax"]]))
      x <- ifelse(bracket, data[["xmax"]], data[["xmin"]])
      x0 <- c(data[["xmax"]], data[["xmax"]], data[["xmax"]])
      y0 <- c(data[["ymin"]], data[["ymax"]], data[["ymin"]])
      x1 <- c(data[["xmin"]], data[["xmin"]], data[["xmax"]])
      y1 <- c(data[["ymin"]], data[["ymax"]], data[["ymax"]])
    }
    if (parse) {
      data[["label"]] <- utils::getFromNamespace("parse_safe", "ggplot2")(data[["label"]])
    }
    text_grob <- grid::textGrob(
      x = x, y = y,
      rot = data[["angle"]], label = data[["label"]],
      hjust = data[["hjust"]], vjust = data[["vjust"]] - 1,
      gp = grid::gpar(
        fontfamily = data[["family"]],
        fontface = data[["fontface"]],
        lineheight = data[["lineheight"]],
        fontsize = data[["fontsize"]] * ggplot2::.pt,
        col = ggplot2::alpha(data[["colour"]], data[["alpha"]]))
    )
    segment_grob <- (
      if (bracket) {
        grid::gList(
          grid::segmentsGrob(
            x0 = x0[[3]],
            y0 = y0[[3]],
            x1 = x1[[3]],
            y1 = y1[[3]],
            default.units = "native",
            gp = grid::gpar(lty = data[["linetype"]],
                            lwd = data[["linewidth"]] * ggplot2::.pt,
                            col = ggplot2::alpha(data[["colour"]], data[["alpha"]]))
          ),
          grid::segmentsGrob(
            x0 = x0[1:2],
            y0 = y0[1:2],
            x1 = x1[1:2],
            y1 = y1[1:2],
            arrow = arrow,
            default.units = "native",
            gp = grid::gpar(lty = data[["linetype"]],
                            lwd = data[["linewidth"]] * ggplot2::.pt,
                            col = ggplot2::alpha(data[["colour"]], data[["alpha"]]))
          )
        )
      } else {
        grid::nullGrob()
      }
    )
    return(grid::gList(text_grob, segment_grob))
  }
)


StatCompare <- ggplot2::ggproto(
  "StatCompare",
  ggplot2::Stat,
  required_aes = c("x", "y", "group"),
  extra_params = c("na.rm", "cutoff", "breaks", "labels", "panel_indep", "orientation"),
  setup_params = function(self, data, params) {
    if ("nudge" %in% names(params)) {
      nudge <- params[["nudge"]]
      stopifnot(is.numeric(nudge) && length(nudge) == 1 && ! is.na(nudge))
    }
    if ("start" %in% names(params)) {
      start <- params[["start"]]
      stopifnot(is.null(start) || (is.numeric(start) && length(start) == 1 && ! is.na(start)))
    }
    if ("breaks" %in% names(params)) {
      breaks <- params[["breaks"]]
      stopifnot(is.null(breaks) || (is.numeric(breaks) && min(breaks) <= 0 && max(breaks) >= 1))
    }
    if ("labels" %in% names(params)) {
      labels <- params[["labels"]]
      stopifnot(is.null(labels) || (is.character(labels) && length(labels) == length(breaks) - 1))
    }
    if ("cutoff" %in% names(params)) {
      cutoff <- params[["cutoff"]]
      stopifnot(is.null(cutoff) || (is.numeric(cutoff) && length(cutoff) == 1 && ! is.na(cutoff)))
    }
    if ("method" %in% names(params)) {
      method <- params[["method"]]
      if (! is.null(method)) { params[["method"]] <- match.fun(method) }
    }
    if ("ref_group" %in% names(params)) {
      ref_group <- params[["ref_group"]]
      stopifnot(is.null(ref_group) || ((is.numeric(ref_group) || is.character(ref_group)) && length(ref_group) == 1 && ! is.na(ref_group)))
    }
    if ("tip_length" %in% names(params)) {
      tip_length <- params[["tip_length"]]
      stopifnot(is.numeric(tip_length) && length(tip_length) == 1 && ! is.na(tip_length))
    }
    if ("parametric" %in% names(params)) {
      parametric <- params[["parametric"]]
      stopifnot(is.logical(parametric) && length(parametric) == 1 && ! is.na(parametric))
    }
    if ("correction" %in% names(params)) {
      correction <- params[["correction"]]
      params[["correction"]] <- match.arg(correction, stats::p.adjust.methods)
    }
    if ("panel_indep" %in% names(params)) {
      panel_indep <- params[["panel_indep"]]
      stopifnot(is.logical(panel_indep) && length(panel_indep) == 1 && ! is.na(panel_indep))
    }
    if ("method_args" %in% names(params)) {
      method_args <- params[["method_args"]]
      stopifnot(is.null(method_args) || is.list(method_args))
    }
    if ("comparisons" %in% names(params)) {
      comparisons <- params[["comparisons"]] <- unique(params[["comparisons"]])
      stopifnot(is.null(comparisons) || (is.list(comparisons) && length(comparisons) >= 1))
      if (! all(vapply(comparisons, \(x) (is.numeric(x) || is.character(x)) && length(x) == 2 , logical(1)))) {
        stop("'comparisons' must be a list of character or numeric vectors of length 2 ...")
      }
    }
    if ("step_increase" %in% names(params)) {
      step_increase <- params[["step_increase"]]
      stopifnot(is.numeric(step_increase) && length(step_increase) == 1 && ! is.na(step_increase))
    }
    if (inherits(data[["x"]], "mapped_discrete") == inherits(data[["y"]], "mapped_discrete")) {
      stop("Can only handle data with groups plotted on either the x-axis or y-axis, but not both ...")
    }
    params[["flipped"]] <- ggplot2::has_flipped_aes(data, params)
    if (is.null(params[["ref_group"]] %||% params[["comparisons"]])) {
      data <- ggplot2::flip_data(data, params[["flipped"]])
      counts <- vapply(split(data, ~ x + PANEL), \(x) { length(unique(x[["group"]])) }, numeric(1))
      if (any(counts > 1)) {
        params[["global"]] <- FALSE
        params[["multiple"]] <- max(counts) > 2
      } else {
        params[["global"]] <- TRUE
        params[["multiple"]] <- max(vapply(split(data[["group"]], data[["PANEL"]]), \(x) { length(unique(x)) }, numeric(1))) > 2
      }
    }
    return(params)
  },
  compute_layer = function(self, data, params, layout) {
    data <- ggplot2::flip_data(data, params[["flipped"]])
    constant_aes <- split(data[, setdiff(colnames(data), c("y", "group")), drop = FALSE], ~ x + PANEL) |>
      lapply(\(x) { as.data.frame(lapply(x[, setdiff(colnames(x), c("x", "PANEL")), drop = FALSE], \(y) { length(unique(stats::na.omit(y))) })) }) |>
      (\(x) { do.call(rbind, args = x) })()
    constant_aes <- unique(data[, c("x", "PANEL", colnames(constant_aes)[vapply(constant_aes, \(x) { all(x == 1) }, logical(1))]), drop = FALSE])
    data <- ggproto_parent(Stat, self)$compute_layer(data, params, layout)
    constant_aes <- constant_aes[, union(c("x", "PANEL"), setdiff(colnames(constant_aes), colnames(data))), drop = FALSE]
    if ("x" %in% colnames(data)) {
      data <- merge(data, constant_aes, by = c("x", "PANEL"), all.x = TRUE)
      data <- data[, setdiff(colnames(data), "x"), drop = FALSE]
    }
    data <- transform(data, flipped_aes = params[["flipped"]])
    data <- ggplot2::flip_data(data, flip = params[["flipped"]])
    if (! params[["panel_indep"]]) {
      data[["q"]] <- stats::p.adjust(data[["p"]], method = params[["correction"]])
    }
    breaks <- params[["breaks"]]
    if (is.null(breaks)) {
      data[["label"]] <- ifelse(data[["q"]] < .Machine[["double.eps"]], sprintf("p < %.2e", .Machine[["double.eps"]]), sprintf("%.2g", data[["q"]]))
    } else {
      labels <- params[["labels"]]
      if (is.null(labels)) {
        labels <- c(vapply(rev(seq_len(length(breaks) - 2)), \(x) { paste(rep("*", x), collapse = "") }, character(1)), "ns")
      }
      data[["label"]] <- as.character(cut(data[["q"]], breaks = breaks, labels = labels, include.lowest = TRUE))
    }
    cutoff <- params[["cutoff"]]
    if (! is.null(cutoff)) {
      data[["label"]] <- ifelse(data[["q"]] > cutoff, NA, data[["label"]])
    }
    if (length(which(is.na(data[["label"]]) & ! is.na(data[["p"]]))) > 0) {
      data <- split(data, data[["PANEL"]]) |>
        lapply(\(x) {
          x <- x[order(x[["ymin"]], decreasing = TRUE), , drop = FALSE]
          for (i in which(is.na(x[["label"]]) & ! is.na(x[["p"]]))) {
            x[seq_len(i), "ymin"] <- x[seq_len(i), "ymin"] - x[seq_len(i), "space"]
            x[seq_len(i), "ymax"] <- x[seq_len(i), "ymax"] - x[seq_len(i), "space"]
          }
          return(x[order(x[["group"]]), , drop = FALSE])
        }) |>
        (\(x) { do.call(rbind, args = x) })()
    }
    data <- data[, setdiff(colnames(data), "space"), drop = FALSE]
    return(data)
  },
  compute_panel = function(data, scales, nudge = 0, start = NULL, global = FALSE, method = NULL, flipped = FALSE, multiple = FALSE, ref_group = NULL, correction = "none", parametric = FALSE, method_args = NULL, comparisons = NULL, step_increase = 0.1, tip_length = 0.03) {
    scales <- ggplot2::flip_data(scales, flipped)
    scale_range <- diff(scales[["y"]][["range"]][["range"]])
    bracket_spacing <- ifelse(step_increase == 0, 0, scale_range * step_increase)
    if (is.null(start)) { start <- scales[["y"]][["range"]][["range"]][2] + nudge * scale_range }
    if (is.null(ref_group) && is.null(comparisons)) {
      .compare <- function(data) {
        return(
          tryCatch({
            do.call(method %||% ifelse(multiple,
                                       ifelse(parametric, stats::oneway.test, stats::kruskal.test),
                                       ifelse(parametric, stats::t.test, stats::wilcox.test)),
                    args = c(list(formula = y ~ group, data = data), method_args)) |>
              (\(x) { data.frame(p = x[["p.value"]], q = NA, method = ifelse(is.na(x[["p.value"]]), NA, x[["method"]])) })()
          }, error = \(e) { warning(e[["message"]]); data.frame(p = NA, q = NA, method = NA) })
        )
      }
      data <- (
        if (global) {
          data.frame(.compare(data), xmin = min(data[["x"]]) - 0.45, xmax = max(data[["x"]]) + 0.45, ymin = start, ymax = start + tip_length * scale_range, space = 0, group = 0)
        } else {
          lapply(split(data, data[["x"]]), \(x) { data.frame(.compare(x), x = x[["x"]][[1]], xmin = x[["x"]][[1]] - 0.45, xmax = x[["x"]][[1]] + 0.45, ymin = start, ymax = start + tip_length * scale_range, space = 0, group = x[["x"]][[1]]) }) |>
            (\(x) { do.call(rbind, args = x) })()
        }
      )
    } else {
      groups <- sort(unique(data[["x"]]))
      if (length(groups) <= 1) {
        comparisons <- NULL
      } else {
        if (is.null(comparisons)) {
          ref_group <- scales[["x"]][["map"]](ref_group)
          if (ref_group %in% groups) {
            comparisons <- lapply(sort(setdiff(groups, ref_group)), \(x) { c(x, ref_group) })
          }
        } else {
          comparisons <- lapply(comparisons, \(x) { scales[["x"]][["map"]](x) })
        }
      }
      i <- 0
      data <- lapply(comparisons, \(comp) {
        compare <- tryCatch({
          do.call(method %||% ifelse(parametric, stats::t.test, stats::wilcox.test), c(list(formula = y ~ x, data = data[data[["x"]] %in% comp, c("x", "y")]), method_args)) |>
            (\(x) { data.frame(p = x[["p.value"]], q = NA, method = ifelse(is.na(x[["p.value"]]), NA, x[["method"]])) })()
        }, error = \(e) { warning(e[["message"]]); data.frame(p = NA, q = NA, method = NA) })
        bracket_start <- start + i * bracket_spacing
        annotation_start <- bracket_start + tip_length * scale_range
        if (! all(is.na(compare))) { i <<- i + 1 }
        return(
          data.frame(
            compare,
            x = comp[[1]],
            xmin = min(comp),
            xmax = max(comp),
            ymin = bracket_start,
            ymax = annotation_start,
            space = bracket_spacing,
            group = paste(sort(comp), collapse = "-")
          )
        )
      }) |>
        (\(x) { do.call(rbind, args = x) })()
    }
    data[["q"]] <- p.adjust(data[["p"]], method = correction)
    return(data)
  }
)
