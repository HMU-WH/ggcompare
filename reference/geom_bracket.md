# Add Brackets with Labels to a ggplot

Add brackets with labeled annotations to a ggplot.

## Usage

``` r
geom_bracket(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  ...,
  arrow = NULL,
  parse = FALSE,
  bracket = TRUE,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html). If
  specified and `inherit.aes = TRUE` (the default), it is combined with
  the default mapping at the top level of the plot. You must supply
  `mapping` if there is no plot mapping.

- data:

  The data to be displayed in this layer. There are three options:

  If `NULL`, the default, the data is inherited from the plot data as
  specified in the call to
  [`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

  A `data.frame`, or other object, will override the plot data. All
  objects will be fortified to produce a data frame. See
  [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
  for which variables will be created.

  A `function` will be called with a single argument, the plot data. The
  return value must be a `data.frame`, and will be used as the layer
  data. A `function` can be created from a `formula` (e.g.
  `~ head(.x, 10)`).

- stat:

  The statistical transformation to use on the data for this layer. When
  using a `geom_*()` function to construct a layer, the `stat` argument
  can be used the override the default coupling between geoms and stats.
  The `stat` argument accepts the following:

  - A `Stat` ggproto subclass, for example `StatCount`.

  - A string naming the stat. To give the stat as a string, strip the
    function name of the `stat_` prefix. For example, to use
    [`stat_count()`](https://ggplot2.tidyverse.org/reference/geom_bar.html),
    give the stat as `"count"`.

  - For more information and other ways to specify the stat, see the
    [layer
    stat](https://ggplot2.tidyverse.org/reference/layer_stats.html)
    documentation.

- position:

  A position adjustment to use on the data for this layer. This can be
  used in various ways, including to prevent overplotting and improving
  the display. The `position` argument accepts the following:

  - The result of calling a position function, such as
    [`position_jitter()`](https://ggplot2.tidyverse.org/reference/position_jitter.html).
    This method allows for passing extra arguments to the position.

  - A string naming the position adjustment. To give the position as a
    string, strip the function name of the `position_` prefix. For
    example, to use
    [`position_jitter()`](https://ggplot2.tidyverse.org/reference/position_jitter.html),
    give the position as `"jitter"`.

  - For more information and other ways to specify the position, see the
    [layer
    position](https://ggplot2.tidyverse.org/reference/layer_positions.html)
    documentation.

- ...:

  additional arguments passed on to
  [layer()](https://ggplot2.tidyverse.org/reference/layer.html).

- arrow:

  `arrow`, the arrow appear at both ends of the brackets, created by
  [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html).

- parse:

  `logical`, whether to parse the labels as expressions and displayed as
  described in [plotmath](https://rdrr.io/r/grDevices/plotmath.html).

- bracket:

  `logical`, whether to display the bracket. If `FALSE`, only the label
  will be displayed.

- inherit.aes:

  If `FALSE`, overrides the default aesthetics, rather than combining
  with them. This is most useful for helper functions that define both
  data and aesthetics and shouldn't inherit behaviour from the default
  plot specification, e.g.
  [`borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).

## Value

`LayerInstance`, a layer that can be added to a ggplot.

## Aesthetics

- required: `xmin`, `xmax`, `ymin`, `ymax`, `label`

- optional: `angle`, `alpha`, `hjust`, `vjust`, `colour`, `family`,
  `fontsize`, `fontface`, `linetype`, `linewidth`, `lineheight`

## Author

HMU-WH

## Examples

``` r
if (FALSE) { # interactive()
library(ggplot2)

ggplot(mpg, aes(class, displ)) +
  geom_boxplot() +
  annotate("bracket", xmin = 2, xmax = 4, ymin = 4.5, ymax = 5, label = "label",
           arrow = grid::arrow(type = "closed", length = unit(0.1, "inches")))
}
```
