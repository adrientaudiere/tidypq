# Diagnostic plot for a `contam_tbl`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A lightweight overview of which detectors flagged how many taxa. For
publication-grade contamination figures, richer plotting lives in the
`ggplotpq` package and consumes the same `contam_tbl`.

## Usage

``` r
# S3 method for class 'contam_tbl'
autoplot(object, ...)

# S3 method for class 'contam_tbl'
plot(x, ...)
```

## Arguments

- object, x:

  (contam_tbl, required) A detector output.

- ...:

  Currently ignored.

## Value

A `ggplot` object (`autoplot`) or `x` invisibly (`plot`).

## Author

Adrien Taudière
