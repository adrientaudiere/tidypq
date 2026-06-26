# Combine several `contam_tbl` objects

Row-binds detector outputs into a single `contam_tbl`. Evidence columns
missing from one detector are filled with `NA`.

## Usage

``` r
# S3 method for class 'contam_tbl'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  (contam_tbl) One or more `contam_tbl` objects.

- deparse.level:

  Ignored; present for [`rbind()`](https://rdrr.io/r/base/cbind.html)
  compatibility.

## Value

A `contam_tbl` holding the union of the flagged taxa.

## Author

Adrien Taudière
