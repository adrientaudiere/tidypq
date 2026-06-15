# Print a wrapped pheatmap object

Draws the `gtable` stored in a `wrapped_pheatmap` object (as returned in
the `heatmap_by_category` slot of
[`neg_control_classify_pq()`](https://adrientaudiere.github.io/tidypq/reference/neg_control_classify_pq.md)).

## Usage

``` r
# S3 method for class 'wrapped_pheatmap'
print(x, ...)
```

## Arguments

- x:

  A `wrapped_pheatmap` object.

- ...:

  Ignored, present for S3 compatibility.

## Value

`x`, invisibly. Called for the side effect of drawing the heatmap.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Author

Adrien Taudière
