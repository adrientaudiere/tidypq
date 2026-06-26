# Known reagent- and laboratory-contaminant genera

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A curated vector of bacterial genera repeatedly reported as reagent, kit
and laboratory contaminants in low-biomass sequencing studies. The list
is compiled from Salter et al. (2014), the `micRoclean` workflow and the
GRIMER `contaminants.yml` / `human-related.yml` resources. It is
intended as a starting blocklist for
[`identify_contam_blocklist_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_blocklist_pq.md);
users should adapt it to their own substrate and primers.

## Usage

``` r
known_contaminant_genera()
```

## Value

A character vector of genus names.

## References

Salter SJ, Cox MJ, Turek EM, et al. (2014) Reagent and laboratory
contamination can critically impact sequence-based microbiome analyses.
BMC Biology 12:87.
[doi:10.1186/s12915-014-0087-z](https://doi.org/10.1186/s12915-014-0087-z)

## Author

Adrien Taudière

## Examples

``` r
head(known_contaminant_genera())
#> [1] "Afipia"        "Aquabacterium" "Asticcacaulis" "Aurantimonas" 
#> [5] "Beijerinckia"  "Bosea"        
length(known_contaminant_genera())
#> [1] 83
```
