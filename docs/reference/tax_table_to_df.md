# Convert a `tax_table` into a data.frame with numeric and list columns

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The `tax_table` slot of a phyloseq object is stored as a **character**
matrix, so numeric information (e.g. a bootstrap support, a trait value)
is silently coerced to text and list-valued cells cannot be represented
at all. This helper lifts those limitations by returning a tibble in
which columns that look numeric (or logical) are converted to their
natural type, and selected columns can be split into genuine list
columns.

It does **not** modify the phyloseq object; it produces a stand-alone
tibble for downstream (tidyverse) manipulation.

## Usage

``` r
tax_table_to_df(
  physeq,
  convert = TRUE,
  split = NULL,
  sep = "/",
  taxa_names_col = "taxon"
)
```

## Arguments

- physeq:

  (phyloseq or taxonomyTable, required) A phyloseq object with a
  `tax_table` slot, or a `tax_table` object directly.

- convert:

  (logical, default TRUE) If TRUE, each column is passed through
  [`utils::type.convert()`](https://rdrr.io/r/utils/type.convert.html)
  (with `as.is = TRUE`) so that columns holding only numbers become
  numeric/integer and `"TRUE"`/`"FALSE"` columns become logical. Columns
  with any non-numeric value are kept as character.

- split:

  (character, optional) Names of columns to turn into list columns by
  splitting each cell on `sep`. Useful when a rank cell packs several
  values (e.g. multiple GBIF matches).

- sep:

  (character, default `"/"`) Separator used to split the columns named
  in `split`. Passed to
  [`base::strsplit()`](https://rdrr.io/r/base/strsplit.html) with
  `fixed = TRUE`.

- taxa_names_col:

  (character or NULL, default `"taxon"`) Name of the leading column
  holding the taxa names. If NULL, taxa names are kept as row names
  instead of a column.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with one row per taxon and one column per taxonomic rank, plus (unless
`taxa_names_col` is NULL) a leading column of taxa names.

## See also

[`phyloseq::tax_table()`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
library(MiscMetabar)
# Numeric-looking ranks are converted to numeric
df <- tax_table_to_df(data_fungi_mini)
df
#> # A tibble: 45 × 13
#>    taxon Domain Phylum Class Order Family Genus Species Trophic.Mode Guild Trait
#>    <chr> <chr>  <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>        <chr> <chr>
#>  1 ASV7  Fungi  Basid… Agar… Russ… Stere… NA    NA      Saprotroph   Wood… NULL 
#>  2 ASV8  Fungi  Basid… Agar… Russ… Stere… Ster… ostrea  Saprotroph   Unde… Whit…
#>  3 ASV12 Fungi  Basid… Agar… Hyme… Schiz… Xylo… radulo… Saprotroph   Unde… Whit…
#>  4 ASV18 Fungi  Basid… Agar… Russ… Stere… Ster… ostrea  Saprotroph   Unde… Whit…
#>  5 ASV25 Fungi  Basid… Agar… Agar… Lyoph… Ossi… lachno… Saprotroph   Wood… Brow…
#>  6 ASV26 Fungi  Basid… Agar… Russ… Stere… Ster… hirsut… Saprotroph   Unde… Whit…
#>  7 ASV27 Fungi  Basid… Agar… Poly… Stecc… Antr… brasil… Saprotroph   Wood… Whit…
#>  8 ASV29 Fungi  Basid… Agar… Auri… Exidi… Basi… eyrei   Saprotroph   Unde… NULL 
#>  9 ASV32 Fungi  Basid… Agar… Cant… Canth… Sist… oblong… Saprotroph-… Ecto… Whit…
#> 10 ASV34 Fungi  Basid… Agar… Agar… Entol… Ento… NA      Saprotroph   Unde… NULL 
#> # ℹ 35 more rows
#> # ℹ 2 more variables: Confidence.Ranking <chr>, Genus_species <chr>

# Keep taxa names as row names instead of a column
tax_table_to_df(data_fungi_mini, taxa_names_col = NULL)
#>        Domain        Phylum              Class           Order
#> ASV7    Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV8    Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV12   Fungi Basidiomycota     Agaricomycetes Hymenochaetales
#> ASV18   Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV25   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV26   Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV27   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV29   Fungi Basidiomycota     Agaricomycetes  Auriculariales
#> ASV32   Fungi Basidiomycota     Agaricomycetes  Cantharellales
#> ASV34   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV35   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV41   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV42   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV46   Fungi Basidiomycota Atractiellomycetes  Atractiellales
#> ASV47   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV48   Fungi Basidiomycota     Agaricomycetes  Auriculariales
#> ASV49   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV50   Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV53   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV54   Fungi Basidiomycota     Agaricomycetes  Auriculariales
#> ASV58   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV59   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV61   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV62   Fungi Basidiomycota     Agaricomycetes            <NA>
#> ASV63   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV64   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV67   Fungi Basidiomycota     Agaricomycetes Hymenochaetales
#> ASV68   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV71   Fungi Basidiomycota    Tremellomycetes            <NA>
#> ASV72   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV75   Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV77   Fungi Basidiomycota     Agaricomycetes      Agaricales
#> ASV82   Fungi Basidiomycota     Agaricomycetes  Auriculariales
#> ASV83   Fungi Basidiomycota    Tremellomycetes     Tremellales
#> ASV85   Fungi Basidiomycota     Agaricomycetes Hymenochaetales
#> ASV91   Fungi Basidiomycota     Agaricomycetes  Auriculariales
#> ASV93   Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV94   Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV99   Fungi Basidiomycota     Agaricomycetes     Polyporales
#> ASV100  Fungi Basidiomycota               <NA>            <NA>
#> ASV101  Fungi Basidiomycota     Agaricomycetes     Corticiales
#> ASV104  Fungi Basidiomycota     Agaricomycetes      Russulales
#> ASV105  Fungi Basidiomycota     Agaricomycetes Hymenochaetales
#> ASV107  Fungi Basidiomycota     Agaricomycetes Hymenochaetales
#> ASV108  Fungi Basidiomycota     Agaricomycetes  Auriculariales
#>                                    Family             Genus       Species
#> ASV7                           Stereaceae              <NA>          <NA>
#> ASV8                           Stereaceae           Stereum        ostrea
#> ASV12                      Schizoporaceae           Xylodon    raduloides
#> ASV18                          Stereaceae           Stereum        ostrea
#> ASV25                       Lyophyllaceae        Ossicaulis     lachnopus
#> ASV26                          Stereaceae           Stereum      hirsutum
#> ASV27                     Steccherinaceae       Antrodiella  brasiliensis
#> ASV29                          Exidiaceae    Basidiodendron         eyrei
#> ASV32   Cantharellales_fam_Incertae_sedis        Sistotrema oblongisporum
#> ASV34                      Entolomataceae          Entocybe          <NA>
#> ASV35                        Polyporaceae             Fomes   fomentarius
#> ASV41                    Tricholomataceae            Mycena        renati
#> ASV42                       Lyophyllaceae        Ossicaulis     lachnopus
#> ASV46   Atractiellales_fam_Incertae_sedis       Helicogloea     pellucida
#> ASV47                         Pterulaceae       Radulomyces       molaris
#> ASV48                         Aporpiaceae          Elmerina        caryae
#> ASV49                   Phanerochaetaceae     Phanerochaete     livescens
#> ASV50       Russulales_fam_Incertae_sedis Gloeohypochnicium      analogum
#> ASV53                        Polyporaceae             Fomes   fomentarius
#> ASV54                     Auriculariaceae       Auricularia          <NA>
#> ASV58                        Polyporaceae             Fomes   fomentarius
#> ASV59                    Hyphodermataceae        Hyphoderma  roseocremeum
#> ASV61                    Hyphodermataceae        Hyphoderma     setigerum
#> ASV62                                <NA>              <NA>          <NA>
#> ASV63                                <NA>              <NA>          <NA>
#> ASV64                        Polyporaceae          Trametes    versicolor
#> ASV67                      Schizoporaceae           Xylodon    raduloides
#> ASV68                       Lyophyllaceae        Ossicaulis     lachnopus
#> ASV71                                <NA>              <NA>          <NA>
#> ASV72                                <NA>              <NA>          <NA>
#> ASV75                      Peniophoraceae        Peniophora   versiformis
#> ASV77                       Lyophyllaceae        Ossicaulis     lachnopus
#> ASV82                          Exidiaceae            Exidia    glandulosa
#> ASV83                                <NA>              <NA>          <NA>
#> ASV85  Hymenochaetales_fam_Incertae_sedis     Peniophorella        pubera
#> ASV91                     Auriculariaceae       Auricularia   mesenterica
#> ASV93                          Stereaceae              <NA>          <NA>
#> ASV94                          Stereaceae           Stereum        ostrea
#> ASV99                        Polyporaceae             Fomes   fomentarius
#> ASV100                               <NA>              <NA>          <NA>
#> ASV101                       Corticiaceae   Marchandiomyces        buckii
#> ASV104                        Hericiaceae          Hericium   coralloides
#> ASV105                     Schizoporaceae           Xylodon    flaviporus
#> ASV107                     Schizoporaceae           Xylodon    raduloides
#> ASV108                         Exidiaceae            Exidia    glandulosa
#>                             Trophic.Mode
#> ASV7                          Saprotroph
#> ASV8                          Saprotroph
#> ASV12                         Saprotroph
#> ASV18                         Saprotroph
#> ASV25                         Saprotroph
#> ASV26                         Saprotroph
#> ASV27                         Saprotroph
#> ASV29                         Saprotroph
#> ASV32             Saprotroph-Symbiotroph
#> ASV34                         Saprotroph
#> ASV35                         Saprotroph
#> ASV41              Pathotroph-Saprotroph
#> ASV42                         Saprotroph
#> ASV46                         Saprotroph
#> ASV47                         Saprotroph
#> ASV48                         Saprotroph
#> ASV49                         Saprotroph
#> ASV50                         Saprotroph
#> ASV53                         Saprotroph
#> ASV54                         Saprotroph
#> ASV58                         Saprotroph
#> ASV59                         Saprotroph
#> ASV61                         Saprotroph
#> ASV62                                  -
#> ASV63                                  -
#> ASV64                         Saprotroph
#> ASV67                         Saprotroph
#> ASV68                         Saprotroph
#> ASV71                                  -
#> ASV72                                  -
#> ASV75              Pathotroph-Saprotroph
#> ASV77                         Saprotroph
#> ASV82             Saprotroph-Symbiotroph
#> ASV83  Pathotroph-Saprotroph-Symbiotroph
#> ASV85                         Saprotroph
#> ASV91                         Saprotroph
#> ASV93                         Saprotroph
#> ASV94                         Saprotroph
#> ASV99                         Saprotroph
#> ASV100                                 -
#> ASV101                        Pathotroph
#> ASV104                        Saprotroph
#> ASV105                        Saprotroph
#> ASV107                        Saprotroph
#> ASV108            Saprotroph-Symbiotroph
#>                                                                      Guild
#> ASV7                                  Wood Saprotroph-Undefined Saprotroph
#> ASV8                                                  Undefined Saprotroph
#> ASV12                                                 Undefined Saprotroph
#> ASV18                                                 Undefined Saprotroph
#> ASV25                                                      Wood Saprotroph
#> ASV26                                                 Undefined Saprotroph
#> ASV27                                                      Wood Saprotroph
#> ASV29                                                 Undefined Saprotroph
#> ASV32                                      Ectomycorrhizal-Wood Saprotroph
#> ASV34                                                 Undefined Saprotroph
#> ASV35                                                      Wood Saprotroph
#> ASV41  Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph
#> ASV42                                                      Wood Saprotroph
#> ASV46                                                 Undefined Saprotroph
#> ASV47                                                 Undefined Saprotroph
#> ASV48                                                 Undefined Saprotroph
#> ASV49                                                      Wood Saprotroph
#> ASV50                                                 Undefined Saprotroph
#> ASV53                                                      Wood Saprotroph
#> ASV54                                                 Undefined Saprotroph
#> ASV58                                                      Wood Saprotroph
#> ASV59                                                 Undefined Saprotroph
#> ASV61                                                 Undefined Saprotroph
#> ASV62                                                                    -
#> ASV63                                                                    -
#> ASV64                                                      Wood Saprotroph
#> ASV67                                                 Undefined Saprotroph
#> ASV68                                                      Wood Saprotroph
#> ASV71                                                                    -
#> ASV72                                                                    -
#> ASV75                                       Plant Pathogen-Wood Saprotroph
#> ASV77                                                      Wood Saprotroph
#> ASV82                                       Endophyte-Undefined Saprotroph
#> ASV83                                 Fungal Parasite-Undefined Saprotroph
#> ASV85                                                 Undefined Saprotroph
#> ASV91                                                 Undefined Saprotroph
#> ASV93                                 Wood Saprotroph-Undefined Saprotroph
#> ASV94                                                 Undefined Saprotroph
#> ASV99                                                      Wood Saprotroph
#> ASV100                                                                   -
#> ASV101                                                     Lichen Parasite
#> ASV104                                                Undefined Saprotroph
#> ASV105                                                Undefined Saprotroph
#> ASV107                                                Undefined Saprotroph
#> ASV108                                      Endophyte-Undefined Saprotroph
#>                      Trait Confidence.Ranking              Genus_species
#> ASV7                  NULL           Probable                      NA_NA
#> ASV8             White Rot           Probable             Stereum_ostrea
#> ASV12            White Rot           Probable         Xylodon_raduloides
#> ASV18            White Rot           Probable             Stereum_ostrea
#> ASV25            Brown Rot           Probable       Ossicaulis_lachnopus
#> ASV26            White Rot           Probable           Stereum_hirsutum
#> ASV27            White Rot    Highly Probable   Antrodiella_brasiliensis
#> ASV29                 NULL           Probable       Basidiodendron_eyrei
#> ASV32            White Rot           Possible   Sistotrema_oblongisporum
#> ASV34                 NULL           Probable                Entocybe_NA
#> ASV35  Brown Rot-White Rot           Probable          Fomes_fomentarius
#> ASV41                 NULL           Probable              Mycena_renati
#> ASV42            Brown Rot           Probable       Ossicaulis_lachnopus
#> ASV46                 NULL           Probable      Helicogloea_pellucida
#> ASV47            White Rot           Probable        Radulomyces_molaris
#> ASV48                 NULL           Probable            Elmerina_caryae
#> ASV49            White Rot    Highly Probable    Phanerochaete_livescens
#> ASV50            White Rot           Probable Gloeohypochnicium_analogum
#> ASV53  Brown Rot-White Rot           Probable          Fomes_fomentarius
#> ASV54                 NULL           Probable             Auricularia_NA
#> ASV58  Brown Rot-White Rot           Probable          Fomes_fomentarius
#> ASV59            White Rot           Probable    Hyphoderma_roseocremeum
#> ASV61            White Rot           Probable       Hyphoderma_setigerum
#> ASV62                    -                  -                      NA_NA
#> ASV63                    -                  -                      NA_NA
#> ASV64            White Rot    Highly Probable        Trametes_versicolor
#> ASV67            White Rot           Probable         Xylodon_raduloides
#> ASV68            Brown Rot           Probable       Ossicaulis_lachnopus
#> ASV71                    -                  -                      NA_NA
#> ASV72                    -                  -                      NA_NA
#> ASV75            White Rot           Probable     Peniophora_versiformis
#> ASV77            Brown Rot           Probable       Ossicaulis_lachnopus
#> ASV82                 NULL           Probable          Exidia_glandulosa
#> ASV83                 NULL           Possible                      NA_NA
#> ASV85            White Rot           Probable       Peniophorella_pubera
#> ASV91                 NULL           Probable    Auricularia_mesenterica
#> ASV93                 NULL           Probable                      NA_NA
#> ASV94            White Rot           Probable             Stereum_ostrea
#> ASV99  Brown Rot-White Rot           Probable          Fomes_fomentarius
#> ASV100                   -                  -                      NA_NA
#> ASV101                NULL           Probable     Marchandiomyces_buckii
#> ASV104           White Rot           Probable       Hericium_coralloides
#> ASV105           White Rot           Probable         Xylodon_flaviporus
#> ASV107           White Rot           Probable         Xylodon_raduloides
#> ASV108                NULL           Probable          Exidia_glandulosa
# }

if (FALSE) { # \dontrun{
# Split a packed column into a list column
tax_table_to_df(data_fungi_mini, split = "Genus", sep = "_")
} # }
```
