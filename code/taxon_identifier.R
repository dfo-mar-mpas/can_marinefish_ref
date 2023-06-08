# demonstration script of how to get all wikidata identifiers for a taxon via
# gbif lookup

# NOTES This script is not vectorized over species, and could be a lot more
# efficient by reducing the number of queries per taxon from two to one, and by
# only generating the lut once instead of once for every column.

#adapted from --- https://gist.github.com/PietrH/bdbadd18564dedfbb5defb5c89a17160

taxon_identifier <- function(x){

# load libraries ----------------------------------------------------------

library(rgbif)
library(WikidataQueryServiceR)
library(dplyr)
library(purrr)

# lookup gbif id for taxon string -----------------------------------------

# has only one value for it's identifiers
# taxon_name_query <- "Verticillium dahliae" 

# has mutliple iNaturalist taxon id's
taxon_name_query="Polypera greeni"

species_key <- rgbif::name_backbone(taxon_name_query)[[1]]


# query wikidata for all external identifiers -----------------------------

# NOTE species with multiple values for a single identifier will return that
# identifier twice, to prevent this use SELECT DISTINCT, I chose to keep these
# because I want to warn when this occurs

sparql_result <-
  WikidataQueryServiceR::query_wikidata(sprintf(
    'SELECT ?species ?propertyLabel ?property WHERE {
    ?species wdt:P846 "%i";
    ?prop ?object.
    ?property wikibase:claim ?prop;
    rdfs:label ?propertyLabel ;
    wikibase:propertyType wikibase:ExternalId.
    FILTER (lang(?propertyLabel) = "en")
  }',
    species_key
  ))

# we map a number of properties for easier mapping later on
wikidata_object <- pull(sparql_result, species) %>%
  basename() %>%
  unique()

# getting the labels for all returned values (keeping duplicates)
wikidata_properties_all <- pull(sparql_result, property) %>%
  basename()

# removing duplicates
wikidata_properties <- unique(wikidata_properties_all)

# get property values via spaql query -------------------------------------

# this function converts the wikidata property to it's label via lookup table we
# generate from the sparql query result

prop_to_label <- function(prop_string) {
  lut <- sparql_result %>%
    transmute(property = basename(property), propertyLabel)
  return(unique(lut[[2]][lut[[1]] == prop_string]))
}


out_table <- sprintf(
  "SELECT %s WHERE {
  %s
}",
  paste0("?", wikidata_properties, collapse = " "),
  paste0(
    "wd:",
    wikidata_object,
    " wdt:",
    wikidata_properties,
    " ?",
    wikidata_properties,
    ".",
    collapse = "\n"
  )
) %>%
  WikidataQueryServiceR::query_wikidata() %>%
  rename_with(function(x) map_chr(x, prop_to_label)) %>%
  mutate(key = species_key) # add species key for easy disambiguation

if (nrow(out_table) > 1) {
  warning(
    taxon_name_query,
    " has multiple identifiers for the property ",
    prop_to_label(wikidata_properties_all[duplicated(wikidata_properties_all)])
  )
}

# print the output table --------------------------------------------------

return(out_table)

}
