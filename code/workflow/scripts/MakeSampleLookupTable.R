#' @title Create sample ID look up table
#' @author Chao Dai
#' @description Create a lookup table for sample ID, get data source, and tissue type info.



library(tidyverse)
library(data.table)


DEBUG_MODE = F

if (!DEBUG_MODE) {
  print("### Working in snakemake script mode!")

  gtex.file = snakemake@input[['gtex']]
  tcga.file = snakemake@input[['tcga']]
  srav.file = snakemake@input[['srav']]
  all_perind.file = snakemake@input[['all_perind']]
  output.file = snakemake@output[[1]]
} else {
  print("### Debug mode!")

  gtex.file = "resources/snaptron-yil/scripts/data/samples_gtex.tsv"
  tcga.file = "resources/snaptron-yil/scripts/data/samples_tcga.tsv"
  srav.file = "resources/snaptron-yil/scripts/data/samples_srav2.tsv"
  all_perind.file = "resources/All_perind_RailIDs.txt"
}


gtex = fread(gtex.file)
tcga = fread(tcga.file)
srav = fread(srav.file)
all_perind = fread(all_perind.file, header = F, col.names = c("rail_id", "counts_col_id"))

# gtex columns: rail_id, Histological_Type, ScientificName (organism)
# tcga columns: rail_id, gdc_cases.project.primary_site
# srav columns: rail_id, sample_attribute, taxon_id

gtex.out = select(gtex, rail_id, Histological_Type) %>%
  unique %>%
  filter(!is.na(Histological_Type)) %>%
  add_column(ds = "GTEx") %>%
  rename("tissue" = "Histological_Type")

tcga.out = select(tcga, rail_id, gdc_cases.project.primary_site) %>%
  unique %>%
  filter(!is.na(gdc_cases.project.primary_site)) %>%
  add_column(ds = "TCGA") %>%
  rename("tissue" = "gdc_cases.project.primary_site")

srav.out = select(srav, rail_id, sample_attribute) %>%
  unique %>%
  filter(!is.na(sample_attribute)) %>%
  add_column(ds = "srav") %>%
  rename("tissue" = "sample_attribute")


all_perind_gtex = inner_join(all_perind, gtex.out, by = "rail_id")
all_perind_tcga = inner_join(all_perind, tcga.out, by = "rail_id")
all_perind_srav = inner_join(all_perind, srav.out, by = "rail_id")

all_perind_dict = do.call(rbind, list(all_perind_gtex, all_perind_tcga, all_perind_srav))

setdiff(all_perind$rail_id, all_perind_dict$rail_id) %>% length %>%
  paste0("!!! ", ., " rail_ids dropped due to the lack of tissue type information.")

print(paste0("### Save ", nrow(all_perind_dict), " records with tissue type info."))
write_tsv(all_perind_dict, output.file)
print("Done!")
