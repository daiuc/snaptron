library(tidyverse)
library(data.table)

DEBUG_MODE = F

if (!DEBUG_MODE) {
  print("### Working in R script mode!")
  args = commandArgs(trailingOnly = T)
  input.num.file = args[1]
  output.file = args[2]

  print(paste0("input: ", input.num.file))
  print(paste0("output: ", output.file, ".gz"))
} else {
  print("### Debug mode!")
  input.num.file = "results/SubsetCounts/GTEx/Bladder/perind_num.counts.gz"
}


# remove clusters that have 0 total reads
print("Remove clusters with total reads <= 5.")

perind_nums = fread(input.num.file, sep = " ", header = F, skip = 1)
header = read_lines(input.num.file, n_max = 1) %>% scan(text = ., what = character())
colnames(perind_nums) = c("chrom", header)

perind_nums[, cluster := str_extract(chrom, "clu_[0-9]+_[-+]")]
df = data.table(cluster = perind_nums$cluster,
                row_sum = select(perind_nums, where(is.numeric)) %>% rowSums
                )

keep_cluster = df[, .(cluster_sum = sum(row_sum)), by = cluster][
  cluster_sum > 5
]$cluster

perind_nums[cluster %in% keep_cluster][, -c("cluster")] %>% fwrite(output.file, sep = " ")
