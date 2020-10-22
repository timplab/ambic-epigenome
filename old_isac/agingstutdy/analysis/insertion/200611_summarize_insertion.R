library(tidyverse)

fp <- "/kyber/Data/Nanopore/projects/ambic/sigma/insertion/insertion_point_candidates.txt"

dat <- read_tsv(fp,col_names = c("lab","rep","qname","qstart","qend","chrom","strand","start","end"))

# add cell label
# and get rounded coordinate (to the nearest kb from the center)
dat <- dat %>%
  mutate(cell = ifelse(grepl("Stable",lab),"Stable","Unstable"),
    coord = round((start + end)/2/1e3) * 1e3
  )

# summarize number per cell type
insert.candidates <- dat %>%
  group_by(cell,chrom,coord) %>%
  summarize(n = n()) %>%
  filter (n >= 2) %>%
  spread(cell,n) %>%
  replace(is.na(.),0)

out_path <- "/kyber/Data/Nanopore/projects/ambic/sigma/insertion/insert_candidates_summarized.txt"
write_tsv(insert.candidates,out_path)

