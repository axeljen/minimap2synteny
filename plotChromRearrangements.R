library(tidyverse)
library(ggrepel)
library(scales)


#### Read in the files from command args
args <- commandArgs(trailingOnly = TRUE)

# fetch genomes file
#genomes <- read_tsv("/Users/axeljensen/Dropbox/syntenymaps_testing/test_genomes.txt")
genomes <- read_tsv(args[1],
                      col_names = T)
# syntenymap
#synmap <- read_tsv("/Users/axeljensen/Dropbox/syntenymaps_testing/test_syntenymap.txt")
synmap <- read_tsv(args[2])
# rearrangements
#rearrangements <- read_tsv("/Users/axeljensen/Dropbox/syntenymaps_testing/out_rearrangements.txt")
rearrangements <- read_tsv(args[3])

# split the ref and query genome into different df
refgen <- genomes[genomes$type == "ref",]
qgen <- genomes[genomes$type == "query",]

# gap to add between chroms reference chroms
gap <- 50000000

refgen <- refgen %>%
  mutate(cumstart = cumsum(length + gap) - length - gap)

# calculate the gaplength in the query genome to match that in the refgenome
ngapsref <- nrow(refgen) - 1
totgaplen <- ngapsref * gap
ngapsq <- nrow(qgen) - 1
gapq <- as.integer(totgaplen / ngapsq)

# and add cumpos to qgenome too
qgen <- qgen %>%
  mutate(cumstart = cumsum(length + gapq) - length - gapq)

# bind them together and add ytack
gens <- rbind(refgen, qgen) %>%
  mutate(ytrack = if_else(type == 'query',1,2))

# add cumulative position to the syntenymap
synmap <- synmap %>%
  # first get rid of nones
  filter(qstart != "None") %>%
  # and make sure that qpos is numeric
  mutate(qstart = as.numeric(qstart)) %>%
  mutate(refcum = tpos + gens$cumstart[match(tchrom,gens$chrom)]) %>%
  mutate(qcum = qstart + gens$cumstart[match(qchrom,gens$chrom)]) %>%
  # add col with ycords 
  mutate(ymin = 1.2, ymax = 1.8) 
  #mutate(cumstart = 0, length = 0, ytrack = 0)

# plot syntenymap as a base
p <- ggplot() +
geom_curve(synmap, mapping = aes(x = refcum, xend = refcum,
                                       y = 2, yend = ymax, color = tchrom),
                 alpha = .01,
           curvature = -.1) +
  geom_curve(synmap, mapping = aes(x = qcum, xend = qcum,
                                     y = 1, yend = ymin, color = tchrom),
               alpha = .01,
               curvature = .1) +
  geom_curve(synmap, mapping = aes(x = refcum, xend = qcum,
                                    y = ymax, yend = ymin,color = tchrom),
                 alpha = .01,
            curvature = .1) +
  theme(legend.position = "none") +
  geom_segment(data = gens, aes(x = cumstart, xend = cumstart + length,
                                                                 y = ytrack, yend = ytrack),
                                                linewidth = 3, lineend = "round") +
  geom_segment(data = gens, aes(x = cumstart, xend = cumstart + length,
                                y = ytrack, yend = ytrack),
               linewidth = 1.5, lineend = "round", color = "white") +
  theme_void() +
  theme(legend.position = "none")
p
# add cumpos to rearrangments
rearr <- rearrangements %>%
  mutate(refstart_1 = as.numeric(refstart_1),
         refstart_2 = as.numeric(refstart_2),
         refend_1 = as.numeric(refend_1),
         refend_2 = as.numeric(refend_2),
         querystart_1 = as.numeric(querystart_1),
         querystart_2 = as.numeric(querystart_2),
         queryend_1 = as.numeric(queryend_1),
          queryend_2 = as.numeric(queryend_2)) %>%
  mutate(cumstart_1 = refstart_1 + gens$cumstart[match(refchrom_1,gens$chrom)],
         cumend_1 = refend_1 + gens$cumstart[match(refchrom_1,gens$chrom)],
          cumstart_2 = refstart_2 + gens$cumstart[match(refchrom_2,gens$chrom)],
          cumend_2 = refend_2 + gens$cumstart[match(refchrom_2,gens$chrom)]) %>%
  mutate(x_1 = cumstart_1 + ((cumend_1 - cumstart_1) / 2),
         x_2 = cumstart_2 + ((cumend_2 - cumstart_2) / 2))

# plot fissions as vertical lines on top of ref genome
p <- p + geom_segment(data = rearr[rearr$event == "fission",],
                 aes(y = 2.05, yend = 1.95,x = x_1, xend = x_1,
                 color = 'red')) +
  ylim(c(0,3)) +
  # add labels specifying the breakpoints
  geom_label_repel(data = rearr[rearr$event == "fission",],
                   aes(x = x_1,
                       y = 2.05, 
                       label = paste0(comma(refstart_1, scale = .000001,
                                            accuracy = 0.01), "-", 
                                      comma(refend_1, scale = .000001,
                                            accuracy = 0.01), "Mb")),
                  direction = "y",
                  nudge_y = 1,
                  min.segment.length = 0,
                  force_pull = 0
                  )

p

# plot translocations
p <- p + geom_segment(data = rearr[rearr$event == "translocation",],
              aes(x = x_1, xend = x_2, y = 1.99, yend = 1.99),
              linewidth = 1.5, color = "darkblue") +
  # and ynversions
  geom_segment(data = rearr[rearr$event == "inversion",],
               aes(x = x_1, xend = x_2, y = 2.01, yend = 2.01),
               linewidth = 1.5, color = "orange")
p <- p +
  geom_curve(data = rearr[rearr$event == "fusion",],
            aes(x = x_1, xend= x_2, y = 1.99, yend = 1.99),
            color = "red", curvature = 0.2)
# save the plot
ggsave(sub(".txt$","_plot.svg",args[3]),width = 12, height = 5)

