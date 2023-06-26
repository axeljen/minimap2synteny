# import tidyverse
library(tidyverse)

# define args as input argument vector
args <- commandArgs(trailingOnly = TRUE)

# chromcords as first as syntney map as second argument
chromcords <- read.table(args[1], header = TRUE, sep = "\t")
syntenyMap <- read.table(args[2], header = TRUE, sep = "\t")
labels <- strsplit(args[3], ",")


# now for testing purposes, we will use the example data
#chromcords <- read.table("genomes.chroms.txt", header = TRUE, sep = "\t")
#syntenyMap <- read.table("genomes.syntenyMap.txt", header = TRUE, sep = "\t")

# first thing is to plot the chromosomes
# set the ymin and ymax based on the genome, lowest integer having the highest y-value

# number of genomes, and number of tracks
ngenomes <- length(unique(chromcords$genome))
# a list storing the ymin and ymaxes for each genome
gencords <- list()
for (i in 1:ngenomes){
  gencords[[as.character(i)]] <- c((ngenomes +1) - i + 0.05, (ngenomes +1) - i - 0.05)
  # add these to the chromcords df
  chromcords$ymin[chromcords$genome == i] <- gencords[[as.character(i)]][1]
  chromcords$ymax[chromcords$genome == i] <- gencords[[as.character(i)]][2]
}

# add them also to the syntenymap
for (i in 2:ngenomes){
  syntenyMap$ymin[syntenyMap$query_genome == i] <- gencords[[as.character(i)]][1]
  syntenyMap$ymax[syntenyMap$query_genome == i] <- gencords[[as.character(i - 1)]][2]
}

p <- ggplot() +
  geom_rect(data = chromcords, aes(xmin = cumpos, xmax = cumpos + length, ymin = ymin, ymax = ymax),
            linejoin = "round",color = "black", fill = "lightgray") +
  #scale_fill_manual(values = c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "grey", "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  labs(x = "Chromosome position", y = "Genome", title = "Synteny map") +
  theme(plot.title = element_text(hjust = 0.5))
p

# now add the synteny map as lines between the chromosomes
# remove Nones
syntenyMap <- syntenyMap %>%
  filter(ref_pos != "None") %>%
  # and make integers of all positions
  mutate(query_pos = as.double(query_pos), 
         ref_pos = as.double(ref_pos))

# add them to p as segments
p <- p +
geom_segment(data = syntenyMap, aes(x = query_pos, xend = ref_pos, y = ymin, yend = ymax, color = refgenome_chrom),
             alpha = 0.02,)
p
# add x and y-cordinates to labels
labels <- as.data.frame(labels)
colnames(labels) <- 'label'
for (i in 1:nrow(labels)){
  ycord <- chromcords$ymax[chromcords$genome == i][1] + (chromcords$ymin[chromcords$genome == i][1] - chromcords$ymax[chromcords$genome == i][1]) / 2
  print(ycord)
  labels$ypos[i] <- ycord
}

# add a scale bar in the bottom left corner
p <- p + 
  geom_segment(aes(x = 0, xend = 100000000, y = 0.5, yend = 0.5)) +
  ylim(c(0, ngenomes + 1)) +
  # add whiskers
  geom_segment(aes(x = 0, xend = 0, y = .55,yend = .45)) +
  geom_segment(aes(x = 100000000, xend = 100000000, y = .55,yend = .45)) +
  # remove axes
  theme_void() +
  xlim(c(-100000000,max(chromcords$cumpos) + max(chromcords$length))) +
  theme(legend.position = "none") +
  labs(title = "") +
  # add labels to each chromosome
  geom_text(data=labels, aes(x = -10000000, y = ypos, label = label), hjust = 1) + # and label
  geom_text(aes(x = 50000000, y = .2, label = "100Mb"), size = 3)

#ggsave(sub(".txt$","_plot.svg","genomes.chroms.txt"),width = 12, height = ngenomes)
ggsave(sub(".txt$","_plot.svg",args[1]),width = 12, height = ngenomes)