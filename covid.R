# Load data
var <- read.csv("variants_long_table.csv", sep = ";")

#Load packages
library(tidyverse)
library("RColorBrewer")

# Transform the data frame into a tibble
var_tb <- as_tibble(var)
var_tb

# Progressing in variants exploration

# We are going to plot:
# The distribution of DP values per chromosome and per sample
# The variants effect per sample
# The genes showing variants effects
#The read depth per position

# The distribution of DP values per chromosome and per sample
ggplot(data = var_tb, aes(x = CHROM, y = DP, fill = SAMPLE)) +
  geom_boxplot() +
  ylim(0, 10000) +
  scale_fill_brewer(palette = "RdYlBu") +
  facet_grid(~SAMPLE ) + 
  theme(legend.position = "bottom") +
  labs(title = "DP per Chromosome", x = "Chromosome", y = "DP") 


# Adding options
ggplot(data = var_tb, aes(x = CHROM, y = DP, fill = SAMPLE)) +
  ylim(0, 10000) +
  scale_fill_brewer(palette = "RdYlBu") +
  labs(title = "DP per Chromosome", x = "Chromosome", y = "DP") +
  geom_violin(trim=FALSE) +
  facet_grid(~SAMPLE ) +
  geom_boxplot(width=0.1) +
  theme_bw()

# Variants effects per sample
# Plotting the variants effects
# Count number of different effects per sample
ggplot(data = var_tb, aes(x = EFFECT, fill = SAMPLE)) +
  scale_fill_brewer(palette = "RdBu") +
  ggtitle("Effect per Sample") +
  theme(legend.position = "bottom") + 
  geom_bar() 

# Flipping the axis
ggplot(data = var_tb, aes(x = EFFECT, fill = SAMPLE)) +
  scale_fill_brewer(palette = "RdBu") +
  ggtitle("Effect per Sample") +
  theme(legend.position = "bottom") + 
  geom_bar() +
  coord_flip() 


# Counting and extracting all effects for all genes
# Counting the effects per gene
var_tb %>% count(EFFECT, GENE, sort = TRUE)

# Counting and extracting variants specific effects for all genes
# counting the variants modifying the stop codons

#Filtering option 1 to select for effect on stop
filter(var_tb, EFFECT == "stop_lost" | EFFECT == "stop_gained")

# Option 2
filter(var_tb, EFFECT %in% c("stop_lost", "stop_gained"))

# Filtering effect on selected column
filter(var_tb, EFFECT %in% c("stop_lost", "stop_gained")) %>% select(SAMPLE, CHROM, GENE, EFFECT)

# DP per position
# Define your variable

p_DP_POS <- ggplot(data = var_tb, aes(x=POS, y=DP, fill= SAMPLE)) + scale_fill_brewer(palette="RdBu") + labs(title="DP_per_Position") + theme(legend.position="bottom")

# Plot

p_DP_POS + geom_point(shape = 21, size = 5)

# Plot with transparency options

p_DP_POS + geom_point(shape = 21, size = 5, alpha = 0.7)

