#Rcode

library(tidyverse)

# set working directory
setwd("/home/user06/ngs-course-final-assignment-25-26-Dimple-bit-14")
dir.create("results")
# read the file
d <- read_tsv("cols-all.tsv",
              col_types = cols(
                #CHROM = col_character(),
                POS = col_integer(),
                ID = col_character(),
                REF = col_character(),
                ALT = col_character(),
                QUAL = col_double(),
                DP = col_double(),
                SNP.INDEL = col_character(),
                Mutation_type = col_character()
              ))

# check data
glimpse(d)
head(d)

# Barplot: SNP vs INDEL
p1 <- d %>%
  ggplot(aes(x = SNP.INDEL)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "SNP vs INDEL count",
       x = "Variant Type",
       y = "Count")
ggsave("results/SNP_vs_INDEL.png", p1, width = 6, height = 4)

# Barplot: Transition vs Transversion
p2 <- d %>%
  filter(Mutation_type != "NA") %>%
  ggplot(aes(x = Mutation_type)) +
  geom_bar(fill = "darkgreen") +
  theme_minimal() +
  labs(title = "Transition vs Transversion",
       x = "Mutation type",
       y = "Count")
ggsave("results/Transition_vs_Transversion.png", p2, width = 6, height = 4)

# Histogram: DP distribution
p3 <- d %>%
  ggplot(aes(x = DP)) +
  geom_histogram(bins = 50, fill = "purple") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Depth distribution",
       x = "DP (log scale)",
       y = "Frequency")
ggsave("results/DP_distribution.png", p3, width = 6, height = 4)

# Boxplot: DP by variant type
p4 <- d %>%
  ggplot(aes(x = SNP.INDEL, y = DP)) +
  geom_boxplot(fill = "orange") +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "DP by variant type",
       x = "Variant type",
       y = "DP (log scale)")
ggsave("results/DP_by_variant_type.png", p4, width = 6, height = 4)

# Histogram: QUAL distribution per variant type
p5 <- d %>%
  filter(QUAL < 900) %>%
  ggplot(aes(x = QUAL)) +
  geom_histogram(bins = 50, fill = "red") +
  scale_x_log10() +
  facet_wrap(~ SNP.INDEL) +
  theme_minimal() +
  labs(title = "QUAL distribution",
       x = "QUAL",
       y = "Count")
ggsave("results/QUAL_distribution.png", p5, width = 8, height = 4)

#Barplot: Variants per chromosome
p6 <- d %>%
  count(`#CHROM`) %>%
  ggplot(aes(x = `#CHROM`, y = n)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Variants per chromosome",
       x = "Chromosome",
       y = "Count")
ggsave("results/Variants_per_chromosome.png", p6, width = 8, height = 4)

#Calculate Ti/Tv ratio
ti <- d %>% filter(Mutation_type == "Transition") %>% nrow()
tv <- d %>% filter(Mutation_type == "Transversion") %>% nrow()
ti_tv_ratio <- ti / tv

# Save Ti/Tv ratio as a text file
writeLines(paste0("Transition/Transversion ratio: ", round(ti_tv_ratio, 3)),
           "results/Ti_Tv_ratio.txt")

# Print to console
ti_tv_ratio
