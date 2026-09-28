# Installing Data packages for the data analysis

install.packages(c("readxl", "tidyverse", "pheatmap"))

# Loading the installed packages

library(readxl)
library(tidyverse)
library(pheatmap)

# Chooosing file and calling the data to make sure the right file has ben called
data <- read_excel(file.choose())
head(data)

# Start Data analysis - Renaming the columns to understand easily
data_plot <- data %>%
  select(
    Protein = Description,
    CTRL_Adapted = `Abundances (Grouped): CTRL, Adapted`,
    CTRL_STRESS = `Abundances (Grouped): CTRL, STRESS`,
    PEG10_Adapted = `Abundances (Grouped): 10% PEG, Adapted`,
    PEG10_STRESS = `Abundances (Grouped): 10% PEG, STRESS`,
    PEG20_Adapted = `Abundances (Grouped): 20% PEG, Adapted`,
    PEG20_STRESS = `Abundances (Grouped): 20PEG, STRESS`
  )

# Calculating the protein abundance changes between STRESS and ADAPTED for every protein
data_plot <- data_plot %>%
  mutate(
    CTRL_log2FC = log2(CTRL_Adapted / CTRL_STRESS),
    PEG10_log2FC = log2(PEG10_Adapted / PEG10_STRESS),
    PEG20_log2FC = log2(PEG20_Adapted / PEG20_STRESS)
  )

# Putting top 10 protein to the data - putting all protein to the data will be difficult to read
top10 <- data_plot %>%
  arrange(desc(abs(CTRL_log2FC))) %>%
  slice_head(n = 10)

#Creating the bar graph - configuring the title, x and y axis name, etc.
ggplot(top10, aes(x = reorder(Protein, CTRL_log2FC), y = CTRL_log2FC)) +
  geom_col() +
  coord_flip() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Top 10 Protein Abundance Changes: CTRL",
    x = "Protein",
    y = "log₂ fold change (Adapted / STRESS)"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "← Higher in STRESS",
    hjust = -0.05,
    vjust = 1.5
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Higher in ADAPTED →",
    hjust = 1.05,
    vjust = 1.5
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

