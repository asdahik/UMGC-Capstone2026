# Install packages once if needed
install.packages(c("readxl", "tidyverse"))

# Load packages
library(readxl)
library(tidyverse)

# Import your Excel file
data <- read_excel(file.choose())

# Select the protein abundance columns
data_plot <- data %>%
  select(
    Protein = Description,
    PEG10_Adapted = `Abundances (Grouped): 10% PEG, Adapted`,
    PEG10_STRESS = `Abundances (Grouped): 10% PEG, STRESS`
  )

# Convert the data into a format for plotting
plot_data <- data %>%
  select(
    Protein = Description,
    Adapted = `Abundances: 10% PEG, Adapted`,
    STRESS = `Abundances: 10% PEG, STRESS`
  ) %>%
  pivot_longer(
    cols = c(Adapted, STRESS),
    names_to = "Condition",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(Abundance), is.finite(Abundance))

# Select the 10 proteins with the highest abundance
ctrl_data <- data %>%
  select(
    Protein = Description,
    Adapted = `Abundances : CTRL, Adapted`,
    STRESS = `Abundances : CTRL, STRESS`
  )

top10 <- ctrl_data %>%
  mutate(Max_Abundance = pmax(Adapted, STRESS, na.rm = TRUE)) %>%
  arrange(desc(Max_Abundance)) %>%
  distinct(Protein, .keep_all = TRUE) %>%
  slice_head(n = 10)

# Keep only those 10 proteins
top10_data <- plot_data %>%
  filter(Protein %in% top10$Protein) %>%
  mutate(
    Protein = factor(
      Protein,
      levels = rev(as.character(top10$Protein))
    )
  )

# Horizontal bar graph
ggplot(top10_data, aes(x = Protein, y = Abundance, fill = Condition)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Protein Abundance: 10% PEG",
    x = "Protein",
    y = "Protein Abundance",
    fill = "Condition"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

