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
    CTRL_Adapted = `Abundances (Grouped): CTRL, Adapted`,
    CTRL_STRESS = `Abundances (Grouped): CTRL, STRESS`
  )

# Convert the data into a format for plotting
plot_data <- data_plot %>%
  pivot_longer(
    cols = c(CTRL_Adapted, CTRL_STRESS),
    names_to = "Condition",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(Abundance), is.finite(Abundance))

# Select the 10 proteins with the highest abundance
top10 <- plot_data %>%
  group_by(Protein) %>%
  summarise(
    Max_Abundance = max(Abundance, na.rm = TRUE)
  ) %>%
  arrange(desc(Max_Abundance)) %>%
  slice_head(n = 10)

# Keep only those 10 proteins
top10_data <- plot_data %>%
  filter(Protein %in% top10$Protein)

# Horizontal bar graph
ggplot(top10_data, aes(x = Protein, y = Abundance, fill = Condition)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Protein Abundance: STRESS vs Adapted",
    x = "Protein",
    y = "Protein Abundance",
    fill = "Condition"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )
