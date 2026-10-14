
# Install packages once if needed
install.packages(c("readxl", "tidyverse"))

# Load packages
library(readxl)
library(tidyverse)

# Import your Excel file
data <- read_excel(file.choose())

# Select protein abundance under 20% PEG drought
drought_data <- data %>%
  select(
    Protein = Description,
    Adapted = `Abundances (Grouped): 20% PEG, Adapted`,
    STRESS = `Abundances (Grouped): 20PEG, STRESS`
  )

# Find the 10 most abundant proteins under 20% PEG
top10 <- drought_data %>%
  mutate(
    Max_Abundance = pmax(Adapted, STRESS, na.rm = TRUE)
  ) %>%
  filter(is.finite(Max_Abundance)) %>%
  arrange(desc(Max_Abundance)) %>%
  distinct(Protein, .keep_all = TRUE) %>%
  slice_head(n = 10)

# Convert selected proteins into plotting format
top10_data <- top10 %>%
  select(Protein, Adapted, STRESS) %>%
  pivot_longer(
    cols = c(Adapted, STRESS),
    names_to = "Condition",
    values_to = "Abundance"
  ) %>%
  filter(is.finite(Abundance)) %>%
  mutate(
    Protein = factor(
      Protein,
      levels = rev(as.character(top10$Protein))
    )
  )

# Create horizontal bar graph
ggplot(top10_data,
       aes(x = Protein, y = Abundance, fill = Condition)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Top 10 Protein Abundances: 20% PEG Drought",
    x = "Protein",
    y = "Protein Abundance",
    fill = "Condition"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )
