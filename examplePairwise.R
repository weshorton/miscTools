### Example plot
library(data.table)
library(ggplot2)

set.seed(1)

## Make data
data_dt <- data.table("Gene" = paste0("Gene", 1:100), "A" = sample(300:1000, 100), "B" = sample(100:800, 100))
meta_dt <- data.table("Gene" = paste0("Gene", 1:100), "Class" = c(rep("good", 50), rep("bad", 50)))

## Format data
melt_dt <- melt(data_dt, id.vars = "Gene")
plot_dt <- merge(melt_dt, meta_dt, by = "Gene", sort = F)

## Format colors
plot_dt$Class <- factor(plot_dt$Class, levels = c("good", "bad"))
colors_v <- c("good" = "green", "bad" = "red")

ggplot(data = plot_dt, aes(x = variable, y = value)) +
  geom_line(aes(group = Gene, color = Class)) +
  geom_point(aes(group = Gene, color = Class)) +
  scale_color_manual(values = colors_v, labels = names(colors_v))