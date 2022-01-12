library(ggplot2)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: plot_heatmaps.R MARKERS OUTPATH.n", call.=FALSE)
}
atlas.path <- args[1]
out.path <- args[2]

load.atlas <- function(atlas.path) {
  df <- as.data.frame(fread(atlas.path))
  #df <- df[order('target'),]
  name <- paste0(df$chr, ':', df$start, '-', df$end)
  dd <- df[, 9:ncol(df)]
  dd$name <- name
  dd$position <- seq(1, nrow(df))
  melted_df <- reshape2::melt(dd, id.vars = c('position', 'name'))
  melted_df
}

data <- load.atlas(atlas.path)

ggplot(data=data, aes(y=position, x=variable, fill=value)) + 
  geom_tile() +
  theme_bw(base_size=18) +
  labs(x="Cell types", y="Markers") + 
  scale_y_reverse(expand=c(0, 0)) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=18, hjust=1),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18, angle=90),
        plot.margin = margin(20, 20, 20, 20),
  ) +
  scale_fill_gradient2(low="yellow", high="blue", mid="white", 
                       midpoint=0.5,  limit=c(0,1), 
                       name="Average\nmethylation") 

ggsave(out.path, width = 16.1, height = 11.2)

