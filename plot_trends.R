###Load required modules
library(ggplot2)
library(ggthemes)
library(reshape2)
library(data.table)

###Open main dataset
dataset <- read_delim("Blood_vs_Spleen_merged_datasets.csv", ",", escape_double = F, trim_ws = T)
dataset <- data.frame(lapply(dataset, function(v) {if (is.character(v)) return(toupper(v))else return(v)}))

###Define function to add multiple plots in a single graph

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###Create a vector with gene of interest
genes <- c('IL10', 'IL21', 'SLAMF8', 'IFNG')

###Extract genes from main dataframe
genesplot <- dataset[dataset$genes %in% genes, ]

###Melt dataframe
toplot <- melt(genesplot)

###Split data by tissue
blood <- toplot[toplot$variable %like% "BD", ]
spleen <- toplot[toplot$variable %like% "SD", ]

###Create groups per timepoint
spleen$variable <- gsub(".*SD0_.*", "SD0", spleen$variable)
spleen$variable <- gsub(".*SD2_.*", "SD2", spleen$variable)
spleen$variable <- gsub(".*SD4_.*", "SD4", spleen$variable)
spleen$variable <- gsub(".*SD6_.*", "SD6", spleen$variable)
spleen$variable <- gsub(".*SD8_.*", "SD8", spleen$variable)
spleen$variable <- gsub(".*SD10_.*", "SD10", spleen$variable)
spleen$variable <- gsub(".*SD12_.*", "SD12", spleen$variable)

blood$variable <- gsub(".*BD0_.*", "BD0", blood$variable)
blood$variable <- gsub(".*BD2_.*", "BD2", blood$variable)
blood$variable <- gsub(".*BD4_.*", "BD4", blood$variable)
blood$variable <- gsub(".*BD6_.*", "BD6", blood$variable)
blood$variable <- gsub(".*BD8_.*", "BD8", blood$variable)
blood$variable <- gsub(".*BD10_.*", "BD10", blood$variable)
blood$variable <- gsub(".*BD12_.*", "BD12", blood$variable)

###PLOT

p1 <- ggplot(data = spleen, aes(x = variable, y = value, colour = genes)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~ genes) + scale_x_discrete(limits=c('SD0', 'SD2', 'SD4', 'SD6', 'SD8', 'SD10', 'SD12')) + theme_economist() + scale_fill_economist() + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.key.size = unit(1, "cm"), plot.title = element_text(family="Tahoma"), text = element_text(family = "Tahoma"), axis.title = element_text(size = 12), legend.text = element_text(size = 9), legend.title=element_text(face = "bold", size = 9))

p2 <- ggplot(data = blood, aes(x = variable, y = value, colour = genes)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~ genes) + scale_x_discrete(limits=c('BD0', 'BD2', 'BD4', 'BD6', 'BD8', 'BD10', 'BD12')) + theme_economist() + scale_fill_economist() + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.key.size = unit(1, "cm"), plot.title = element_text(family="Tahoma"), text = element_text(family = "Tahoma"), axis.title = element_text(size = 12), legend.text = element_text(size = 9), legend.title=element_text(face = "bold", size = 9))

multiplot(p1, p2, cols = 2)





