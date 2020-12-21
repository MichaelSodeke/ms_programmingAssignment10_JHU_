# [my_theme]

# reference: https://github.com/tidyverse/ggplot2/blob/master/R/theme-elements.r
#            https://ggplot2.tidyverse.org/reference/theme.html  

library(ggplot2)
library(gridExtra)

my_theme = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%

    theme(
      # Specify axis options
      axis.line.x = element_line(size = 0.1, colour = "grey20"), 
      axis.line.y = element_line(size = 0.1, colour = "grey20"),  
      axis.ticks.x = element_line(size=0.1),
      axis.ticks.y = element_line(size=0.1),
      axis.ticks.length = unit(0.3, "lines"),
      axis.text.x = element_text(color="grey20", face ="plain", size=12, family="Segoe UI"), 
      axis.text.y = element_text(color="grey20", face ="plain", size=12, family="Segoe UI"),  
      axis.title.x=element_text(color="grey20", face="plain", size=15, margin = margin(0, 10, 10, 0)),
      axis.title.y=element_text(color="grey20", face="plain", size=15, margin = margin(0, 10, 0, 0), angle = 90),  
      # Specify legend options
      legend.box = "horizontal",
      legend.box.just = c("right", "top"), 
      # Specify panel options
      panel.background = element_rect(fill = "white", color=NULL),
      panel.border = element_rect(fill = NA, color = "white", linetype = "blank"),  
      panel.grid = element_blank(),  
      panel.spacing = unit(0.5, "lines"),  
      # Specify facetting options
      strip.background = element_rect(fill = "grey20", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(face = "bold", size = base_size*1.5, color = "grey20", vjust = 1, hjust=0),  
      plot.subtitle = element_text(size = base_size*1.0, color = "grey20", vjust = 1, hjust=0), 
      plot.caption = element_text(size = base_size*1.0, color = "grey20", vjust = 1, hjust=0),    
      plot.margin = unit(rep(1, 4), "lines")   
    ) 
}