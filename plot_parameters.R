## Parameters and functions used for plotting and calculations as of April 4th

cell_mass <- c(2.227e-7,8.05e-8,1.52e-8,4.30e-8,1.96e-7,3.76e-6,4.77e-9,6.9e-8,9.68e-8)
productivity <- function(mat){
  prod <- sum(mat * cell_mass)
  return (prod)
}

fontfam = 'sans'
def_label_fontface = 'plain'
grid_labels = c('a)', 'b)', 'c)', 'd)')
grid_labels_9 = c('a)', 'b)', 'c)', 'd)')
def_label_size = 12

# global_theme = theme_tufte() + theme(axis.line=element_line(color='black'), 
#                                   text=element_text(size=16)) # old default

scaleFUN <- function(x) sprintf("%.2f", x)

global_theme = theme_gray() + theme(panel.background=element_rect(fill="white"), panel.border=element_rect(linetype="solid", fill=NA, size=0.3),
                             axis.ticks=element_line(size=0.3), axis.ticks.length=unit(0.2, "cm"),
                             axis.title.y=element_text(size=rel(1.2)), axis.title.x=element_text(size=rel(1.2)),
                             axis.text.y=element_text(size=rel(1.1)), axis.text.x=element_text(size=rel(1.1)),
                             strip.text=element_text(size=rel(1)), strip.background=element_rect(fill="white"),
                             strip.text.x=element_text(angle = 0, hjust = 0))
