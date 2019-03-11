require(vegan)
require(dplyr)
require(reshape2)
require(ggplot2)
require(ggthemes)
require(cowplot)
require(grid)
require(codyn)

### in order to run this the sesim model has to have been run already and the output be in local memory
### you also need to have the microcosm data saved in a data/ folder in the working directory

cell_mass <- c(2.227e-7,8.05e-8,1.52e-8,4.30e-8,1.96e-7,3.76e-6,4.77e-9,6.9e-8,9.68e-8)
productivity <- function(mat){
  prod <- sum(mat * cell_mass)
  return (prod)
}

df <- read.csv("data/Microcosm_RawData_CF_2017.csv")
names(df)[1] <- "treatment"
df <- df[-13] # import the raw data

# setup of metric storage for lab data 
metrics <- data.frame(matrix(NA, nrow = 240, ncol = 7))
colnames(metrics) <- c("treatment", "replicate", "time", "div", "rich", "eve", "prod")
metrics[1:3] <- df[1:3]
metrics$treatment <- factor(metrics$treatment, levels = c("low", "intermediate", "high"))

# calculate each metric
for(i in unique(metrics$time)){
  for(j in unique(metrics$replicate)){
    metrics[metrics$replicate == j & metrics$time == i,]$div =
      apply(df[df$replicate == j & df$time == i,-(1:3)], 1, vegan::diversity)
    metrics[metrics$replicate == j & metrics$time == i,]$rich =
      apply(df[df$replicate == j & df$time == i,-(1:3)], 1, vegan::specnumber)
    metrics[metrics$replicate == j & metrics$time == i,]$prod =
      apply(df[df$replicate == j & df$time == i,-(1:3)], 1, productivity)
    metrics[metrics$replicate == j & metrics$time == i,]$eve =
      metrics[metrics$replicate == j & metrics$time == i,]$div /
      log(metrics[metrics$replicate == j & metrics$time == i,]$rich)
  }
}

### HELPER FUNCTIONS FOR PLOTTING ###
add_labs = function (plot, xlab='', ylab='') {
  ylab_obj = textGrob(ylab, gp=gpar(fontfamily='serif', fontsize=16), rot=90)
  xlab_obj = textGrob(xlab, gp=gpar(fontfamily='serif', fontsize=16))
  labelled_plot = grid.arrange(arrangeGrob(plot, left=ylab_obj, bottom=xlab_obj))
  return (labelled_plot)
}

reshape_dat <- function(df, subset_name, subset_col='treatment', melt_vars=c('time', 'replicate'))  {
  sub <- df[df[subset_col] == subset_name,]
  sub_drop <- dplyr::select(sub, -subset_col)
  sub_tran <- sub_drop
  sub_tran[species_names] <- sqrt(sub_tran[species_names] + 1.0)
  melted_subset <- melt(sub_tran,id.vars=melt_vars,na.rm=TRUE)
  return (melted_subset)
}

### FIGURE 2 ###
plot_individual_box <- function(frame, yname, ylab, xname='treatment', scientific_y=FALSE, theme=theme_tufte()){
  plot = ggplot(data=frame, aes_string(x=xname, y=yname)) + 
    geom_boxplot() + xlab('') + ylab(ylab) + 
    theme + theme(axis.line = element_line(color = 'black'), text = element_text(size=16))
  if (scientific_y)
    plot = plot + scale_y_continuous(labels = function(y) format(y, scientific=TRUE))
  return (plot)
}

boxplot_data = metrics[metrics$time==30,]
div = plot_individual_box(boxplot_data, 'div', 'Shannon diversity')
rich = plot_individual_box(boxplot_data, 'rich', 'Species richness')
eve = plot_individual_box(boxplot_data, 'eve', 'Pielou\'s evenness')
prod = plot_individual_box(boxplot_data, 'prod', 'Productivity', scientific_y=TRUE)

fig_2_no_lab = plot_grid(div, rich, eve, prod, nrow=2, align='hv', labels=c('a)', 'b)', 'c)', 'd)'))
fig_2 = add_labs(fig_2_no_lab, xlab='Treatment')
rm(fig_2_no_lab)
ggsave('plots/fig_2.pdf', fig_2, 
       width=6, height=6, units='in' ,dpi=300)

### FIGURE 4 ###
assembly_data <- na.omit(df)
species_names <- c('B', 'E', 'C', 'Pa', 'Pb', 'S', 'T', 'L', 'V')
assembly_data[species_names] <- sqrt(assembly_data[species_names] + 1.0)

assembly.trajectory <- function(data, level, col, limitx = c(-0.3,0.5), limity = c(-0.4,0.2), theme = theme_tufte()){
  dist <- vegdist(data[-(1:3)],method="bray")
  disp <- betadisper(dist,group=data$treatment,sqrt.dist=F,type="centroid")
  temp <- cbind(data[1:3],scores(disp)$sites)
  plot <- ggplot(data=temp[temp$treatment==level,],aes(x=PCoA1,y=PCoA2,group=replicate,by=time)) +
    geom_path(colour=col,size=1) + geom_point(data=temp[temp$treatment==level & temp$time==30,],size=5) +
    geom_point(data=temp[1,], size = 5, shape = 17) + xlim(limitx) + ylim(limity) +
    theme + theme(axis.line = element_line(color = 'black'), text = element_text(size=16)) +
    xlab('') + ylab('')
  return (plot)
}    

low = assembly.trajectory(assembly_data, 'low', 'darkslategray2')
int = assembly.trajectory(assembly_data, 'intermediate', 'darkslategray')
high = assembly.trajectory(assembly_data, 'high', 'black')

fig_4_no_lab = plot_grid(low, int, high, nrow=1, labels=c('a)', 'b)', 'c)'))
fig_4 = add_labs(fig_4_no_lab, 'PCoA1', 'PCoA2')
rm(fig_4_no_lab)
ggsave('plots/fig_4.pdf', fig_4,
       width=8, height=3, units='in', dpi=300)

### FIGURE 5 ###
calculate_synchrony <- function(df, treatments){
  all_synchrony_output <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("replicate", "synchrony", "treatment"))
  for (treatment in treatments){
    synch_data <- reshape_dat(df, treatment)
    synch_output <- synchrony(synch_data,time.var="time",species.var="variable",abundance.var="value",
                              metric="Loreau",replicate.var="replicate")
    synch_output$treatment <- treatment
    all_synchrony_output <- rbind(all_synchrony_output, synch_output)
  }
  return(all_synchrony_output)
}

treatment_levels = c('low', 'intermediate', 'high')
community_synchrony <- calculate_synchrony(df, treatment_levels)
# anova(lm(synchrony ~ treatment,data=community_synchrony))
# TukeyHSD(aov(lm(synchrony ~ treatment,data=community_synchrony)))

plot_synchrony <- function(synch_data, xname='treatment', yname='synchrony', xlab='Treatment', ylab='Loreau\'s community synchrony', theme=theme_tufte()){
  synchrony_plot <- ggplot(synch_data, aes_string(x=xname, y=yname)) + geom_boxplot() +
    theme + theme(axis.line = element_line(color = 'black'), text = element_text(size=16)) +
    xlab(xlab) + ylab(ylab)
  return (synchrony_plot)
}

fig_5 <- plot_synchrony(community_synchrony)
ggsave('plots/fig_5.pdf', fig_5, 
       width=6, height=6, units='in', dpi=300)

### FIGURE 9 ###
last_step <- max(rep_samples$time)
grouping_metric_cols_replicate <- c('replicate', 'patch', 'connectivity', 'quality')
grouping_metric_cols <- c('patch', 'connectivity', 'quality')

find_abundance_cv <- function(df, timestep=last_step, metric_cols_rep=grouping_metric_cols_replicate, 
                              metric_cols=grouping_metric_cols, species_cols=species) {
  cv_data <- df[df$time == timestep,]
  cv_data <- cv_data[c(metric_cols_rep, species_cols)]
  cv_data_tran <- cv_data
  cv_data_tran[species_cols] <- decostand(cv_data_tran[species_cols], 'hellinger', na.rm=TRUE)
  cv_sumdata <- cv_data_tran %>%
    group_by(.dots=metric_cols_rep) %>%
    transform(mean=rowMeans(cv_data_tran[species_cols])) %>%
    transform(sd=apply(cv_data_tran[species_cols], 1, sd))
  cv_sumdata <- dplyr::select(cv_sumdata, -species_cols)
  cv_sumdata$cv <- cv_sumdata$sd / cv_sumdata$mean
  plotdata <- cv_sumdata %>% 
    group_by(.dots=metric_cols) %>%
    summarize(mean=mean(cv), sd=sd(cv))
  plotdata$ymax <- plotdata$mean + plotdata$sd
  plotdata$ymin <- plotdata$mean - plotdata$sd
  plotdata[metric_cols] <- lapply(plotdata[metric_cols], factor)
  return (plotdata)
}

plotdata <- find_abundance_cv(rep_samples)

plot_lines <- function(df, xname='patch', yname='mean', groupname='quality', ymin='ymin', ymax='ymax', theme=theme_tufte(),
                       get_legend=FALSE, legend_title='Quality') {
  if (get_legend==FALSE){
    cv_plot <- ggplot(df, aes_string(x=xname, y=yname, group=groupname, colour=groupname)) + 
      geom_line(size=1, linetype='twodash') + geom_ribbon(aes_string(ymin=ymin, ymax=ymax, fill=groupname), alpha=0.3, size=0) + 
      theme + theme(axis.line = element_line(color = 'black'), text = element_text(size=16), legend.position='None') +
      xlab('') + ylab('')
    return (cv_plot)  
  }
  if (get_legend==TRUE){
    temp_plot <- ggplot(df, aes_string(x=xname, y=yname, group=groupname, colour=groupname)) + geom_line(size=2) + 
      labs(color=legend_title)  + theme
    legend <- cowplot::get_legend(temp_plot)
    return (legend)  
  }
}

absent <- plot_lines(plotdata[plotdata$connectivity==1,])
linear <- plot_lines(plotdata[plotdata$connectivity==2,])
circular <- plot_lines(plotdata[plotdata$connectivity==3,])
global <- plot_lines(plotdata[plotdata$connectivity==4,])
legend <- plot_lines(plotdata[plotdata$connectivity==1,], get_legend=TRUE)

fig_9_no_lab = plot_grid(absent, linear, legend, circular, global, nrow=2, labels=c('a)', 'b)', '', 'c)', 'd)'),
                         rel_widths=c(1,1,0.15,1,1))
fig_9 = add_labs(fig_9_no_lab, xlab='Patch', ylab='Coefficient of variation of abundance')
rm(fig_9_no_lab)
ggsave('plots/fig_9.pdf', fig_9, 
       width=6, height=6, units='in', dpi=300)





