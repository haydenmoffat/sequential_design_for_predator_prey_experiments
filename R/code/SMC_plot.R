#-----------------------------------------------------------------------------------------------#
# Script Name: SMC_plot                                                                         #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This R script plots the marginal posterior distributions for each of the paramters after the  #
# sequential experimental design algorithm is complete.                                         #
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#


# Plot marginal posterior distributions
plot_data <- theta[,,M] %>% 
  exp %>% 
  as.data.frame()

names(plot_data) <- c('a', 'Th', 'lambda')

# Plot for a
figure_A <- ggplot(plot_data, aes(a)) +
  geom_density(lwd = 1.1, fill = 'blue', alpha = 0.3) +
  xlab(TeX('$a$')) +
  theme(axis.title=element_text(size=20,face="bold"))

# Plot for Th
figure_B <- ggplot(plot_data, aes(Th)) +
  geom_density(lwd = 1.1, fill = 'red', alpha = 0.3) +
  xlab(TeX('$T_{h}$')) +
  theme(axis.title=element_text(size=20,face="bold"))

# Plot for lambda
figure_C <- ggplot(plot_data, aes(lambda)) +
  geom_density(lwd = 1.1, fill = 'green', alpha = 0.3) +
  xlab(expression(lambda)) +
  theme(axis.title=element_text(size=20,face="bold"))

if (models[M] %in% c(1,2)){
  plot <- ggarrange(figure_A, figure_B, figure_C, nrow = 1) %>% 
    annotate_figure(top = text_grob(paste("Marginal posterior distributions for ", modeltype[[models[M]]]), 
                                    color = "black", face = "bold", size = 14))
} else{
  plot <- ggarrange(figure_A,figure_B, nrow = 1) %>% 
    annotate_figure(top = text_grob(paste("Marginal posterior distributions for", modeltype[[models[M]]]), 
                                    color = "black", face = "bold", size = 14))
}

print(plot)