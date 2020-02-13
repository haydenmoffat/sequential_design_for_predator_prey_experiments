#-----------------------------------------------------------------------------------------------#
# Script Name: utility_plot                                                                     #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This R script plots the utilities for the different designs.                                  #                           
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#

# Plot utility curve
util_plot <- ggplot() + 
  geom_point(aes(x = 1:length(utility),y = utility),colour="black", shape=21, lwd=2, fill = "blue") +
  xlab("Design")+
  ylab("Utility value")+
  theme_pander()+
  ggtitle(paste("Utility Curve - Iteration ", i))

print(util_plot)