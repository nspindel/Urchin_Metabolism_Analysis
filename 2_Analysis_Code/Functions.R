###############################################################################
# Utility functions                                                           #
#                                                                             #
# Authors: N.B. Spindel & D.K. Okamoto                                        #
###############################################################################

# Function to set up consistent ggplot aesthetics:
gg_options <-  function(){
  theme_bw(base_size = 16)+theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.background =  element_blank(),
    legend.background =  element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank())}

# Function to format a number to include three decimal places to the right of zero when printing
round_for_print <- function(x) print(format(round(x, digits = 3), nsmall = 3), quote = FALSE)