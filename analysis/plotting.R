library(ggplot2)

scatterplot_log_log <- function(data, xname, yname, filename = "") {
  f <- ggplot(data, aes_string(x = xname, y = yname)) + 
         scale_y_log10() + scale_x_log10() + 
         geom_point(alpha = 0.45, color = "red");

  f <- f + theme_bw() + theme(legend.position = "none");
  #  theme(axis.text =  element_text(size = rel(0.85)))
  f <- f + theme(panel.border = element_rect(colour="white"))
  #f <- f + theme(strip.background = element_rect(colour="white", fill="white"))
  f;

  #if(length(filename) > 0) {
  #  ggsave(filename, width = 7, height = 4, paper = 'special')
  #}
  #f;
}

scatterplot_lin_log <- function(data, xname, yname, filename = "") {
  f <- ggplot(data, aes_string(x = xname, y = yname)) + 
         scale_y_log10() + 
         geom_point(alpha = 0.45, color = "red");

  f <- f + theme_bw() + theme(legend.position = "none");
  #  theme(axis.text =  element_text(size = rel(0.85)))
  f <- f + theme(panel.border = element_rect(colour="white"))
  #f <- f + theme(strip.background = element_rect(colour="white", fill="white"))
  f;

  #if(length(filename) > 0) {
  #  ggsave(filename, width = 7, height = 4, paper = 'special')
  #}
  #f;
}
