# R script to generate simple performance charts of msprime VS simbac
# create pdf figure with command line:
# Rscript plot_performance_gene_conversion.R


df = read.csv(file = "../../data/simbac.csv")

library(ggplot2)


pdf(file = "gc_perf.pdf", width = 6, height = 3)
ggplot(data = df, aes(x = sample_size, y = time, color = tool)) +
  geom_point() + 
  labs(x = "sample size") +
  xlim(0,500)

    #coord_trans(y="log2") + 
  #scale_x_continuous("sample size", breaks = c(10,50,100,500,1000,4000,10000), labels = c("10","50","100", "500","1000", "4000", "10000"), c(0,10000), trans = "log2")
  #scale_x_continuous("sample size", breaks = c(10,50,100,500,1000), labels = c("10","50","100", "500","1000"), c(0,1000))

  
dev.off()


pdf(file = "gc_perf_log.pdf", width = 6, height = 3)
ggplot(data = df_long, aes(x = sample_size, y = time, color = tool)) +
  geom_point() +
  labs(x = "sample size") +
  xlim(0,500) +
  #coord_trans(y="log2")
  scale_y_continuous("time", breaks = c(10,60,120,240,600,1200,2400,3600), trans = "log2" )

    #coord_trans(y="log2") + 
  #scale_x_continuous("sample size", breaks = c(10,50,100,500,1000,4000,10000), labels = c("10","50","100", "500","1000", "4000", "10000"), c(0,10000), trans = "log2")
  #scale_x_continuous("sample size", breaks = c(10,50,100,500,1000), labels = c("10","50","100", "500","1000"), c(0,1000))


dev.off()

