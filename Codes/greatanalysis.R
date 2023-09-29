library(rGREAT)
library(ggplot2)
library(data.table)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(EnvStats)
library(ggforce)
library(stringr)
library(ggpubr)
library(tidytext)
library(dplyr)
theme_set(theme_pubr())

#Main directory of AID-related results
main = paste0("~/rprojects/results/immunochip_res_new/temp/new_generated_res/")

#List of traits
trait_list = c("IBD","JIA","MS","PSO","RA","T1D","CEL")

final_obj = list()
final_obj1 = list()
final_data = list()

for (i in c(1:length(trait_list))){
  plot_obj = c()
  print(i)
  main_dir = paste0(main,trait_list[i],"/")
  
  #The table comprising AID-specific results
  direct = paste0(main_dir,"final_total_list.txt")
  
  data = read.table(direct, sep = "\t", header = TRUE)
  data = select(data, c("chr","start","end"))
  data = unique(data)
  
  job = submitGreatJob(data)
  
  tb = getEnrichmentTables(job)
  
  head(tb[[2]])
  
  interest_obj = tb[[2]]
  
  plot_obj = data.frame(matrix(nrow = 20, ncol = 2))
  colnames(plot_obj) = c("bio_function","Pvalue")
  plot_obj$bio_function = interest_obj$name[1:20]
  plot_obj$Pvalue = interest_obj$Binom_Adjp_BH[1:20]
  final_obj[[i]] = ggtexttable(plot_obj, rows = NULL, theme = ttheme("mOrange"))
  tab_add_title(final_obj[[i]], trait_list[i])
  
  plot_obj$enrichment = -log10(plot_obj$Pvalue)
  thresh_index = plot_obj$enrichment > (-log10(0.05))
  plot_obj = plot_obj[thresh_index,]
  final_obj1[[i]] = ggplot(plot_obj, aes(y=bio_function, x=enrichment)) + 
    geom_bar(stat = "identity") +
    ylab("biological function") + xlab("-log10(Binomial adjusted P_value)")+
    ggtitle(paste0("great pathway results for: ",trait_list[i]))
    
  plot_obj["trait"] = rep(trait_list[i],length(plot_obj$Pvalue))
  
  final_data[[i]] = plot_obj
}


#Generating Figure 6
final_data = do.call(rbind,final_data)
index = grep("^vesicle",final_data$bio_function)
final_data[index,]$bio_function = "vesicle fusion with endoplasmic membrane"
final_data["Trait"] = final_data$trait

#The directory of the figure 6
ggplot_dir = "~/rprojects/results/immunochip_res_new/temp/new_generated_res/great.tiff"
file.remove(ggplot_dir)
ggplot(data = final_data, aes(y = reorder_within(bio_function,enrichment,Trait), x = enrichment, fill = Trait)) + 
  geom_bar(stat = "identity")+ 
  scale_y_reordered()+
  ggforce::facet_col(vars(trait),scales = "free",space = "free",strip.position = "right")+
  xlab("-log10(P value)") + ylab("Biological Pathways")+
  ggtitle("Great Pathway Results")+ #theme_bw()+
  theme(strip.text.y = element_text(size=20),axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),legend.title = element_text(size = 20),
        legend.text =  element_text(size = 18),
        strip.text.x = element_text(size=16),
        axis.text.y = element_text(size = 18),title = element_text(size = 20))+
  
  scale_fill_manual(breaks = c("IBD","JIA","MS","PSO","RA","T1D","CEL"),
                    values = c("#7fc97f","#beaed4","#fdc086","#fee391","#386cb0","#f0027f","#bf5b17"))
ggsave(ggplot_dir, width = 30, height = 32, limitsize = FALSE, dpi = 300)


