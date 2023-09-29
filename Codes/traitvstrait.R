library(ggplot2)
library(data.table)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(EnvStats)
library(hrbrthemes)
library(stringr)
library(ggpubr)
theme_set(theme_pubr())

My_Theme = theme(
  axis.title.x = element_text(size = 10, hjust = 0.5),
  axis.title.y = element_text(size = 10, hjust = 0.5))

#Address of the final dataframe having all results
data_dir = "~/rprojects/results/imm_res_back/temp/new_generated_res/final_dataframe.txt"
data = read.table(data_dir, sep = "\t", header = TRUE)

trait_list = unique(data$trait)
#List of traits
trait_list = trait_list = c("ATD","CEL","PBC","IBD","JIA","MS","PSO","RA","T1D")
motif_list = list()

#Getting the list of all affected motifs
for (i in c(1:length(trait_list))){
  motif_temp = data[data$trait == trait_list[i],]$motif
  motif_temp = unique(motif_temp)
  motif_list[[i]] = motif_temp
}

heat_matrix = matrix(0L ,nrow = length(trait_list), ncol = length(trait_list))

for (m in c(1:length(trait_list))){
  for (n in c(1:length(trait_list))){
    if (m == n){
      heat_matrix[m,n] = 0
    }else{
      heat_matrix[m,n] = length(intersect(motif_list[[m]],motif_list[[n]]))    
    }
  }
  
}

rownames(heat_matrix) = trait_list
colnames(heat_matrix) = trait_list

data_plot = melt(heat_matrix)
names(data_plot) = c("trait1","trait2","overlap_number")

plot_o = ggplot(data_plot, aes(trait1, trait2, fill= overlap_number)) + 
  geom_tile()+
  xlab("Trait1")+ylab("Trait2")+ggtitle("Number of Common Motifs")+
  #scale_fill_gradient(low="white", high="red")+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_fill_discrete(name = "Number of Overlapping Motifs")+
  scale_fill_gradientn(breaks = c(0,1,2),colours = c("light blue","dark blue"))+
  #scale_fill_distiller(palette = "heat") +
  theme_ipsum() + My_Theme

save_file =  "~/rprojects/results/imm_res_back/temp/new_generated_res/traitvstrait_motif.tiff"

file.remove(save_file)

#Generating the figure showing the number of common affected motifs 
#between traits (Figure 4.A) 
ggplot(data_plot, aes(trait1, trait2, fill= overlap_number)) + 
  geom_tile()+
  xlab("Trait1")+ylab("Trait2")+ggtitle("Number of Common Motifs")+
  #scale_fill_gradient(low="white", high="red")+
  #scale_fill_continuous(breaks = c(0,1,2))+
  scale_fill_gradient(low = "#c6dbef",high = "#08306b",guide = "colorbar",breaks = c(2,1,0))+
  
  guides(fill=guide_legend(title='Number of Overlapping Motifs'))+
  #scale_fill_distiller(palette = "heat") +
  theme_ipsum() + My_Theme

ggsave(save_file,dpi = 1000)



sample_list = list()

for (i in c(1:length(trait_list))){
  sample_temp = data[data$trait == trait_list[i],]$sample_type
  sample_temp = unique(sample_temp)
  sample_list[[i]] = sample_temp
}
#sample_list[[2]] = character(0)
heat_matrix = matrix(0L ,nrow = length(trait_list), ncol = length(trait_list))

for (m in c(1:length(trait_list))){
  for (n in c(1:length(trait_list))){
    if (m == n){
      heat_matrix[m,n] = 0
    }else{
      heat_matrix[m,n] = length(intersect(sample_list[[m]],sample_list[[n]]))    
    }
  }
  
}

rownames(heat_matrix) = trait_list
colnames(heat_matrix) = trait_list

data_plot = melt(heat_matrix)
names(data_plot) = c("trait1","trait2","overlap_number")

plot_o = ggplot(data_plot, aes(trait1, trait2, fill= overlap_number)) + 
  geom_tile()+
  xlab("Trait1")+ylab("Trait2")+ggtitle("Number of common samples")+
  #scale_fill_gradient(low="white", high="red")+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  #scale_fill_distiller(palette = "heat") +
  theme_ipsum() + My_Theme

save_file =  "~/rprojects/results/imm_res_back/temp/new_generated_res/traitvstrait_sample.tiff"

file.remove(save_file)

#Generating the figure showing the number of common affected samples 
#between traits 

ggplot(data_plot, aes(trait1, trait2, fill= overlap_number)) + 
  geom_tile()+
  xlab("Trait1")+ylab("Trait2")+ggtitle("Number of common samples")+
  scale_fill_gradient(low = "#c6dbef",high = "#08306b",guide = "colorbar")+
  
  #scale_fill_gradient(low="white", high="red")+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  #scale_fill_distiller(palette = "heat") +
  theme_ipsum() + My_Theme

ggsave(save_file,dpi = 1000)



cell_list = list()

for (i in c(1:length(trait_list))){
  cell_temp = data[data$trait == trait_list[i],]$cell_type
  cell_temp = unique(cell_temp)
  cell_list[[i]] = cell_temp
}
#cell_list[[2]] = character(0)
heat_matrix = matrix(0L ,nrow = length(trait_list), ncol = length(trait_list))

for (m in c(1:length(trait_list))){
  for (n in c(1:length(trait_list))){
    if (m == n){
      heat_matrix[m,n] = 0
    }else{
      heat_matrix[m,n] = length(intersect(cell_list[[m]],cell_list[[n]]))    
    }
  }
  
}

rownames(heat_matrix) = trait_list
colnames(heat_matrix) = trait_list

data_plot = melt(heat_matrix)
names(data_plot) = c("trait1","trait2","overlap_number")

plot_o = ggplot(data_plot, aes(trait1, trait2, fill= overlap_number)) + 
  geom_tile()+
  xlab("Trait1")+ylab("Trait2")+ggtitle("Number of common cell-types")+
  #scale_fill_gradient(low="white", high="red")+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  #scale_fill_distiller(palette = "heat") +
  theme_ipsum() + My_Theme

save_file =  "~/rprojects/results/imm_res_back/temp/new_generated_res/traitvstrait_cell.tiff"

file.remove(save_file)

#Generating the figure showing the number of common affected cell types 
#between traits (Figure 4.B) 
ggplot(data_plot, aes(trait1, trait2, fill= overlap_number)) + 
  geom_tile()+
  xlab("Trait1")+ylab("Trait2")+ggtitle("Number of Common Cell Types")+
  #scale_fill_gradient(low="white", high="red")+
  #scale_fill_continuous(breaks = c(0,1,2))+
  scale_fill_gradient(low = "#c6dbef",high = "#08306b",guide = "colorbar",breaks = c(6:0))+
  
  guides(fill=guide_legend(title='Number of Overlapping Cell Types'))+
  #scale_fill_distiller(palette = "heat") +
  theme_ipsum() + My_Theme

ggsave(save_file,dpi = 1000)
