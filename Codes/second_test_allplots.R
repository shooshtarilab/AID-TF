library(rGREAT)
library(ggplot2)
library(data.table)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(EnvStats)
library(hrbrthemes)
library(stringr)
library(stringi)
library(ggpubr)
library(grid)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_pubr())
My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 5),
  axis.title.y = element_text(size = 20),
  strip.text.y = element_text(size = 10, angle = 0))

#The address of the final dataframe comprising all results
whole_frame_dir = "~/rprojects/review_data/results/new_generated_res_test/final_dataframe_humna.txt"
whole_frame = read.table(whole_frame_dir,sep = "\t", header = TRUE)

#The file comprising motif names
motif_tf_dir = "~/rprojects/review_data/results/new_generated_res_test/motif_list/motif_name_new.txt"
motif_tf_data = read.table(motif_tf_dir, header = TRUE)

#The main directory of results
main = paste0("~/rprojects/review_data/results/new_generated_res_test/")
dir.create(main)

#List of traits
trait_list = c("IBD","JIA","MS","PSO","RA","T1D","CEL")

#List of cell types 
cell_type_dir = "~/rprojects/results/imm_res_back/temp/list.txt"
cell_type_list = read.table(cell_type_dir, sep = "\t", header = FALSE)
cell_type_list = cell_type_list$V1
cell_type_list = str_remove(cell_type_list, ".bed")
cell_type_list = str_remove(cell_type_list, "e*")


final_obj = list()
final_data_freq = list()
for (i in c(1:length(trait_list))){
  plot_obj = c()
  print(i)
  
  #The directory of the trait of interest
  main_dir = paste0(main,trait_list[i],"/")
  
  #The file having trait-specific results
  direct = paste0(main_dir,"human_tf_table.txt")
  
  data = read.table(direct, sep = "\t", header = TRUE)
  
  motif_list = data$motif
  motif_list = unique(motif_list)
  tf_list = data$TF
  sample_freq = c()
  cell_freq = c()
  show_list = c()
  
  for (j in c(1:length(motif_list))){
    tf_temp = unique(data[data$motif == motif_list[j],]$TF)
    if (length(tf_temp)>1){
      tf_temp = paste(tf_temp, collapse = ",")
    }else{
      tf_temp = tf_temp
    }
    samplefreq_temp = data[data$motif == motif_list[j],]$sample_frequency
    cellfreq_temp = data[data$motif == motif_list[j],]$cell_frequency
    show_temp = paste0(tf_temp,"(",motif_list[j],")")
    show_list = append(show_list,show_temp)
    sample_freq = append(sample_freq,samplefreq_temp[1])
    cell_freq = append(cell_freq,cellfreq_temp[1])
    
  }
  
  sample_temp = rep(1,length(sample_freq))
  combined_data = c(sample_freq,cell_freq)
  combined_label = c(show_list,show_list)
  Frequency_type = c(rep("Sample_frequency",length(show_list)),rep("Cell_frequency",length(show_list)))
  
  data_plot_freq <- data.frame(combined_data,combined_label,Frequency_type)
  data_plot_freq["Trait"] = trait_list[i]
  
  final_data_freq[[i]] = data_plot_freq

}

final_data_freq1 = do.call(rbind,final_data_freq)

final_data_sample = final_data_freq1[final_data_freq1$Frequency_type == "Sample_frequency",]

for (i in c(1:length(trait_list))){
  new_data = final_data_sample[final_data_sample$Trait == trait_list[i],]
  index = order(new_data$combined_data)
  new_data = new_data[index,]
  final_data_sample[final_data_sample$Trait == trait_list[i],] = new_data
  
}


new_obj = list()

new_obj[[1]] = ggplot(data = final_data_sample, aes(y = reorder(combined_label,combined_data), x = combined_data, fill = Trait)) + 
  geom_bar(stat = "identity")+ 
  facet_grid(Trait~., space = 'free_y', scales = 'free_y')+
  xlab("Frequency") + ylab("Motif-TF")+
  ggtitle("Sample Frequency")+
  theme(strip.text.y = element_text(size=0),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text =  element_text(size = 16), strip.background.y = element_blank())
        #                                                          strip.background.x = element_blank())
  #scale_x_discrete(guide = guide_axis(angle = 90)) +  theme(strip.text.x = element_text(size=0),
  #                                                          strip.background.x = element_blank())


#ggsave(ggplot_dir, width = 10, height = 15, limitsize = FALSE, dpi = 300)


final_data_sample = final_data_freq1[final_data_freq1$Frequency_type == "Cell_frequency",]


#ggplot_dir = paste0(main, "freqplot_sample.tiff")
#file.remove(ggplot_dir)


new_obj[[2]] = ggplot(data = final_data_sample, aes(y = reorder(combined_label,combined_data), x = combined_data, fill = Trait)) + 
  geom_bar(stat = "identity")+ 
  facet_grid(Trait~., space = 'free_y', scales = 'free_y')+
  xlab("Frequency") + ylab("Motif-TF")+
  ggtitle("Cell Frequency")+
  theme(strip.text.y = element_text(size=0),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text =  element_text(size = 16), strip.background.y = element_blank())

#Generating Figure 3
ggplot_dir = paste0(main, "freqplot_sample.tiff")
file.remove(ggplot_dir)

ggarrange(new_obj[[1]],new_obj[[2]],labels = c("A","B"),font.label = list(size = 16))


ggsave(ggplot_dir, width = 20, height = 15, limitsize = FALSE, dpi = 300)




plot_obj = list()
heat_obj = list()

#The address of the final dataframe comprising all results
data_list_dir = "~/rprojects/review_data/results/new_generated_res_test/final_dataframe_humna.txt"

data_list = read.table(data_list_dir, sep = "\t", header = TRUE)
motif_list = unique(data_list$motif)
sample_list = unique(data_list$sample_type)

max_val = max(-log10(data_list$motif_pvalue))

final_matrix = list()
final_celltype = c()
final_trait = c()

final_trait_num = c()

for (i in c(1:length(trait_list))){
#for (i in c(1:1)){ 
  print(i)
  #The directory of the trait of interest
  main_dir = paste0(main,trait_list[i],"/")
  
  direct = paste0(main_dir,"human_tf_table.txt")
  
  data = read.table(direct, sep = "\t", header = TRUE)
  
  #motif_list = unique(data$motif)
  sample_list = unique(data$sample_type)
  
  sample_list = str_sort(sample_list)
  
  print(sample_list)
  
  heat_matrix = matrix(0L ,nrow = length(motif_list), ncol = length(sample_list))
  
  for (m in c(1:length(sample_list))){
    for (n in c(1:length(motif_list))){
      sample_temp = data[data$sample_type == sample_list[m],]
      sample_motif = sample_temp[sample_temp$motif == motif_list[n],]
      if (length(sample_motif$motif) > 0){
        temp = -log10(sample_motif$motif_pvalue)
        heat_matrix[n,m] = temp[1]
      }
    }

  }
  motif_cell_index = rowSums(heat_matrix)>0
  heat_matrix_org = heat_matrix
  heat_matrix = as.matrix(heat_matrix[motif_cell_index,])
  show_list = c()
  for (j in c(1:length(motif_list))){
    tf_temp = unique(data_list[data_list$motif == motif_list[j],]$TF)
    if(length(tf_temp)>1){
      tf_temp = paste(tf_temp, collapse = ",")
    }
    show_temp = paste0(tf_temp," (",motif_list[j],")")
    show_list = append(show_list,show_temp)
    
  }
  #colnames(heat_matrix_org) = sample_list
  rownames(heat_matrix_org) = show_list
  show_list = str_remove(show_list,"^/")[motif_cell_index]
  colnames(heat_matrix) = sample_list
  rownames(heat_matrix) = show_list
  data_plot = melt(heat_matrix)
  names(data_plot) = c("motif","sample","enrichment")
  new_motif = c()
  for (j in c(1:length(data_plot$motif))){
    motif_temp = data_plot[j,]$motif
    index = motif_tf_data$Motif_ID == motif_temp
    tf_temp = motif_tf_data[index,]$Motif_name
    tf_temp = tf_temp[1]
    
    new_motif = append(new_motif,paste0(tf_temp," (",motif_temp,")"))
  }
  new_motif = str_remove(new_motif,"^/")
  data_plot$new_motif = new_motif
  cell_type_index = list()
  data_plot['cell_type'] = rep("other",length(data_plot$motif))

  cell_type_org = rep("others",length(sample_list))
  for (j in c(1:length(cell_type_list))){
    temp_cell = cell_type_list[j]
    
    index = grep(temp_cell, sample_list)
    #cell_type_index[[i]] = index
    cell_type_org[index] = rep(temp_cell,length(index))
  }
  
  df = data.frame(cell_type = cell_type_org)
  ha = HeatmapAnnotation(df = df)
  
  #colnames(heat_matrix) = c(1:length(cell_list))
  
  CD14_list = grep("CD14[^0-9]",sample_list)
  CD19_list = grep("CD19[^0-9]",sample_list)
  CD20_list = grep("CD20[^0-9]",sample_list)
  CD3_list = grep("CD3[^0-9]",sample_list)
  CD34_list = grep("CD34[^0-9]",sample_list)
  CD4_list = grep("CD4[^0-9]",sample_list)
  CD56_list = grep("CD56[^0-9]",sample_list)
  CD8_list = grep("CD8[^0-9]",sample_list)
  Thymus_list = grep("Thymus",sample_list)

  
  df["cell_type_immune"] = rep(0,length(df$cell_type))
  df[CD14_list,]$cell_type_immune = rep("CD14",length(CD14_list))
  df[CD19_list,]$cell_type_immune = rep("CD19",length(CD19_list))
  df[CD20_list,]$cell_type_immune = rep("CD20",length(CD20_list))
  df[CD3_list,]$cell_type_immune = rep("CD3",length(CD3_list))
  df[CD34_list,]$cell_type_immune = rep("CD34",length(CD34_list))
  df[CD4_list,]$cell_type_immune = rep("CD4",length(CD4_list))
  df[CD56_list,]$cell_type_immune = rep("CD56",length(CD56_list))
  df[CD8_list,]$cell_type_immune = rep("CD8",length(CD8_list))
  df[Thymus_list,]$cell_type_immune = rep("Thymus",length(Thymus_list))
  
  
  
  new_df = data.frame(matrix(ncol = 2,nrow = length(df$cell_type)))
  
  colnames(new_df) = c("Cell_type","dummy")
  
  #new_df$Cell_type = df$cell_type_immune
  new_df$Cell_type = df$cell_type
  
  #if (nrow(heat_matrix) > 1){
  #  new_index = order(new_df$Cell_type)
  #  new_df = new_df[new_index,]
  #  heat_matrix = heat_matrix[1:nrow(heat_matrix),new_index]
  #}
  
  new_df = subset(new_df, select = -c(dummy))
  col_fun = colorRamp2(c(0, max_val), c("white", "red"))
  #ha = HeatmapAnnotation(df = new_df,
  #                       col = list(Cell_type = c("CD14" = "#a6cee3", "CD19" = "#1f78b4", "CD20" = "#b2df8a",
  #                                                  "CD3" = "#33a02c", "CD34" = "#fb9a99", "CD4" = "#e31a1c",
  #                                                  "CD56" = "#fdbf6f", "CD8" = "#ff7f00", "Thymus" = "#cab2d6")))
  
  annot_stting = gpar(fontface = "bold")
  
  if (trait_list[i] == "IBD" | trait_list[i] == "CEL"){
    annot_stting = gpar(fontface = "bold")
  }
  ha = HeatmapAnnotation(df = new_df,
                         col = list(Cell_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                                  "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                                  "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                                  "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                                  "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                                  "Mobilized.CD8" = "#fc4e2a")),
                         #annotation_name_gp = gpar(fontsize = 18,fontface = "bold"),
                         annotation_name_gp = annot_stting,
                         annotation_label = "Cell Type",
                         show_legend = FALSE)
  
  cn = colnames(heat_matrix)
  colnames(heat_matrix) = c()
  
  #Generating trait-specific TF enrichment heatmaps in Supp. Figure 3
  heat_dir = paste0(main,trait_list[i], "/new_heatmap3.png")
  file.remove(heat_dir)
  #png(heat_dir,width = 1900, height = 1900, res = 300)
  #png(heat_dir, width = 600 + 1400*ncol(heat_matrix)/29, height = 200+1900*nrow(heat_matrix)/40,res = 300)
  png(heat_dir, width = 650 + 1400*ncol(heat_matrix)/29, height = 200+1900*nrow(heat_matrix)/40,res = 300)
  #rownames(heat_matrix) = c()
  colnames(heat_matrix) = c()
  
  
  final_trait_num = append(final_trait_num,ncol(heat_matrix))
  heat_obj[[i]] = Heatmap(heat_matrix, top_annotation = ha, column_order = c(1:length(sample_list)),
                        row_order = c(1:nrow(heat_matrix)),
                        column_title_gp = grid::gpar(fontface = "bold"),
                        
                        #bottom_annotation = HeatmapAnnotation(
                        #  text = anno_text(cn, rot = 90, offset = unit(1, "npc"), just = "right"),
                        #  annotation_height = max_text_width(cn)
                        #),
                        #row_names_gp = grid::gpar(fontsize = 12),
                        
                        
                        column_title = trait_list[i],
                        
                        heatmap_legend_param = list(title = paste0("Enrichment(",trait_list[i],")")), col = col_fun,
                        
                        #width = ncol(heat_matrix), 
                        #height = 1,
                        show_heatmap_legend = FALSE
                        )
  draw(heat_obj[[i]])
  #plot(1)
  #date_time<-Sys.time()
  #while((as.numeric(Sys.time()) - as.numeric(date_time))<5){}
  dev.off()
  
  final_matrix[[i]] = heat_matrix_org
  final_trait = append(final_trait, rep(trait_list[i],ncol(heat_matrix)))
  final_celltype = append(final_celltype, df$cell_type)
  
}


lgd1 = Legend(col_fun = col_fun, title = "Enrichment",direction = "vertical",
              #labels_gp = gpar(fontsize = 28),
              title_gp = gpar(fontface = "bold"))

lgd2 = Legend(labels = c("CD14.Primary", "CD19.Primary", "CD20.Primary",
                         "CD3.Primary", "CD34.Primary", "CD4.Primary",
                         "CD56.Primary" , "CD8.Primary" , "Thymus" ,
                         "CD3.Cord" , "Mobilized.CD3" ,
                         "Mobilized.CD34" , "Mobilized.CD4" , "Mobilized.CD56", 
                         "Mobilized.CD8"),
              legend_gp = gpar(fill = c("#a6cee3", "#1f78b4", "#b2df8a",
                          "#33a02c", "#fb9a99","#e31a1c",
                           "#fdbf6f", "#ff7f00", "#cab2d6",
                           "#00441b", "#99d8c9",
                            "#e7298a", "#b30000", "#ffffcc", 
                            "#fc4e2a")),
              #labels_gp = gpar(fontsize = 28),
              title_gp = gpar(fontface = "bold"),
              #legend_gp = gpar(fill = 1:16),
              title = "Cell Type", direction = "vertical",ncol = 1)

pd = packLegend(lgd1, lgd2)

heat_dir = paste0(main, "/anot_heatmap.png")
file.remove(heat_dir)
png(heat_dir,width = 500, height = 1200, res = 300)
draw(pd)
dev.off()



lgd1 = Legend(col_fun = col_fun, title = "Enrichment",direction = "vertical",
              labels_gp = gpar(fontsize = 18),
              title_gp = gpar(fontsize = 18, fontface = "bold"))

lgd2 = Legend(labels = c("CD14.Primary", "CD19.Primary", "CD20.Primary",
                         "CD3.Primary", "CD34.Primary", "CD4.Primary",
                         "CD56.Primary" , "CD8.Primary" , "Thymus" ,
                         "CD3.Cord" , "Mobilized.CD3" ,
                         "Mobilized.CD34" , "Mobilized.CD4" , "Mobilized.CD56", 
                         "Mobilized.CD8"),
              legend_gp = gpar(fill = c("#a6cee3", "#1f78b4", "#b2df8a",
                                        "#33a02c", "#fb9a99","#e31a1c",
                                        "#fdbf6f", "#ff7f00", "#cab2d6",
                                        "#00441b", "#99d8c9",
                                        "#e7298a", "#b30000", "#ffffcc", 
                                        "#fc4e2a")),
              labels_gp = gpar(fontsize = 18),
              title_gp = gpar(fontsize = 18,fontface = "bold"),
              #legend_gp = gpar(fill = 1:16),
              title = "Cell_type", direction = "vertical")

pd = packLegend(lgd1, lgd2)


final_matrix_data = do.call(cbind,final_matrix)



new_frame = data.frame(matrix(nrow = length(final_trait),ncol = 1))
colnames(new_frame) = c("Cell_type")

new_frame$Cell_type = final_celltype

ha = HeatmapAnnotation(df = new_frame,
                       col = list(Cell_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                                "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                                "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                                "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                                "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                                "Mobilized.CD8" = "#fc4e2a")),
                       annotation_name_gp = gpar(fontsize = 18,fontface = "bold"),
                       annotation_label = c("Cell Type"),
                       #show_legend = c(FALSE,FALSE),
                       annotation_legend_param = list(labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20,fontface = "bold"))
                       
                       
)

lgd1 = Legend(labels = c("IBD", "JIA", "MS",
                         "PSO", "RA", "T1D","CEL"),
              legend_gp = gpar(fill = c("#7fc97f", "#beaed4", "#fdc086",
                                        "#ffff99", "#386cb0", "#f0027f","#bf5b17")),
              labels_gp = gpar(fontsize = 18),
              title_gp = gpar(fontsize = 20,fontface = "bold"),
              #legend_gp = gpar(fill = 1:16),
              title = "Trait", direction = "vertical")


lgd2 = Legend(labels = c("CD14.Primary", "CD19.Primary", "CD20.Primary",
                         "CD3.Primary", "CD34.Primary", "CD4.Primary",
                         "CD56.Primary" , "CD8.Primary" , "Thymus" ,
                         "CD3.Cord" , "Mobilized.CD3" ,
                         "Mobilized.CD34" , "Mobilized.CD4" , "Mobilized.CD56", 
                         "Mobilized.CD8"),
              legend_gp = gpar(fill = c("#a6cee3", "#1f78b4", "#b2df8a",
                                        "#33a02c", "#fb9a99","#e31a1c",
                                        "#fdbf6f", "#ff7f00", "#cab2d6",
                                        "#00441b", "#99d8c9",
                                        "#e7298a", "#b30000", "#ffffcc", 
                                        "#fc4e2a")),
              labels_gp = gpar(fontsize = 18),
              title_gp = gpar(fontsize = 20,fontface = "bold"),
              #legend_gp = gpar(fill = 1:16),
              title = "Cell Type", direction = "vertical")

pd = packLegend(lgd1, lgd2)

col_fun = colorRamp2(c(0, max(final_matrix_data)), c("#e0e0e0", "#67001f"))

final_heat = Heatmap(final_matrix_data, top_annotation = ha, column_order = c(1:length(final_trait)),
                        #column_names_gp = grid::gpar(fontsize = 9),
                        column_title_gp = grid::gpar(fontsize = 18,fontface = "bold"),
                        rect_gp = gpar(col= "black"),
                        
                        #bottom_annotation = HeatmapAnnotation(
                        #  text = anno_text(cn, rot = 90, offset = unit(1, "npc"), just = "right"),
                        #  annotation_height = max_text_width(cn)
                        #),
                        row_names_gp = grid::gpar(fontsize = 24),
                        row_names_rot = 0,
                        #Legend(title_gp = gpar(fontsize = 32),labels_gp = gpar(fontsize = 32), col = col_fun,
                        #column_title = trait_list[i],
                        heatmap_legend_param = list(title = paste0("Enrichment"),
                                                    labels_gp = gpar(fontsize = 18),
                                                    title_gp = gpar(fontsize = 20,fontface = "bold")
                                                    ), col = col_fun,
                        column_split = final_trait,
                        column_gap = unit(4, "mm")
                     #show_heatmap_legend = FALSE,
                        #width = 5000000000000000000000, 
                        #height = 1,
                     )


#Generating Figure 2
heat_dir = paste0(main, "/final_heatmap_newnew.png")
file.remove(heat_dir)
png(heat_dir,width = 16800, height = 8600, res = 300)
#png(heat_dir,width = 12800, height = 8600, res = 300)
draw(final_heat)

dev.off()

