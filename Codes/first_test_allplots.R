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
library(ComplexHeatmap)
library(circlize)
library(Biobase)

theme_set(theme_pubr())

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 0.5),
  axis.title.y = element_text(size = 20),
  strip.text.x = element_text(size = 10, angle = 90))


#Main directory of sample enrichment figures
main = paste0("~/rprojects/results/imm_res_back/temp/first_test1/first_test/")

#List of traits
trait_list = c("ATD","CEL","IBD","JIA","MS","PBC","PSO","RA","T1D")


final_obj = list()
table_obj = list()
enrich_list = list()

for (i in c(1:length(trait_list))){
  plot_obj = c()
  print(i)
  
  #Directory of trait-specific sample enrichment results
  main_dir = paste0(main,trait_list[i],"/")
  
  cell_enrich = readRDS(paste0(main_dir,"cell_enrich.rds"))
  cell_list = readRDS(paste0(main_dir,"cell_list.rds"))
  #test_snp_number = readRDS(paste0(main_dir, "test_snp_number.rds"))
  #ref_snp_number = readRDS(paste0(main_dir, "ref_snp_number.rds"))
  
  
  real_pval = cell_enrich
  list = cell_list
  
  list = str_remove(list, ".bed")
  
  
  #motif_list = str_remove(motif_list, ".bed")
  list = str_remove(list, "EpiUwRmapDNase")
  
  
  cell_list = list
  
  
  gr1 = grep("CD",list)
  gr3 = grep("Thymus",list)
  
  
  gr = union(gr1,gr3)
  group = character(length(real_pval))
  
  group[gr] = "Immune"
  group[-gr] = "Non-immune"
  

  
  plot_frame = data.table(
    Group = group,
    Name = c(1:length(real_pval)),
    #Name = cell_list,
    Value = -log10(real_pval),
    key = "Group"
  )

  
  ColorsDT <-  data.table(Group= c("Immune","Non-immune"),
                          Color=c("#333BFF", "#CC6600"), key="Group")
  plot_frame[ColorsDT, Color := i.Color]
  
  
  if(trait_list[i] == "RA" | trait_list[i] == "JIA"){
    obj = ggplot(plot_frame ,aes(x= reorder(Name, Value), y=Value, fill=Group)) +
      geom_bar(stat="identity") + ylab("-log10(P value)") + xlab("Samples") + 
      ggtitle(trait_list[i])+
      ylim(c(0,8))+
      scale_x_discrete(guide = guide_axis(angle = 90))+
      geom_hline(yintercept=-log10(0.05/length(list)), linetype = "dashed")+
      theme(#axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  }else{
    obj = ggplot(plot_frame ,aes(x= reorder(Name, Value), y=Value, fill=Group)) +
      geom_bar(stat="identity") + ylab("-log10(P value)") + xlab("Samples") + 
      ggtitle(trait_list[i])+
      ylim(c(0,6))+
      scale_x_discrete(guide = guide_axis(angle = 90))+
      geom_hline(yintercept=-log10(0.05/length(list)), linetype = "dashed")+
      theme(#axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }
  final_obj[[i]] = obj + scale_fill_manual(values=c("orange", "gray"))

  enrich_list[[i]] = cell_enrich
  
  new_index = order(cell_enrich)
  cell_enrich_new = cell_enrich[new_index]
  cell_list_new = cell_list[new_index]
  
  table_frame = data.frame(matrix(nrow = 10, ncol = 2))
  colnames(table_frame) = c("sample_type", "enrichment")
  table_frame$sample_type = cell_list_new[1:10]
  table_frame$enrichment = format(round(-log10(cell_enrich_new[1:10]),2),nsmall = 2)
  
  
  table_obj[[i]] = ggtexttable(table_frame, rows = NULL, theme = ttheme("mOrange"),
                               cols = c("sample_ID","-log10(Pvalue of Enrichment)"))
  tab_add_title(table_obj[[i]], trait_list[i])
  
}


heat_matrix = matrix(0L ,nrow = length(trait_list), ncol = length(enrich_list[[1]]))

for (m in c(1:length(trait_list))){
  for (n in c(1:length(enrich_list[[1]]))){
    
    heat_matrix[m,n] = -log10(enrich_list[[m]][n])    
  
  }
  
}

rownames(heat_matrix) = trait_list
colnames(heat_matrix) = cell_list

data_plot = melt(heat_matrix)
names(data_plot) = c("trait","sample","enrichment")

#Generating Supp. Figure 1
ggplot_dir = paste0(main, "enrichment_pval1.tiff")
file.remove(ggplot_dir)

ggarrange(final_obj[[1]],final_obj[[2]],final_obj[[3]],final_obj[[4]], final_obj[[5]],
          final_obj[[6]], final_obj[[7]], final_obj[[8]], final_obj[[9]], common.legend = TRUE,
          labels = c("A","B","C","D","E","F","G","H","I"),
          nrow = 3,ncol = 3)

ggsave(ggplot_dir,dpi = 300,width = 10, height = 10, limitsize = FALSE)


#generating enrichment tables of samples
ggplot_dir = paste0(main, "enrichment_table.tiff")
file.remove(ggplot_dir)

p = ggarrange(table_obj[[3]],table_obj[[4]], table_obj[[5]],
          table_obj[[7]], table_obj[[8]], table_obj[[9]], 
          labels = c(trait_list[3],trait_list[4],trait_list[5],trait_list[7],
                     trait_list[8],trait_list[9]),
          nrow = 3,ncol = 2)

annotate_figure(p, top = text_grob(paste0("Enrichemnt Threshould is ",format(round(-log10(0.05/length(cell_list)), 2), nsmall = 2)), 
                                         color = "black", face = "bold", size = 14))
ggsave(ggplot_dir,dpi = 300,width = 15, height = 11, limitsize = FALSE)

#List of all samples
cell_type_dir = "~/rprojects/results/imm_res_back/temp/list.txt"
cell_type_list = read.table(cell_type_dir, sep = "\t", header = FALSE)
cell_type_list = cell_type_list$V1
cell_type_list = str_remove(cell_type_list, ".bed")
cell_type_list = str_remove(cell_type_list, "e*")

cell_type_index = list()
data_plot['cell_type'] = rep("other",length(data_plot$trait))
for (i in c(1:length(cell_type_list))){
  temp_cell = cell_type_list[i]
  
  index = grep(temp_cell, data_plot$sample)
  cell_type_index[[i]] = index
  data_plot[index,]$cell_type = rep(temp_cell,length(index))
}

cell_type_org = rep("others",length(cell_list))

for (i in c(1:length(cell_type_list))){
  temp_cell = cell_type_list[i]
  
  index = grep(temp_cell, cell_list)
  #cell_type_index[[i]] = index
  cell_type_org[index] = rep(temp_cell,length(index))
}

df = data.frame(cell_type = cell_type_org)
ha = HeatmapAnnotation(df = df)

immune_list1 = grep("CD",df$cell_type)
immune_list2 = grep("Thymus",df$cell_type)
immune_list = union(immune_list1,immune_list2)

nerv_list1 = grep("Brain", df$cell_type)
nerv_list2 = grep("Neuronal", df$cell_type)
nerv_list3 = grep("Spinal.Cord", df$cell_type)
nerv_list11 = union(nerv_list1,nerv_list2)
nerv_list = union(nerv_list11, nerv_list3)


heart_list = grep("Heart", df$cell_type)

intestine_list = grep("Intestine", df$cell_type)

kidney_list = grep("Kidney", df$cell_type)

lung_list = grep("Lung", df$cell_type)

muscle_list = grep("Muscle", df$cell_type)

stem_list = grep("Stem", df$cell_type)

annot_list = c()

annot_list = union(annot_list,immune_list)
annot_list = union(annot_list,nerv_list)
annot_list = union(annot_list,heart_list)
annot_list = union(annot_list,intestine_list)
annot_list = union(annot_list,kidney_list)
annot_list = union(annot_list,lung_list)
annot_list = union(annot_list,muscle_list)
annot_list = union(annot_list,stem_list)
all_list = c(1:length(df$cell_type))
others_list = all_list[- annot_list]


df["cell_type_new"] = rep(0,length(df$cell_type))
df[immune_list,]$cell_type_new = rep("Immune_cells",length(immune_list))
df[nerv_list,]$cell_type_new = rep("Nervous_system",length(nerv_list))
df[heart_list,]$cell_type_new = rep("Heart",length(heart_list))
df[intestine_list,]$cell_type_new = rep("Intestine",length(intestine_list))
df[kidney_list,]$cell_type_new = rep("Kidney",length(kidney_list))
df[lung_list,]$cell_type_new = rep("Lung",length(lung_list))
df[muscle_list,]$cell_type_new = rep("Muscle",length(muscle_list))
df[stem_list,]$cell_type_new = rep("Stem_cells",length(stem_list))
df[others_list,]$cell_type_new = rep("Other",length(others_list))


nonimmune_list = all_list[-immune_list]

df["cell_type_immune"] = rep(0,length(df$cell_type))

df[immune_list,]$cell_type_immune = df[immune_list,]$cell_type
df[nonimmune_list,]$cell_type_immune = rep("non_immune",length(nonimmune_list))

new_df = data.frame(matrix(nrow = length(df$cell_type),ncol = 2))
colnames(new_df) = c("Cell_type","Immune_type")

new_df$Cell_type = df$cell_type_new
new_df$Immune_type = df$cell_type_immune

new_index = order(new_df$Immune_type)
new_df = new_df[new_index,]
heat_matrix = heat_matrix[1:9,new_index]

new_index = order(new_df$Cell_type)
new_df = new_df[new_index,]
heat_matrix = heat_matrix[1:9,new_index]


col_fun = colorRamp2(c(0, max(heat_matrix)), c("white", "red"))

new_df1 = subset(new_df, select = -c(Cell_type,Immune_type))

ha = HeatmapAnnotation(df = new_df,
                       col = list(Cell_type = c("Immune_cells" = "#e41a1c", "Nervous_system" = "#377eb8",
                                "Heart" = "#4daf4a", "Intestine" = "#984ea3", "Kidney" = "#ff7f00",
                               "Lung" = "#ffff33", "Muscle" = "#a65628", "Stem_cells" = "#f781bf",
                               "Other" = "#999999","mamad" = "white"),
                               Immune_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                               "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                               "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                               "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                               "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                               "Mobilized.CD8" = "#fc4e2a")),
                        annotation_label = c("Cell Type","Immune Type")
                       
                       
                       )

#Generating heatmap A in Supp. Figure 2
colnames(heat_matrix) = c()
heat_dir = paste0(main, "new_heatmap.png")
file.remove(heat_dir)
png(heat_dir, width = 4000, height = 1800, res = 400)
draw(Heatmap(heat_matrix, top_annotation = ha, column_order = c(1:length(cell_list)),
        heatmap_legend_param = list(title = "Enrichment"), col = col_fun, 
        cell_fun = function(j,i,x,y,w,h,fill){
          if(heat_matrix[i,j] > -log10(0.05/length(list))){
            grid.text("*", x, y)
          }
        }))
dev.off()



immune_index = new_df$Cell_type == "Immune_cells"
new_df1 = new_df[immune_index,]
heat_matrix1 = heat_matrix[1:nrow(heat_matrix),immune_index]

new_df1 = subset(new_df1, select = -c(Cell_type))
ha = HeatmapAnnotation(df = new_df1,
                       col = list(Immune_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                                  "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                                  "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                                  "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                                  "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                                  "Mobilized.CD8" = "#fc4e2a")
                                  ),
                       annotation_label = "Immune Type"
)


#Generating heatmap B in Supp. Figure 2
colnames(heat_matrix1) = c()
heat_dir = paste0(main, "new_heatmap_immune.png")
file.remove(heat_dir)
png(heat_dir, width = 4000, height = 1800, res = 400)
draw(Heatmap(heat_matrix1, top_annotation = ha, column_order = c(1:sum(immune_index)),
             heatmap_legend_param = list(title = "Enrichment"), col = col_fun,
             cell_fun = function(j,i,x,y,w,h,fill){
               if(heat_matrix1[i,j] > -log10(0.05/length(list))){
                 grid.text("*", x, y)
               }}))
dev.off()

unique_new_df = unique(new_df)
new_heat_matrix = matrix(nrow = nrow(heat_matrix),ncol = length(unique_new_df$Cell_type))
rownames(new_heat_matrix) = rownames(heat_matrix)

for (i in c(1:length(unique_new_df$Cell_type))){
  found_index = (new_df$Cell_type == unique_new_df[i,]$Cell_type) &
    (new_df$Immune_type == unique_new_df[i,]$Immune_type)
  found_index = which(found_index)
  if (length(found_index) > 1){
    new_heat_matrix[1:nrow(new_heat_matrix),i] = rowMedians(heat_matrix[1:nrow(heat_matrix),found_index])
  }else{
    new_heat_matrix[1:nrow(new_heat_matrix),i] = heat_matrix[1:nrow(heat_matrix),found_index]
  }
}

col_fun = colorRamp2(c(0, max(new_heat_matrix)), c("white", "red"))


ha = HeatmapAnnotation(df = unique_new_df,
                       col = list(Cell_type = c("Immune_cells" = "#e41a1c", "Nervous_system" = "#377eb8",
                                                "Heart" = "#4daf4a", "Intestine" = "#984ea3", "Kidney" = "#ff7f00",
                                                "Lung" = "#ffff33", "Muscle" = "#a65628", "Stem_cells" = "#f781bf",
                                                "Other" = "#999999","mamad" = "white"),
                                  Immune_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                                  "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                                  "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                                  "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                                  "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                                  "Mobilized.CD8" = "#fc4e2a")),
                       annotation_label = c("Cell Type","Immune Type")
)


#Generating heatmap C in Supp. Figure 2
colnames(new_heat_matrix) = c()
heat_dir = paste0(main, "new_heatmap_unique.png")
file.remove(heat_dir)
png(heat_dir, width = 4000, height = 1800, res = 400)
draw(Heatmap(new_heat_matrix, top_annotation = ha, column_order = c(1:ncol(new_heat_matrix)),
             heatmap_legend_param = list(title = "Enrichment"), col = col_fun,
             cell_fun = function(j,i,x,y,w,h,fill){
               if(new_heat_matrix[i,j] > -log10(0.05/length(list))){
                 grid.text("*", x, y)
               }}))
dev.off()



immune_index = unique_new_df$Cell_type == "Immune_cells"
new_df1 = unique_new_df[immune_index,]
heat_matrix1 = new_heat_matrix[1:nrow(new_heat_matrix),immune_index]

new_df1 = subset(new_df1, select = -c(Cell_type))
ha = HeatmapAnnotation(df = new_df1,
                       col = list(Immune_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                                  "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                                  "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                                  "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                                  "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                                  "Mobilized.CD8" = "#fc4e2a")),
                       annotation_label = c("Immune Type")
)


#Generating heatmap D in Supp. Figure 2
colnames(heat_matrix1) = c()
heat_dir = paste0(main, "new_heatmap_immune_unique.png")
file.remove(heat_dir)

rownames(heat_matrix1) = rownames(new_heat_matrix)
png(heat_dir, width = 4000, height = 1800, res = 400)
draw(Heatmap(heat_matrix1, top_annotation = ha, column_order = c(1:sum(immune_index)),
             heatmap_legend_param = list(title = "Enrichment"), col = col_fun,
             cell_fun = function(j,i,x,y,w,h,fill){
               if(heat_matrix1[i,j] > -log10(0.05/length(list))){
                 grid.text("*", x, y)
               }}))
dev.off()

#Generating all qqplots
pict_dir = paste0(main, "all_qqplots.png")
file.remove(pict_dir)
png(pict_dir, height = 1800, width = 1800, res = 300)
par(mfrow=c(3,3))



for (i in c(1:length(trait_list))){
  plot_obj = c()
  print(i)
  
  #Directory of trait-specific sample enrichment results
  main_dir = paste0(main,trait_list[i],"/")
  
  cell_enrich = readRDS(paste0(main_dir,"cell_enrich.rds"))
  cell_list = readRDS(paste0(main_dir,"cell_list.rds"))
  #enrich_ratio = readRDS(paste0(main_dir,"enrich_ratio.rds"))
  #test_snp_number = readRDS(paste0(main_dir, "test_snp_number.rds"))
  #ref_snp_number = readRDS(paste0(main_dir, "ref_snp_number.rds"))
  
  #barplot(sort(cell_enrich))
  
  real_pval = cell_enrich
  list = cell_list
  
  list = str_remove(list, ".bed")
  
  
  #motif_list = str_remove(motif_list, ".bed")
  list = str_remove(list, "EpiUwRmapDNase")
  
  
  cell_list = list
  
  
  gr1 = grep("CD",list)
  #gr2 = grep("placenta",cell_list)
  gr3 = grep("Thymus",list)
  #gr4 = grep("Spleen",cell_list)
  
  
  #gr = union(union(gr1,gr2),union(gr3,gr4))
  gr = union(gr1,gr3)
  group = character(length(real_pval))
  
  group[gr] = "immune"
  group[-gr] = "non-immune"
  
  null_pval = runif(100000, min = 0, max = 1)
  
  
  
  immune_pval = real_pval[gr]
  nonimmune_pval = real_pval[-gr]
  
  lambda = median(-log10(real_pval))/-log10(0.5)
  qqPlot(y = -log10(real_pval), x = -log10(null_pval), equal.axes = TRUE, add.line = TRUE, qq.line.type = "0-1",
         xlab = "expected -log10(p_value)", ylab = "observed -log10(p_value)" , 
         main = paste0(trait_list[i],": lambda =",sprintf("%.3f",lambda)))
  
 

}

dev.off()



final_obj = list()
for (i in c(1:length(trait_list))){
  plot_obj = c()
  print(i)
  
  #Directory of trait-specific sample enrichment results
  main_dir = paste0(main,trait_list[i],"/")
  
  cell_enrich = readRDS(paste0(main_dir,"cell_enrich.rds"))
  cell_list = readRDS(paste0(main_dir,"cell_list.rds"))
  #enrich_ratio = readRDS(paste0(main_dir,"enrich_ratio.rds"))
  #test_snp_number = readRDS(paste0(main_dir, "test_snp_number.rds"))
  #ref_snp_number = readRDS(paste0(main_dir, "ref_snp_number.rds"))
  
  #barplot(sort(cell_enrich))
  
  real_pval = cell_enrich
  list = cell_list
  
  list = str_remove(list, ".bed")
  
  
  #motif_list = str_remove(motif_list, ".bed")
  list = str_remove(list, "EpiUwRmapDNase")
  
  
  cell_list = list
  
  found_list = c()
  
  Breast_gr = grep("Breast",cell_list)
  found_list = append(found_list,Breast_gr)
  
  CD14_gr = grep("CD14.Primary",cell_list)
  found_list = append(found_list,CD14_gr)
  
  CD19_gr = grep("CD19.Primary",cell_list)
  found_list = append(found_list,CD19_gr)
  
  CD20_gr = grep("CD20.Primary",cell_list)
  found_list = append(found_list,CD20_gr)
  
  CD34_gr = grep("CD34.Primary",cell_list)
  found_list = append(found_list,CD34_gr)
  
  mobilized_CD34_gr = grep("Mobilized.CD34.",cell_list)
  found_list = append(found_list,mobilized_CD34_gr)
  
  CD3_gr = grep("CD3.Primary",cell_list)
  found_list = append(found_list,CD3_gr)
  
  mobilized_CD3_gr = grep("Mobilized.CD3.",cell_list)
  found_list = append(found_list,mobilized_CD3_gr)
  
  cord_CD3_gr = grep("CD3.Cord",cell_list)
  found_list = append(found_list,cord_CD3_gr)
  
  CD4_gr = grep("CD4.Primary",cell_list)
  found_list = append(found_list,CD4_gr)
  
  mobilized_CD4_gr = grep("Mobilized.CD4.",cell_list)
  found_list = append(found_list,mobilized_CD4_gr)
  
  CD56_gr = grep("CD56.Primary",cell_list)
  found_list = append(found_list,CD56_gr)
  
  mobilized_CD56_gr = grep("Mobilized.CD56.",cell_list)
  found_list = append(found_list,mobilized_CD56_gr)
  
  CD8_gr = grep("CD8.Primary",cell_list)
  found_list = append(found_list,CD8_gr)
  
  mobilized_CD8_gr = grep("Mobilized.CD8.",cell_list)
  found_list = append(found_list,mobilized_CD8_gr)
  
  Adrenal_gr = grep("Adrenal",cell_list)
  found_list = append(found_list,Adrenal_gr)
  
  Brain_gr = grep("Brain",cell_list)
  found_list = append(found_list,Brain_gr)
  
  Heart_gr = grep("Heart",cell_list)
  found_list = append(found_list,Heart_gr)
  
  Intestine_gr = grep("Intestine",cell_list)
  found_list = append(found_list,Intestine_gr)
  
  Adrenal_gr = grep("Adrenal",cell_list)
  found_list = append(found_list,Adrenal_gr)
  
  Kidney_gr = grep("Kidney",cell_list)
  found_list = append(found_list,Kidney_gr)
  
  Lung_gr = grep("Lung",cell_list)
  found_list = append(found_list,Lung_gr)
  
  Muscle_gr = grep("Muscle",cell_list)
  found_list = append(found_list,Muscle_gr)
  
  
  Placenta_gr = grep("Placenta",cell_list)
  found_list = append(found_list,Placenta_gr)
  
  Thymus_gr = grep("Thymus",cell_list)
  found_list = append(found_list,Thymus_gr)
  
  Skin_gr = grep("Skin",cell_list)
  found_list = append(found_list,Skin_gr)
  
  BMP4_gr = grep("BMP4",cell_list)
  found_list = append(found_list,BMP4_gr)
  
  Stem_gr = grep("Stem",cell_list)
  found_list = append(found_list,Stem_gr)
  
  Neuronal_gr = grep("Neuronal",cell_list)
  found_list = append(found_list,Neuronal_gr)
  
  all_index = c(1:length(cell_list))
  others_gr = all_index[-found_list]
  
  
  
  box_frame = data.frame(matrix(ncol = 2, nrow = length(cell_list)))
  
  colnames(box_frame) = c("cell_type","pval")
  
  box_frame[Breast_gr,]$cell_type = rep("Breast",length(cell_list[Breast_gr]))
  box_frame[Breast_gr,]$pval = -log10(real_pval[Breast_gr])
  
  box_frame[Lung_gr,]$cell_type = rep("Lung",length(cell_list[Lung_gr]))
  box_frame[Lung_gr,]$pval = -log10(real_pval[Lung_gr])
  
  box_frame[Intestine_gr,]$cell_type = rep("Intestine",length(cell_list[Intestine_gr]))
  box_frame[Intestine_gr,]$pval = -log10(real_pval[Intestine_gr])
  
  box_frame[CD14_gr,]$cell_type = rep("CD14",length(cell_list[CD14_gr]))
  box_frame[CD14_gr,]$pval = -log10(real_pval[CD14_gr])
  
  box_frame[CD19_gr,]$cell_type = rep("CD19",length(cell_list[CD19_gr]))
  box_frame[CD19_gr,]$pval = -log10(real_pval[CD19_gr])
  
  box_frame[CD20_gr,]$cell_type = rep("CD20",length(cell_list[CD20_gr]))
  box_frame[CD20_gr,]$pval = -log10(real_pval[CD20_gr])
  
  box_frame[CD34_gr,]$cell_type = rep("CD34",length(cell_list[CD34_gr]))
  box_frame[CD34_gr,]$pval = -log10(real_pval[CD34_gr])
  
  box_frame[mobilized_CD34_gr,]$cell_type = rep("Mobilized.CD34",length(cell_list[mobilized_CD34_gr]))
  box_frame[mobilized_CD34_gr,]$pval = -log10(real_pval[mobilized_CD34_gr])
  
  
  box_frame[CD3_gr,]$cell_type = rep("CD3",length(cell_list[CD3_gr]))
  box_frame[CD3_gr,]$pval = -log10(real_pval[CD3_gr])
  
  box_frame[mobilized_CD3_gr,]$cell_type = rep("Mobilized.CD3",length(cell_list[mobilized_CD3_gr]))
  box_frame[mobilized_CD3_gr,]$pval = -log10(real_pval[mobilized_CD3_gr])
  
  
  box_frame[cord_CD3_gr,]$cell_type = rep("Cord.CD3",length(cell_list[cord_CD3_gr]))
  box_frame[cord_CD3_gr,]$pval = -log10(real_pval[cord_CD3_gr])
  
  
  box_frame[CD4_gr,]$cell_type = rep("CD4",length(cell_list[CD4_gr]))
  box_frame[CD4_gr,]$pval = -log10(real_pval[CD4_gr])
  
  box_frame[mobilized_CD4_gr,]$cell_type = rep("Mobilized.CD4",length(cell_list[mobilized_CD4_gr]))
  box_frame[mobilized_CD4_gr,]$pval = -log10(real_pval[mobilized_CD4_gr])
  
  
  box_frame[CD56_gr,]$cell_type = rep("CD56",length(cell_list[CD56_gr]))
  box_frame[CD56_gr,]$pval = -log10(real_pval[CD56_gr])
  
  box_frame[mobilized_CD56_gr,]$cell_type = rep("Mobilized.CD56",length(cell_list[mobilized_CD56_gr]))
  box_frame[mobilized_CD56_gr,]$pval = -log10(real_pval[mobilized_CD56_gr])
  
  box_frame[CD8_gr,]$cell_type = rep("CD8",length(cell_list[CD8_gr]))
  box_frame[CD8_gr,]$pval = -log10(real_pval[CD8_gr])
  
  
  box_frame[mobilized_CD8_gr,]$cell_type = rep("Mobilized.CD8",length(cell_list[mobilized_CD8_gr]))
  box_frame[mobilized_CD8_gr,]$pval = -log10(real_pval[mobilized_CD8_gr])
  
  box_frame[Adrenal_gr,]$cell_type = rep("Adrenal",length(cell_list[Adrenal_gr]))
  box_frame[Adrenal_gr,]$pval = -log10(real_pval[Adrenal_gr])
  
  box_frame[Kidney_gr,]$cell_type = rep("Kidney",length(cell_list[Kidney_gr]))
  box_frame[Kidney_gr,]$pval = -log10(real_pval[Kidney_gr])
  
  box_frame[Brain_gr,]$cell_type = rep("Brain",length(cell_list[Brain_gr]))
  box_frame[Brain_gr,]$pval = -log10(real_pval[Brain_gr])
  
  
  box_frame[Heart_gr,]$cell_type = rep("Heart",length(cell_list[Heart_gr]))
  box_frame[Heart_gr,]$pval = -log10(real_pval[Heart_gr])
  
  
  box_frame[Muscle_gr,]$cell_type = rep("Muscle",length(cell_list[Muscle_gr]))
  box_frame[Muscle_gr,]$pval = -log10(real_pval[Muscle_gr])
  
  box_frame[Placenta_gr,]$cell_type = rep("Placenta",length(cell_list[Placenta_gr]))
  box_frame[Placenta_gr,]$pval = -log10(real_pval[Placenta_gr])
  
  box_frame[Stem_gr,]$cell_type = rep("Stem",length(cell_list[Stem_gr]))
  box_frame[Stem_gr,]$pval = -log10(real_pval[Stem_gr])
  
  
  box_frame[Skin_gr,]$cell_type = rep("Skin",length(cell_list[Skin_gr]))
  box_frame[Skin_gr,]$pval = -log10(real_pval[Skin_gr])
  
  box_frame[Neuronal_gr,]$cell_type = rep("Neuronal",length(cell_list[Neuronal_gr]))
  box_frame[Neuronal_gr,]$pval = -log10(real_pval[Neuronal_gr])
  
  
  box_frame[Thymus_gr,]$cell_type = rep("Thymus",length(cell_list[Thymus_gr]))
  box_frame[Thymus_gr,]$pval = -log10(real_pval[Thymus_gr])
  
  box_frame[others_gr,]$cell_type = rep("Others",length(cell_list[others_gr]))
  box_frame[others_gr,]$pval = -log10(real_pval[others_gr])
  
  
  box_frame = box_frame[complete.cases(box_frame),]
  
  ggplot_dir = paste0(main_dir, "enrichment_boxplot.tiff")
  file.remove(ggplot_dir)
  
  box_frame["Immune_type"] = rep("Non_immune",length(box_frame$cell_type))
  
  cd_index = grep("CD",box_frame$cell_type)
  thymus_index = grep("Thymus", box_frame$cell_type)
  
  immune_index = union(cd_index,thymus_index)
  
  box_frame[immune_index,]$Immune_type = "Immune"
  
  box_frame["Cell Type"] = box_frame$Immune_type
  if(trait_list[i] == "RA" | trait_list[i] == "JIA"){
  final_obj[[i]] = ggplot(box_frame, aes(x=reorder(cell_type,pval,FUN = median), y=pval, fill = `Cell Type`)) +
    geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90))+
    xlab("Cell Types")+
    ylim(c(0,8))+
    ggtitle(trait_list[[i]])+
    geom_hline(yintercept=-log10(0.05/length(list)), linetype = "dashed", color = "gray")+
    geom_hline(yintercept= -log10(0.05), linetype = "dotted", color = "gray")+
    #guides(color=FALSE)+
    ylab("-log10(P value)") + scale_fill_manual(values=c("orange", "gray"))
  }else{
    final_obj[[i]] = ggplot(box_frame, aes(x=reorder(cell_type,pval,FUN = median), y=pval, fill = `Cell Type`)) +
      geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90))+
      xlab("Cell Types")+
      ylim(c(0,6))+
      ggtitle(trait_list[[i]])+
      geom_hline(yintercept=-log10(0.05/length(list)), linetype = "dashed", color = "gray")+
      geom_hline(yintercept= -log10(0.05), linetype = "dotted", color = "gray")+
      #guides(color=FALSE)+
      ylab("-log10(P value)") + scale_fill_manual(values=c("orange", "gray"))
  }
}

#Generating Figure 1
ggplot_dir = paste0(main, "enrichment_boxplot.tiff")
file.remove(ggplot_dir)

ggarrange(final_obj[[1]],final_obj[[2]],final_obj[[3]],final_obj[[4]], final_obj[[5]],
          final_obj[[6]], final_obj[[7]], final_obj[[8]], final_obj[[9]],common.legend = TRUE,
          labels = c("A","B","C","D","E","F","G","H","I"),
          nrow = 3,ncol = 3)

ggsave(ggplot_dir,dpi = 300,width = 16, height = 10, limitsize = FALSE)
