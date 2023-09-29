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
library("scales")
library(patchwork)
library(reshape2)
library(circlize)
library(cowplot)

theme_set(theme_pubr())

My_Theme = theme(
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  strip.text.x = element_text(size = 30, angle = 90),
  title = element_text(size = 40),
  legend.text = element_text(size = 25),
  legend.title = element_text(size = 30))

#List of motifs
motif_tf_dir = "~/rprojects/motif_list/motif_name.txt"
motif_tf_data = read.table(motif_tf_dir, header = TRUE)

#Main directory of trait-specific results
main_all = paste0("~/rprojects/results/imm_res_back/temp/new_generated_res/")

#List of traits
trait_list = c("IBD","JIA","MS","PSO","RA","T1D","CEL")

final_obj = list()
final_data_obj = list()

final_matrix_alltraits = list()
final_celltype_alltraits = list()
final_trait_alltraits = list()
final_trait_number = c()

for (i in c(1:length(trait_list))){

  #Directory of each trait
  main_dir = paste0(main_all,trait_list[i],"/")
  
  cell_type_list = c("CD3.Primary","CD8.Primary","CD56.Primary","Thymus","CD14.Primary",
                     "CD4.Primary","CD34.Primary","CD19.Primary","CD20.Primary")
  cell_type_list1 = str_remove(cell_type_list, ".Primary")
  final_data = list()
  
  final_matrix = c()
  final_celltype = c()
  
  for (j in c(1:length(cell_type_list))){
    
    #ChromHMM result of each cell type for a trait
    main = paste0(main_dir, cell_type_list[j],"/",cell_type_list1[j],".txt")
    
    
    if (file.exists(main)){
      data = read.table(main, sep = "\t", header = TRUE)
      data = subset(data,select = -Genome..)
      data_raw = data
      
      heat_matrix = matrix(0L ,nrow = nrow(data_raw), ncol = ncol(data_raw) - 1)
      
      for (m in c(1:nrow(data_raw))){
        for (n in c(1:ncol(data_raw)-1)){
          
          if (data_raw[m,n+1] > 2){
            heat_matrix[m,n] = log(as.numeric(data_raw[m,n+1])+1)
          }else{
            heat_matrix[m,n] = 0
          }
          
        }
        
      }
      motif_list = colnames(data_raw)[2:ncol(data_raw)]
      motif_list = str_remove(motif_list,".bed")
      show_list = c()
      for (s in c(1:length(motif_list))){
        motif_temp = motif_list[s]
        index = motif_tf_data$Motif_ID == motif_temp
        tf_temp = motif_tf_data[index,]$Motif_name
        tf_temp = tf_temp[1]
        
        show_list = append(show_list,paste0(tf_temp," (",motif_temp,")"))
        
      }
      show_list = str_remove(show_list,"^/")
      colnames(heat_matrix) = show_list
      rownames(heat_matrix) = data_raw$State..User.order.
      data = melt(data)
      data['cell_type'] = rep(cell_type_list[j],length(data$value))
      final_data[[j]] = data
      final_matrix = cbind(final_matrix,heat_matrix)
      final_celltype = append(final_celltype,rep(cell_type_list[j],ncol(data_raw)-1))
    }
    
  }
  
  final_matrix_alltraits[[i]] = final_matrix
  final_celltype_alltraits[[i]] = final_celltype
  final_trait_alltraits[[i]] = rep(trait_list[i],length(final_celltype))
  final_trait_number = append(final_trait_number,length(final_celltype))
  data_plot = do.call(rbind,final_data)
  data_plot["enrichment"] = data_plot$value
  
  data_plot$variable = str_remove(as.character(data_plot$variable),".bed")
  
  state_string = str_remove(data_plot$State..User.order.,"[0-9]*")
  state_string = str_remove(state_string,"_*")
  data_plot$State..User.order. = state_string
  
  data_plot$enrichment = log10(data_plot$enrichment + 1)
  
  new_motif = c()
  for (j in c(1:length(data_plot$variable))){
    motif_temp = data_plot[j,]$variable
    index = motif_tf_data$Motif_ID == motif_temp
    tf_temp = motif_tf_data[index,]$Motif_name
    tf_temp = tf_temp[1]
    
    new_motif = append(new_motif,paste0(tf_temp," (",motif_temp,")"))
    
  }
  
  data_plot$variable = new_motif
  
  data_plot['trait']  = rep(trait_list[i],length(data_plot$value))
  
  final_data_obj[[i]] = data_plot
  
  final_obj[[i]] = ggplot(data_plot, aes(variable, State..User.order., fill= enrichment)) + 
    geom_tile()+ 
    scale_fill_gradient(low = "#c6dbef",high = "#08306b",guide = "colorbar",name = "log10(enrichment_ratio)")+
    xlab("TF(Motif)")+ylab("Functional Region")+ggtitle(paste0("Enrichmnet ratio, ",trait_list[[i]]))+
    scale_x_discrete(guide = guide_axis(angle = 90))+
    facet_grid(.~ cell_type, space = 'free_x', scales = 'free_x') +
    theme(legend.title = "enrichment ratio")+
    theme_classic(base_size = 14, base_family = 'mono') +
    theme(panel.grid.minor.x = element_blank()) + 
    # remove facet spacing on x-direction
    theme(panel.spacing.x = unit(0,"line")) +
    # switch the facet strip label to outside 
    # remove background color
    theme(strip.placement = 'outside',
          strip.background.x = element_blank()) + My_Theme
  

}


max_val = max(do.call(cbind,final_matrix_alltraits))

heat_obj = list()

#Main directory of final figure
main_dir_heat = "~/rprojects/results/imm_res_back/temp/new_generated_res/"
for (i in c(1:length(final_matrix_alltraits))){
  
  final_small_matrix = list()
  final_small_celltype = list()
  final_small_trait = list()
  
  final_small_matrix[[1]] = final_matrix_alltraits[[i]]
  final_small_celltype[[1]] = final_celltype_alltraits[[i]]
  final_small_trait[[1]] = final_trait_alltraits[[i]]

  
  
  final_plot_matrix = do.call(cbind,final_small_matrix)

  final_plot_celltype = c()
  for (j in c(1:length(final_small_celltype))){
    final_plot_celltype = append(final_plot_celltype,final_small_celltype[[j]])
  }
  
  final_plot_trait = c()
  for (j in c(1:length(final_small_trait))){
    final_plot_trait = append(final_plot_trait,final_small_trait[[j]])
  }
  
  
  
  new_df = data.frame(matrix(nrow = length(final_plot_celltype),ncol = 1))
  colnames(new_df) = c("Cell_type")
  
  new_df$Cell_type = final_plot_celltype
  
  col_fun = colorRamp2(c(0,  max_val), c("#c6dbef", "#08306b"))
  
  
  ha = HeatmapAnnotation(df = new_df,
                         col = list(Cell_type = c("CD14.Primary" = "#a6cee3", "CD19.Primary" = "#1f78b4", "CD20.Primary" = "#b2df8a",
                                                  "CD3.Primary" = "#33a02c", "CD34.Primary" = "#fb9a99", "CD4.Primary" = "#e31a1c",
                                                  "CD56.Primary" = "#fdbf6f", "CD8.Primary" = "#ff7f00", "Thymus" = "#cab2d6",
                                                  "non_immune" = "white", "CD3.Cord" = "#00441b", "Mobilized.CD3" = "#99d8c9",
                                                  "Mobilized.CD34" = "#e7298a", "Mobilized.CD4" = "#b30000", "Mobilized.CD56" = "#ffffcc", 
                                                  "Mobilized.CD8" = "#fc4e2a")),
                                    annotation_name_gp = gpar(fontface = "bold"),
                         annotation_label = "Cell Type",
                         show_legend = FALSE
                         
                         
  )
  
  
  heat_obj[[i]] = Heatmap(final_plot_matrix, top_annotation = ha, column_order = c(1:ncol(final_plot_matrix)),
                                row_order = c(1:nrow(final_plot_matrix)),
                                heatmap_legend_param = list(title = "Enrichment"), col = col_fun,
                                column_split = final_plot_trait,
                                show_heatmap_legend = FALSE,column_title_gp = grid::gpar(fontface = "bold")#,width = final_trait_number[i])
  )
  heat_file  = paste0(main_dir_heat,trait_list[i],"/new_heatmap_chrom.png")
  file.remove(heat_file)
  png(heat_file, width = 3600*(final_trait_number[i]/52) + (6/52)*3600, height = 1800, res = 300)
  draw(heat_obj[[i]])
  
  dev.off()
  
}

lgd1 = Legend(col_fun = col_fun, title = "log10(fold-enrichment)",direction = "horizontal",
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
              title = "Cell Type", direction = "vertical")

pd = packLegend(lgd1, lgd2, direction = "horizontal")

pd_dir = paste0(main_dir_heat,"annotation.png")
file.remove(pd_dir)
png(pd_dir, width = 3600*(53/53), height = 1600, res = 400)
draw(pd)
dev.off()

#Generating Figure 5
heat_dir = paste0(main_all, "chrom_flag.png")
file.remove(heat_dir)
png(heat_dir, width = 10000, height = 4800, res = 400)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2,widths = c(53,10))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat_obj[[2]]+heat_obj[[3]]+heat_obj[[5]]+heat_obj[[6]],
     auto_adjust = FALSE,newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(heat_obj[[1]]+heat_obj[[4]]+heat_obj[[7]],auto_adjust = FALSE, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(lgd2)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(lgd1)
upViewport()

dev.off()


