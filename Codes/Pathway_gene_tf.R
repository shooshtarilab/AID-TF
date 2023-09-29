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
library(ape)
library(dplyr)
library(stringr)
library(stringi)
library(xlsx)
theme_set(theme_pubr())

#The parent directory of trait files each having affected binding sites 
#as a .bed file
data_dir = "~/rprojects_whole/MS_postprocess/results/pathways/"

#List of traits
trait_list = c("IBD","JIA","MS","PSO","RA","T1D","CEL")

for (trait in trait_list){
  region_dir = paste0(data_dir,trait, "/",trait, "_regions.bed")
  region_data = read.table(file = region_dir, header = FALSE, sep = "\t")
  job = submitGreatJob(region_data)
  tb = getEnrichmentTables(job)
  BiologProcess = tb[[2]]
  BiologProcess = BiologProcess[BiologProcess$Binom_Adjp_BH < 0.05,]
  if (nrow(BiologProcess) > 100){
    BiologProcess = BiologProcess[1:100,]
  }
  BiologProcess[["gene_list"]] = rep("none", nrow(BiologProcess))
  for (go_id in BiologProcess$ID){
    print("gene_started")
    #final_dat = plotRegionGeneAssociationGraphs(job, which_plot = 1:3, ontology = "GO Biological Process",
    #                                            term_id = go_id, request_interval = 10, max_tries = 100, verbose = TRUE)
    final_dat = getRegionGeneAssociations(job, term_id = go_id,
                                          ontology = "GO Biological Process")
    print("gene_ended")
    gene_list = unique(unlist(final_dat$annotated_genes))
    gene_list = knitr::combine_words(gene_list, and = ",")
    BiologProcess[BiologProcess$ID == go_id,]$gene_list = gene_list
    print(go_id)
  }
  
  #Writing AID-related pathways as a table in the folder of that AID
  BiologProcess = select(BiologProcess, c("name", "Binom_Adjp_BH", "gene_list"))
  final_frame_dir = paste0(data_dir, trait, "/", trait, "_final_frame.txt")
  write.table(BiologProcess, file = final_frame_dir, sep = "\t", col.names = TRUE, 
              row.names = FALSE, quote = FALSE)
  print(trait)
}


#The address of the final dataframe comprising all results
main_datframe_dir = "~/rprojects_whole/MS_postprocess/data/final_dataframe.txt"
main_dataframe = read.table(file = main_datframe_dir, header = TRUE, sep = "\t")

trait_list = unique(main_dataframe$trait)


#The main directory of TF-gene pair results 
main_dir = "~/rprojects_whole/MS_postprocess/results/final_results/all_TF_gene/"

for (trait in trait_list){
  trait_dataframe = main_dataframe[main_dataframe$trait == trait,]

  trait_data = select(trait_dataframe, c("chr","start","end","cell_type","TF"))
  trait_data = unique(trait_data)
  trait_data = trait_data[order(as.integer(sub("chr","",trait_data$chr))),]

  trait_data_region = unique(select(trait_data, c("chr","start","end")))
  job = submitGreatJob(trait_data_region)
  tb = getEnrichmentTables(job)
  
  final_dat = getRegionGeneAssociations(job)
  final_trait_dataframe = list()
  #final_dat = plotRegionGeneAssociations(job, which_plot = 1:3,
  #                                            request_interval = 10, max_tries = 100, verbose = TRUE)
  for (i in c(1:nrow(trait_data))){
    #found_index = ((seqnames(final_dat) == trait_data$chr) &
    #                 (start(final_dat) == trait_data$start) &
    #                 (end(final_dat) == trait_data$end))
    
    trait_data_granges = makeGRangesFromDataFrame(trait_data[i,])
    found_index = unique(subjectHits(findOverlaps(trait_data_granges, final_dat)))
    frame_length = length(unlist(final_dat$annotated_genes[found_index]))
    temp_frame = as.data.frame(matrix(0, nrow = frame_length, ncol = 7))
    colnames(temp_frame) = c("chr","start","end","cell_type","TF","gene","distTSS")
    temp_frame$chr = rep(trait_data$chr[i], nrow(temp_frame))
    temp_frame$start = rep(trait_data$start[i], nrow(temp_frame))
    temp_frame$end = rep(trait_data$end[i], nrow(temp_frame))
    temp_frame$cell_type = rep(trait_data$cell_type[i], nrow(temp_frame))
    temp_frame$TF = rep(trait_data$TF[i], nrow(temp_frame))
    temp_frame$gene = unlist(final_dat$annotated_genes[found_index])
    temp_frame$distTSS = unlist(final_dat$dist_to_TSS[found_index])
    
    final_trait_dataframe[[i]] = unique(temp_frame)
  }
  final_trait_dataframe = do.call(rbind, final_trait_dataframe)
  final_trait_dataframe["trait"] = trait
  final_trait_dataframe = final_trait_dataframe[order(as.integer(sub("chr","",final_trait_dataframe$chr)),
                                                      final_trait_dataframe$start),]
  #Writing TF-gene results for each trait
  save_dir = paste0(main_dir, trait, ".txt")
  write.table(final_trait_dataframe, file = save_dir, sep = "\t", col.names = TRUE, 
              row.names = FALSE, quote = FALSE)
  print(trait)
  
}


#List of traits
trait_list = c("IBD","JIA","MS","PSO","RA","T1D")

#The directory of AID-related pathways' tables
main_dir = "~/rprojects_whole/MS_postprocess/results/pathways/"
#The directory of trait-specific TF-gene pairs
main_dir_gene = "~/rprojects_whole/MS_postprocess/results/final_results/all_TF_gene/"
final_frame = list()
final_frame_gene = list()
count = 1

for (trait in trait_list){
  temp_file = paste0(main_dir, "/", trait, "/", trait, "_final_frame.txt")
  temp_data = read.table(file = temp_file, header = TRUE, sep = "\t")
  
  if (nrow(temp_data) > 100){
    temp_data = temp_data[1:100,]
  }else{
    temp_data = temp_data
  }
  temp_gene_list = str_split(temp_data$gene_list, pattern = ",",
                             simplify = FALSE)
  temp_gene_list = stri_remove_empty(unique(str_trim(unlist(temp_gene_list))))
  temp_data[["trait"]] = rep(trait, nrow(temp_data))
  temp_data$gene_list = stri_remove_empty(str_trim(temp_data$gene_list))
  final_frame[[count]] = temp_data
  
  temp_gene_file = paste0(main_dir_gene, "/", trait, ".txt") 
  temp_gene = read.table(file = temp_gene_file, header = TRUE, sep = "\t")  
  temp_gene = temp_gene[temp_gene$gene %in% temp_gene_list,]
  final_frame_gene[[count]] = temp_gene
  
  count = count + 1
}

final_frame = do.call(rbind, final_frame)
final_frame_gene = do.call(rbind, final_frame_gene)

gene_tf = unique(paste0(final_frame_gene$TF,
                        ",", final_frame_gene$gene))

final_frame_gene[["gene_tf"]] = paste0(final_frame_gene$TF,
                                       ",", final_frame_gene$gene)

final_frame_gene_nocell = subset(final_frame_gene, select = -c(cell_type))
final_frame_gene_nocell = unique(final_frame_gene_nocell)

final_frame_gene_nocell[["cell_types"]] = NA

for (temp_dat in gene_tf){
  cell_list = unique(final_frame_gene$cell_type[final_frame_gene$gene_tf == temp_dat])
  cell_list = paste(cell_list, collapse = ", ")
  
  found_index = final_frame_gene_nocell$gene_tf == temp_dat
  final_frame_gene_nocell$cell_types[found_index] = cell_list
}

final_frame_gene_nocell = subset(final_frame_gene_nocell,
                                 select = -c(gene_tf))

#The output directory
main_dir = "~/rprojects_whole/MS_postprocess/results/final_results/"
#The address of the table having TF-gene results for all traits
gene_dir = paste0(main_dir, "gene_tf_table.xlsx")
#The address of the table showing affected pathways and 
#genes involved in them for all traits
pathway_dir = paste0(main_dir, "pathway_table.xlsx")

write.xlsx(final_frame, file = pathway_dir, row.names = FALSE)
write.xlsx(final_frame_gene_nocell, file = gene_dir, row.names = FALSE)

