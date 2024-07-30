library(GenomicRanges)
library(stringr)
library(rtracklayer)
library(dplyr)
library(liftOver)


#mm10 to hg19 chain file address
chain_file = "~/rprojects/review_data/Codes/mm10ToHg19.over.chain"
chain_file_hg38 = "~/rprojects/review_data/Codes/hg38ToHg19.over.chain"

grg_liftover_mm10TOhg19 = function(in_grg,chain_in){
  chain_data = import.chain(chain_in)
  seqlevelsStyle(in_grg) = "UCSC" 
  out_grg = liftOver(in_grg, chain_data)
  out_grg = unlist(out_grg)
  return(out_grg)
}


#Directory of the dataframe comprising all results
full_table_dir = "~/rprojects/review_data/results/new_generated_res_test/final_dataframe_humna.txt"

full_table = read.table(file = full_table_dir, sep = "\t", header = TRUE)

#Directory of the ChIP-seq data
chip_dir = "~/rprojects/review_data/ChIP-seq/"

#Getting the list of TFs in the ChIP-seq data
chip_file_list = list.files(chip_dir)

#Getting the list of TFs with ChIP-seq data available
chip_tfs <- str_extract(chip_file_list, "(?<=^Oth\\.Bld\\.\\d{2}\\.)[^\\.]+(?=\\.AllCell\\.bed$)")

#Creating a dataframe from TF names and their ChIP file
chip_frame = as.data.frame(matrix(0, nrow = length(chip_tfs),
                                  ncol = 2))

colnames(chip_frame) = c("TF","ChIP_file")

chip_frame$TF = chip_tfs
chip_frame$ChIP_file = chip_file_list

#Loading the ChIP-seq data of all TFs
chip_data_list = list()

for (tf in chip_frame$TF){
  tf_chip_temp_dir = paste0(chip_dir, chip_frame$ChIP_file[chip_frame$TF == tf])
  tf_chip_temp = import.bed(tf_chip_temp_dir)
  chip_data_list[[tf]] = GRanges(tf_chip_temp)
  print(tf)
}

#Directory of the immune ChIP-seq data
chip_dir_immune = "~/rprojects/review_data/ChIP_immue/"

#Getting the list of TFs in the ChIP-seq data
chip_file_list_immune = list.files(chip_dir_immune)

#Getting the list of TFs with ChIP-seq data available
chip_tfs_immune <- str_split(chip_file_list_immune, pattern = ".csv",
                             simplify = TRUE)[,1]

#Creating a dataframe from TF names and their ChIP file
chip_frame_immune = as.data.frame(matrix(0, nrow = length(chip_tfs_immune),
                                         ncol = 2))

colnames(chip_frame_immune) = c("TF","ChIP_file")

chip_frame_immune$TF = chip_tfs_immune
chip_frame_immune$ChIP_file = chip_file_list_immune

#Loading the ChIP-seq data of all TFs
chip_data_list_immune = list()

for (tf in chip_frame_immune$TF){
  tf_chip_temp_dir = paste0(chip_dir_immune, chip_frame_immune$ChIP_file[chip_frame_immune$TF == tf])
  tf_chip_temp = read.table(file = tf_chip_temp_dir, 
                            sep = ",", 
                            header = TRUE)
  tf_chip_temp = GRanges(ranges = IRanges(start = tf_chip_temp$OCR_start,
                                          end = tf_chip_temp$OCR_end),
                         seqnames = tf_chip_temp$Chr)
  temp_grg = GRanges(tf_chip_temp)
  temp_grg = grg_liftover_mm10TOhg19(in_grg = temp_grg, chain_in = chain_file)
  chip_data_list_immune[[tf]] = temp_grg
  print(tf)
}

#Directory of the ChIP-seq data remap
chip_dir_remap = "~/rprojects/review_data/remap_chip_human/"

#Getting the list of TFs in the ChIP-seq data
chip_file_list_remap = list.files(chip_dir_remap)

#Getting the list of TFs with ChIP-seq data available
chip_tfs_remap <- gsub("remap2022_|_all_macs2_hg38_v1_0.bed", "", chip_file_list_remap)

#Creating a dataframe from TF names and their ChIP file
chip_frame_remap = as.data.frame(matrix(0, nrow = length(chip_tfs_remap),
                                  ncol = 2))

colnames(chip_frame_remap) = c("TF","ChIP_file")

chip_frame_remap$TF = chip_tfs_remap
chip_frame_remap$ChIP_file = chip_file_list_remap

#Loading the ChIP-seq data of all TFs
chip_data_list_remap = list()

for (tf in chip_frame_remap$TF){
  tf_chip_temp_dir = paste0(chip_dir_remap, chip_frame_remap$ChIP_file[chip_frame_remap$TF == tf])
  tf_chip_temp = import.bed(tf_chip_temp_dir)
  tf_chip_temp = grg_liftover_mm10TOhg19(tf_chip_temp, 
                                         chain_in = chain_file_hg38)
  chip_data_list_remap[[tf]] = GRanges(tf_chip_temp)
  print(tf)
}



#Getting the list of TFs with at least one ChIP data
all_chip_tf_list = unique(c(chip_frame$TF, chip_frame_immune$TF,
                            chip_frame_remap$TF))

chip_data_all = list()

for (tf in all_chip_tf_list){
  tf_data = c()
  if (tf %in% names(chip_data_list)){
    temp_data = chip_data_list[[tf]]
    mcols(temp_data) = NULL
    tf_data = c(tf_data, temp_data)
    #mcols(tf_data) = NULL
  }
  if(tf %in% names(chip_data_list_immune)){
    temp_data = chip_data_list_immune[[tf]]
    mcols(temp_data) = NULL
    tf_data = c(tf_data, temp_data)
    #mcols(tf_data) = NULL
  }
  if(tf %in% names(chip_data_list_remap)){
    temp_data = chip_data_list_remap[[tf]]
    mcols(temp_data) = NULL
    tf_data = c(tf_data, temp_data)
    #mcols(tf_data) = NULL
  }
  tf_data = do.call(c, tf_data)
  chip_data_all[[tf]] = tf_data
}

# for(tf in all_chip_tf_list){
#   if (tf %in% chip_frame$TF){
#     mcols(chip_data_list[[tf]]) = NULL
#     if (tf %in% chip_frame_immune$TF){
#       temp_grg = c(chip_data_list[[tf]], chip_data_list_immune[[tf]])
#     }
#     else{
#       temp_grg = chip_data_list[[tf]]
#     }
#   }else{
#     temp_grg = chip_data_list_immune[[tf]]
#   }
#   chip_data_all[[tf]] = temp_grg
# }

#Adding ChIP and conservation columns to the main table
full_table[["Chip_conf"]] = "No"
full_table[["Mice_Conserved"]] = "No"

#Confirming TF binding sites using ChIP-seq data
for (tf in names(chip_data_all)){
  tf_temp_data = full_table[full_table$TF == tf,]
  chip_tf_data = chip_data_all[[tf]]
  
  tf_temp_grg = GRanges(ranges = IRanges(
    start = tf_temp_data$start,
    end = tf_temp_data$end
  ),seqnames = tf_temp_data$chr)
  
  tf_conf_ind = unique(queryHits(findOverlaps(tf_temp_grg,
                                              chip_tf_data)))
  
  tf_temp_data$Chip_conf[tf_conf_ind] = 'YES'
  full_table[full_table$TF == tf,] = tf_temp_data
}

chip_conf_full_table = full_table[full_table$Chip_conf == 'YES',]

tf_conf_data = as.data.frame(matrix(0, nrow = length(names(chip_data_all)),
                                    ncol = 4))

colnames(tf_conf_data) = c("TF","Bind_num","ChiP_bind_num",
                           "ChIP_conf_ratio")
tf_conf_data$TF = names(chip_data_all)

for (tf in names(chip_data_all)){
  tf_temp_data = full_table[full_table$TF == tf,]
  chip_conf_tf_temp_data = chip_conf_full_table[chip_conf_full_table$TF == tf,]
  #conf_perc = nrow(chip_conf_tf_temp_data)/nrow(tf_temp_data)
  #tf_conf_data[tf_conf_data$TF == tf,]$Bind_num = nrow(tf_temp_data)
  #tf_conf_data[tf_conf_data$TF == tf,]$ChiP_bind_num = nrow(chip_conf_tf_temp_data)
  #tf_conf_data[tf_conf_data$TF == tf,]$ChIP_conf_perc = conf_perc
  chip_conf_tf_temp_data = unique(chip_conf_tf_temp_data[,c("chr","start","end")])
  tf_temp_data = unique(tf_temp_data[,c("chr","start","end")])
  conf_perc = nrow(chip_conf_tf_temp_data)/nrow(tf_temp_data)
  tf_conf_data[tf_conf_data$TF == tf,]$Bind_num = nrow(tf_temp_data)
  tf_conf_data[tf_conf_data$TF == tf,]$ChiP_bind_num = nrow(chip_conf_tf_temp_data)
  tf_conf_data[tf_conf_data$TF == tf,]$ChIP_conf_ratio = conf_perc
}

#The output table of ChIP-seq analysis
tf_conf_file = "~/rprojects/review_data/results/new_generated_res_test/aff_tf_conf.txt"
write.table(file = tf_conf_file, tf_conf_data,sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

#Conservation file
cons_file = "~/rprojects/review_data/Conserved_hg19/hg19.LECIFv1.1.bw"
cons_data = import.bw(cons_file)

cons_data_filt = cons_data[cons_data$score>0.2]


full_grg = GRanges(ranges = IRanges(start = full_table$start,
                                    end = full_table$end),
                   seqnames = full_table$chr)

df_cons_idx = unique(queryHits(findOverlaps(full_grg,
                                            cons_data_filt)))

full_table$Mice_Conserved[df_cons_idx] = "YES"

conserv_conf_full_table = full_table[full_table$Mice_Conserved == 'YES',]

tf_conserv_data = as.data.frame(matrix(0, nrow = length(unique(full_table$TF)),
                                    ncol = 4))

colnames(tf_conserv_data) = c("TF","Bind_num","Conserved_bind_num",
                           "Conserved_ratio")
tf_conserv_data$TF = unique(full_table$TF)

for (tf in unique(full_table$TF)){
  tf_temp_data = full_table[full_table$TF == tf,]
  chip_conf_tf_temp_data = conserv_conf_full_table[conserv_conf_full_table$TF == tf,]
  #conf_perc = nrow(chip_conf_tf_temp_data)/nrow(tf_temp_data)
  #tf_conf_data[tf_conf_data$TF == tf,]$Bind_num = nrow(tf_temp_data)
  #tf_conf_data[tf_conf_data$TF == tf,]$ChiP_bind_num = nrow(chip_conf_tf_temp_data)
  #tf_conf_data[tf_conf_data$TF == tf,]$ChIP_conf_perc = conf_perc
  chip_conf_tf_temp_data = unique(chip_conf_tf_temp_data[,c("chr","start","end")])
  tf_temp_data = unique(tf_temp_data[,c("chr","start","end")])
  conf_perc = nrow(chip_conf_tf_temp_data)/nrow(tf_temp_data)
  tf_conserv_data[tf_conserv_data$TF == tf,]$Bind_num = nrow(tf_temp_data)
  tf_conserv_data[tf_conserv_data$TF == tf,]$Conserved_bind_num = nrow(chip_conf_tf_temp_data)
  tf_conserv_data[tf_conserv_data$TF == tf,]$Conserved_ratio = conf_perc
}

#The output table of conservation analysis
tf_conf_file = "~/rprojects/review_data/results/new_generated_res_test/aff_tf_conserved.txt"
write.table(file = tf_conf_file, tf_conserv_data,sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
