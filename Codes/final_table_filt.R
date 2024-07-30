library(GenomicRanges)
library(stringr)
library(rtracklayer)
library(dplyr)
library(liftOver)
library(xlsx)
library(tidyr)


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
all_chip_tf_list = unique(c(chip_frame$TF, chip_frame_immune$TF))

chip_data_all = list()

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

#Conservation file
cons_file = "~/rprojects/review_data/Conserved_hg19/hg19.LECIFv1.1.bw"
cons_data = import.bw(cons_file)

cons_data_filt = cons_data[cons_data$score>0.2]

human_list = c("HAP1","MYC","ARNT","TLX1","NFIC","CTF1",
               "SPI-B","ELF5","GABP-BETA","ELF4","ELF2",
               "SPI1","CREB","FHL1","CREM","E2F-1","AHR",
               "NKX3-1","OSF2")

full_table_dir = "~/rprojects/review_data/results/new_generated_res_test/Supp-Tables/Supp.TableS4.xlsx"
full_table = read.xlsx(file = full_table_dir, 1, header = TRUE)

full_table = full_table[!is.na(full_table$start),]

df_separated <- full_table %>%
  separate_rows(TF, sep = "/")
df_separated <- df_separated %>%
  separate_rows(TF, sep = ":")
df_separated = df_separated[df_separated$TF != "",]

new_tf_names = toupper(str_split(df_separated$TF, pattern = "_RC",
                                 simplify = TRUE)[,1])
df_separated$TF = new_tf_names

df_separated = df_separated[df_separated$TF %in% human_list,]

df_separated[["ChIP_confirmed"]] = "N"
df_separated[["Mouse_conserved"]] = "N"

#Confirming TF binding sites using ChIP-seq data
for (tf in names(chip_data_all)){
  tf_temp_data = df_separated[df_separated$TF == tf,]
  chip_tf_data = chip_data_all[[tf]]
  
  tf_temp_grg = GRanges(ranges = IRanges(
    start = tf_temp_data$start,
    end = tf_temp_data$end
  ),seqnames = tf_temp_data$chr)
  
  tf_conf_ind = unique(queryHits(findOverlaps(tf_temp_grg,
                                              chip_tf_data)))
  
  tf_temp_data$ChIP_confirmed[tf_conf_ind] = 'Y'
  df_separated[df_separated$TF == tf,] = tf_temp_data
  print(tf)
}

df_grg = GRanges(ranges = IRanges(start=df_separated$start,
                                  end=df_separated$end),
                 seqnames = df_separated$chr)


df_cons_idx = unique(queryHits(findOverlaps(df_grg,
                                            cons_data_filt)))

df_separated[df_cons_idx,]$Mouse_conserved = "Y"

out_table = "~/rprojects/review_data/results/new_generated_res_test/new_tables/supp/Supp.TableS4_new.txt"
write.table(df_separated, file = out_table, sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)


##########################################################
full_table_dir = "~/rprojects/review_data/results/new_generated_res_test/Supp-Tables/Supp.TableS7.xlsx"
full_table = read.xlsx(file = full_table_dir, 1, header = TRUE)

full_table = full_table[!is.na(full_table$Chromosome),]

df_separated <- full_table %>%
  separate_rows(Transcription.Factors, sep = "/")
df_separated <- df_separated %>%
  separate_rows(Transcription.Factors, sep = ":")
df_separated = df_separated[df_separated$Transcription.Factors != "",]

new_tf_names = toupper(str_split(df_separated$Transcription.Factors, pattern = "_RC",
                                 simplify = TRUE)[,1])
df_separated$Transcription.Factors = new_tf_names

df_separated = df_separated[df_separated$Transcription.Factors %in% human_list,]

df_separated[["ChIP_confirmed"]] = "N"
df_separated[["Mouse_conserved"]] = "N"

#Confirming TF binding sites using ChIP-seq data
for (tf in names(chip_data_all)){
  tf_temp_data = df_separated[df_separated$Transcription.Factors == tf,]
  chip_tf_data = chip_data_all[[tf]]
  
  tf_temp_grg = GRanges(ranges = IRanges(
    start = tf_temp_data$Start.Position.of.the.TF.Binding.Site,
    end = tf_temp_data$End.Position.of.the.TF.Binding.Site
  ),seqnames = tf_temp_data$Chromosome)
  
  tf_conf_ind = unique(queryHits(findOverlaps(tf_temp_grg,
                                              chip_tf_data)))
  
  tf_temp_data$ChIP_confirmed[tf_conf_ind] = 'Y'
  df_separated[df_separated$Transcription.Factors == tf,] = tf_temp_data
  print(tf)
}

df_grg = GRanges(ranges = IRanges(start=df_separated$Start.Position.of.the.TF.Binding.Site,
                                  end=df_separated$End.Position.of.the.TF.Binding.Site),
                 seqnames = df_separated$Chromosome)


df_cons_idx = unique(queryHits(findOverlaps(df_grg,
                                            cons_data_filt)))

df_separated[df_cons_idx,]$Mouse_conserved = "Y"

out_table = "~/rprojects/review_data/results/new_generated_res_test/new_tables/supp/Supp.TableS7_new.txt"
write.table(df_separated, file = out_table, sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

#############################################################
full_table_dir = "~/rprojects/review_data/results/new_generated_res_test/Main-Tables/Table-2.xlsx"
full_table = read.xlsx(file = full_table_dir, 1, header = TRUE)

full_table = full_table[!is.na(full_table$Chromosome),]

df_separated <- full_table %>%
  separate_rows(Transcription.Factors, sep = "/")
df_separated <- df_separated %>%
  separate_rows(Transcription.Factors, sep = ":")
df_separated = df_separated[df_separated$Transcription.Factors != "",]

new_tf_names = toupper(str_split(df_separated$Transcription.Factors, pattern = "_RC",
                                 simplify = TRUE)[,1])
df_separated$Transcription.Factors = new_tf_names

df_separated = df_separated[df_separated$Transcription.Factors %in% human_list,]

df_separated[["ChIP_confirmed"]] = "N"
df_separated[["Mouse_conserved"]] = "N"

#Confirming TF binding sites using ChIP-seq data
for (tf in names(chip_data_all)){
  tf_temp_data = df_separated[df_separated$Transcription.Factors == tf,]
  chip_tf_data = chip_data_all[[tf]]
  
  tf_temp_grg = GRanges(ranges = IRanges(
    start = tf_temp_data$Start.Position.of.the.TF.Binding.Site,
    end = tf_temp_data$End.Position.of.the.TF.Binding.Site
  ),seqnames = tf_temp_data$Chromosome)
  
  tf_conf_ind = unique(queryHits(findOverlaps(tf_temp_grg,
                                              chip_tf_data)))
  
  tf_temp_data$ChIP_confirmed[tf_conf_ind] = 'Y'
  df_separated[df_separated$Transcription.Factors == tf,] = tf_temp_data
  print(tf)
}

df_grg = GRanges(ranges = IRanges(start=df_separated$Start.Position.of.the.TF.Binding.Site,
                                  end=df_separated$End.Position.of.the.TF.Binding.Site),
                 seqnames = df_separated$Chromosome)


df_cons_idx = unique(queryHits(findOverlaps(df_grg,
                                            cons_data_filt)))

df_separated[df_cons_idx,]$Mouse_conserved = "Y"

out_table = "~/rprojects/review_data/results/new_generated_res_test/new_tables/main/Table-2.txt"
write.table(df_separated, file = out_table, sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
