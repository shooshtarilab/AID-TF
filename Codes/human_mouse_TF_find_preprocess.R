library(ape)
library(dplyr)
library(tidyr)
library(stringr)
library(liftOver)

full_table_dir = "~/rprojects/review_data/results/new_generated_res_test/final_dataframe.txt"

full_table = read.table(file = full_table_dir, sep = "\t", 
                        header = TRUE)

df_separated <- full_table %>%
  separate_rows(TF, sep = "/")
df_separated <- df_separated %>%
  separate_rows(TF, sep = ":")
df_separated = df_separated[df_separated$TF != "",]

new_tf_names = toupper(str_split(df_separated$TF, pattern = "_RC",
                                 simplify = TRUE)[,1])
df_separated$TF = new_tf_names

tf_list = unique(df_separated$TF)

human_list = c("HAP1","MYC","ARNT","TLX1","NFIC","CTF1",
               "SPI-B","ELF5","GABP-BETA","ELF4","ELF2",
               "SPI1","CREB","FHL1","CREM","E2F-1","AHR",
               "NKX3-1","OSF2")
non_found_tfs = tf_list[!tf_list %in% human_list]

motif_list = unique(df_separated$motif[df_separated$TF %in% 
                                       human_list])

found_motifs = list.files(path = "~/rprojects/review_data/motif_data/")
found_motifs = str_split(found_motifs, pattern = ".bed",
                         simplify = TRUE)[,1]

missing_motifs = motif_list[!(motif_list %in% found_motifs)]
