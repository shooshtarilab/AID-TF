library(ggplot2)
library(data.table)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(EnvStats)

library(stringr)

#Entering the trait of interest
trait = "RA"

print(trait)

#Address of the file having all of the immunochip SNPs 
#of the trait
ref_dir = paste0("/home/nader95/rprojects/tfproject/immunochip_data/",trait,"/allSNPs-Sorted.bed")
ref_snp = read.table(ref_dir, sep = '\t', header = FALSE)



mod ="all_snps"

#Directory of the immunochip data of the trait
direct = paste0("/home/nader95/rprojects/tfproject/immunochip_data/",trait)

#Name of the file having Ci-SNPs
snp_ext = "/credSNPs-v2-Sorted.bed"

ne = "chr"


dir = paste0(direct,snp_ext)

data = read.table(dir, sep = '\t', header = FALSE)


chr_start = data$V1
snp_start = data$V2


total_snp_number = length(snp_start)
total_ref_number = length(ref_snp$V1)

x_test = data.frame(matrix(ncol = 2, nrow = length(snp_start)))
colnames(x_test) = c('V1','V2')
x_test$V1 = chr_start
x_test$V2 = snp_start

x_ref = data.frame(matrix(ncol = 2, nrow = length(ref_snp$V3)))
colnames(x_ref) = c('V1','V2')
x_ref$V1 = ref_snp$V1
x_ref$V2 = ref_snp$V3

x_range_test = IRanges(start = x_test$V2, end = x_test$V2+1)
x_gr_test = GRanges(seqnames = x_test$V1, ranges = x_range_test)

x_range_ref = IRanges(start = x_ref$V2-1, end = x_ref$V2)
x_gr_ref = GRanges(seqnames = x_ref$V1, ranges = x_range_ref)

x_gr_ref = setdiff(x_gr_ref,x_gr_test)

total_snp_number = length(x_gr_test)
total_ref_number = length(x_gr_ref)


#Address of the file having the list of samples' DHS data
list_dir = "/scratch-deleted-2021-mar-20/nader95/rprojects/tfproject/list/list.txt"
list = read.table(list_dir, sep = "\t", header = FALSE)
list = list$V1
list = list[1:376]


#Address of the file comprising P values of enrichment of
#samples for the trait of interest from the first step 
pval_dir = paste0("/home/nader95/rprojects/tfproject/results/immunochip_res/first_test/",trait,"/cell_enrich.rds")
real_pval = readRDS(pval_dir)
ordered_index = order(real_pval)
real_pval = real_pval[ordered_index]
list = list[ordered_index]
list_index = -log10(real_pval) > -log10(0.05/376)
list = list[list_index]

print("for will start")
print(list)

new_frame = data.frame(matrix(ncol = 10, nrow = 1e5))
colnames(new_frame) = c("cell_type", "motif", "chr" ,"SNP_start","SNP_end","motif_start","motif_end","motif_pvalue","cell_pvalue","SNP_id")


great_frame = data.frame(matrix(ncol = 3, nrow = 1e5))
colnames(great_frame) = c("chr" ,"start","end")



count = 1

for (i in c(1:length(list))){
    
    #The address of the DHS data of a sample
    data_temp_dir = paste0("/scratch-deleted-2021-mar-20/nader95/rprojects/tfproject/data/",list[i])

    data_temp = read.table(data_temp_dir, sep = "\t", header = FALSE)

    show_list = str_remove(list[i],".bed")
    show_list = str_remove(show_list, "EpiUwRmapDNase")

    #Address of the file having the list of affected
    #motifs within the sample for the trait of interest
    motif_list = readRDS(paste0("~/rprojects/tfproject/results/immunochip_res/third_test/",trait,"/",list[i],"/","motif_list.rds"))
    
    #Address of the file having the P values of enrichment of
    #motifs within the sample for the trait of interest
    motif_enrich = readRDS(paste0("~/rprojects/tfproject/results/immunochip_res/third_test/",trait,"/",list[i],"/","cell_enrich.rds"))

    ordered_index = order(motif_enrich)
    motif_enrich = motif_enrich[ordered_index]
    motif_list = motif_list[ordered_index]

    motif_enrich1 = p.adjust(motif_enrich,method = "fdr")

    if (sum(motif_enrich1 < 0.1) == 0){
        next
    }
    
    motif_enrich = motif_enrich[motif_enrich1 < 0.1]
    motif_list = motif_list[motif_enrich1 < 0.1]
                  
    cell_enrich = c()
    ref_num = c()
    test_num = c()

    for (j in c(1:length(motif_list))){
        found_index = data_temp$V4 == motif_list[j]
        data_found = data_temp[found_index,]
        effect_index = data_found$V10 != 0
        #print(effect_index)
        data_temp_filt = data_found[effect_index,]
        data_range = IRanges(start = data_temp_filt$V10 - 1, end = data_temp_filt$V10)
        data_gr = GRanges(seqnames = data_temp_filt$V1, ranges = data_range)
        #data_gr = unique(data_gr)
        motif_range = IRanges(start = data_temp_filt$V2, end = data_temp_filt$V3)
        motif_gr = GRanges(seqnames = data_temp_filt$V1, ranges = motif_range)

        hit_index = findOverlaps(x_gr_test, data_gr)
        snp_index = queryHits(hit_index)
        hit_index = subjectHits(hit_index)
        hit_num = length(hit_index)

        snp_name = data[snp_index,]$V4
        print(snp_name)
        print(paste0("number of regions = ",hit_num))
        #print(hit_num)
        #print(count)

        if (hit_num >0){
            new_frame[count:(count+hit_num-1),]$cell_type = rep(show_list,hit_num)
            new_frame[count:(count+hit_num-1),]$motif = rep(motif_list[j],hit_num)
            new_frame[count:(count+hit_num-1),]$chr = seqnames(data_gr[hit_index,])
            new_frame[count:(count+hit_num-1),]$SNP_id = snp_name            
            
            new_frame[count:(count+hit_num-1),]$SNP_start = start(data_gr[hit_index,])
            new_frame[count:(count+hit_num-1),]$SNP_end = end(data_gr[hit_index,])
            new_frame[count:(count+hit_num-1),]$motif_start = start(motif_gr[hit_index,])
            new_frame[count:(count+hit_num-1),]$motif_end = end(motif_gr[hit_index,])

            new_frame[count:(count+hit_num-1),]$motif_pvalue = rep(motif_enrich[j],hit_num)
            new_frame[count:(count+hit_num-1),]$cell_pvalue = rep(real_pval[i],hit_num)


            great_frame[count:(count+hit_num-1),]$chr = seqnames(data_gr[hit_index,])
            great_frame[count:(count+hit_num-1),]$start = start(motif_gr[hit_index,])
            great_frame[count:(count+hit_num-1),]$end = end(motif_gr[hit_index,])


            count = count + hit_num
        }
        
        print('s')
        }
        
        
        
  
  print(i)
  
}


#Directory of the file in which the main results table
#and the greatpathway table (genomic position of affected 
#regions) are written (for the trait of interest)

main_dir = paste0("~/rprojects/tfproject/results/immunochip_res/new_generated_res/",trait)
dir.create(main_dir)

#Writing the main results table having all of
#the results for the trait of interest (subsections of Supp. table S1)
nonna_index = ! is.na(new_frame$motif)
new_frame = new_frame[nonna_index,]
save_file = paste0(main_dir,"/final_total_list.txt")

file.remove(save_file)
write.table(new_frame,file=save_file,quote=F,sep="\t",col.names=T,row.names=F)

#Writing the table used for great pathway analysis 
#comprising the genomic locations of all affected regions
nonna_index = ! is.na(great_frame$chr)
great_frame = great_frame[nonna_index,]
save_file = paste0(main_dir,"/final_great_list.txt")

file.remove(save_file)
write.table(great_frame,file=save_file,quote=F,sep="\t",col.names=T,row.names=F)



 

