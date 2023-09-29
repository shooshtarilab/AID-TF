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


test_num = c()
ref_num = c()
cell_enrich = c()

#Address of the file having the list of samples' DHS data
list_dir = "/scratch/nader95/rprojects/tfproject/list/list.txt"
list = read.table(list_dir, sep = "\t", header = FALSE)
list = list$V1
list = list[1:376]

#Address of the file having enriched samples
#for the trait of interest
list_dir = paste0("/home/nader95/rprojects/tfproject/results/cell_type_res/",
                  trait,"/enrich_list.rds")
list_index = readRDS(list_dir)

list = list[list_index]

print("for will start")
print(list)

for (i in c(1:length(list))){
    
    #The address of the DHS data of a sample
    data_temp_dir = paste0("/scratch/nader95/rprojects/tfproject/data/",list[i])

    data_temp = read.table(data_temp_dir, sep = "\t", header = FALSE)


    motif_list = unique(data_temp$V4)
    cell_enrich = c()
    ref_num = c()
    test_num = c()

    for (j in c(1:length(motif_list))){
        found_index = data_temp$V4 == motif_list[j]
        data_found = data_temp[found_index,]

        data_temp_filt = data_found
        data_range = IRanges(start = data_temp_filt$V2, end = data_temp_filt$V3)
        data_gr = GRanges(seqnames = data_temp_filt$V1, ranges = data_range)
        data_gr = disjoin(unique(data_gr))

        a = 0
        b = 0
        c = 0
        d = 0
        n = 0
        mat = matrix(0, nrow = 2, ncol = 2)
        res = 0


        hit_index = findOverlaps(x_gr_test, data_gr)
        test_hit_number = length(unique(queryHits(hit_index)))
        
        ref_index = findOverlaps(x_gr_ref, data_gr)
        ref_hit_number = length(unique(queryHits(ref_index)))
        
        test_num = append(test_num, test_hit_number)
        ref_num = append(ref_num, ref_hit_number)
        
        a = test_hit_number
        c = total_snp_number - a
        b = ref_hit_number
        d = total_ref_number - b
        n = total_ref_number + total_snp_number
        
        
        
        mat[1,1] = a
        mat[2,1] = c
        mat[1,2] = b
        mat[2,2] = d
        
        res = fisher.test(mat)
        cell_enrich = append(cell_enrich, res$p.value)

        print("s")
    }

    #Creating a directory for the second test results of
    #the trait of interest
    main_dir1 = paste0("/home/nader95/rprojects/tfproject/results/motif_type_corrected/",trait,"/")
    unlist(main_dir1, recursive = TRUE)
    dir.create(main_dir1)
    
    #Creating a subdirectory for each sample
    main_dir = paste0(main_dir1,list[i],"/")
    unlist(main_dir, recursive = TRUE)
    dir.create(main_dir)

    #Writing the P values of the enrichment of motifs within 
    #the sample for the trait of interest
    saveRDS(cell_enrich, file = paste0(main_dir,"cell_enrich.rds"))
    #Writing the list of motifs within the sample
    saveRDS(motif_list, file = paste0(main_dir,"motif_list.rds"))
    

    print(i)


}







        

        














        
