#DEPENDENCIES - START
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(DEGreport)
library(pheatmap)
library(GGally)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(cowplot)
#DEPENDENCIES - END


###############----------------------------------------------------###########################################################
###############-----------COUNTS_START-----------------------------###########################################################
###############----------------------------------------------------###########################################################

#DIRECTORY ORGANISATION AND READING IN THE DATA

#Data preparation
setwd("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set")
#######SGE READ PREPARATION 
#makes lists of file locations to unzip, then run the output list of commands using termincal
gen.count.locations <- function(exon,sg) {
  paste0("cd /Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set/QUANTS/",exon,"_",sg,"/results/pycroquet
gunzip *.query_to_library_counts.tsv.gz")
}
#Function to apply the above function to listed exons and designated sgRNA_A
a_locs<-lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(y) {gen.count.locations(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
b_locs<-lapply (c("3_2","4","5_2"), function(y) {gen.count.locations(exon=y, sg="B")
})
#make a file of all the files directories to gunzip
lapply(a_locs, write, "./a_locs.txt", append=TRUE)
lapply(b_locs, write, "./b_locs.txt", append=TRUE)

#######PLASMID READ PREPARATION 
#makes lists of file locations to unzip, then run the output list of commands using termincal
gen.count.locations <- function(exon,sg) {
  paste0("cd /Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/data/QUANTS_plasmids/45336/",exon,"_",sg,"/results/pycroquet
gunzip *.query_to_library_counts.tsv.gz")
}
#Function to apply the above function to listed exons and designated sgRNA_A
a_locs<-lapply (c("1_1","1_2","2_1","2_2","3_1","3_2","4","5_1","5_2","6","7","8","9"), function(y) {gen.count.locations(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
b_locs<-lapply (c("1_1","1_2","2_1","2_2","3_1","3_2","4","5_1","5_2","6","7","8","9"), function(y) {gen.count.locations(exon=y, sg="B")
})
#make a file of all the files directories to gunzip
lapply(a_locs, write, "./a_locs_plasmids.txt", append=TRUE)
lapply(b_locs, write, "./b_locs_plasmids.txt", append=TRUE)

#Function to read in separate count.csv files, $PATH will need to be changed. 
read.counts <- function(exon,sg) {
  #directories 
  dir=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set/QUANTS/",exon,"_",sg,"/results/pycroquet/")
  dir_plasmid=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set/QUANTS_plasmids/45336/",exon,"_",sg,"/results/pycroquet/")
  out=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set/data/", "E", exon, "_SG", sg, "_count_frame.csv")
  #makes a list of count file names held in the dir 
  temp = list.files(dir, pattern="*rep_1.query_to_library_counts.tsv|rep_2.query_to_library_counts.tsv|rep_3.query_to_library_counts.tsv|R1.query_to_library_counts.tsv|R2.query_to_library_counts.tsv|R3.query_to_library_counts.tsv")
  temp_plasmid = list.files(dir_plasmid, pattern="*R1.query_to_library_counts.tsv|R2.query_to_library_counts.tsv|R3.query_to_library_counts.tsv")
  #makes a list of full paths to each file
  file_all_count=paste0(dir,temp)
  file_all_count_plasmid=paste0(dir_plasmid,temp_plasmid)
  #reads the data from these paths into a list of dataframes
  myfiles = lapply(file_all_count, read.table, skip = 2, header =T, comment.char = "")
  myfiles_plasmid = lapply(file_all_count_plasmid, read.table, skip = 2, header =T, comment.char = "", fill=TRUE)
  myfiles_plasmid = lapply(myfiles_plasmid,"[",c(1,5))
  #makes a large dataframe with all counts - n.b. be careful of the ordering of the new column names make sure agrees to the order in the lists above - this version is ascending numerical
  df <- myfiles %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("rep")|contains("RAD51C"))
  df_plasmid <- myfiles_plasmid %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("RAD51C"))
  df<- df %>% left_join(df_plasmid, by=c("X.id"))
  df <- df %>% `colnames<-` (c("id", "D14R1","D14R2","D14R3","D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","PLASMID_R1","PLASMID_R2","PLASMID_R3")) %>% replace(is.na(.), 0)
  df$PLASMID_MEAN <- (df$PLASMID_R1 + df$PLASMID_R2 + df$PLASMID_R3) / 3
  #saves dataframe to output location and re-imports it as a named df with names taken from definitions above
  write.csv(df, file=out, row.names = FALSE)
  x=paste0("count_table_", "E", exon, "_SG", sg)
  assign(x,read.csv(out), envir = .GlobalEnv)
}

###############-----------RUN FOR ALL LIBRARIES_STRART-------------###########################################################
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(y) {read.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(y) {read.counts(exon=y, sg="B")
})
###############-----------RUN FOR ALL LIBRARIES_END----------------###########################################################

D4_average <- (count_table_E1_2_SGA$D4R1 + count_table_E1_2_SGA$D4R2 + count_table_E1_2_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E2_1_SGA$D4R1 + count_table_E2_1_SGA$D4R2 + count_table_E2_1_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E2_2_SGA$D4R1 + count_table_E2_2_SGA$D4R2 + count_table_E2_2_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E3_1_SGA$D4R1 + count_table_E3_1_SGA$D4R2 + count_table_E2_2_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E3_2_SGB$D4R1 + count_table_E3_2_SGB$D4R2 + count_table_E3_2_SGB$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E4_SGB$D4R1 + count_table_E4_SGB$D4R2 + count_table_E4_SGB$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E5_1_SGA$D4R1 + count_table_E5_1_SGA$D4R2 + count_table_E5_1_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E5_2_SGB$D4R1 + count_table_E5_2_SGB$D4R2 + count_table_E5_2_SGB$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E6_SGA$D4R1 + count_table_E6_SGA$D4R2 + count_table_E6_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E7_SGA$D4R1 + count_table_E7_SGA$D4R2 + count_table_E7_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E8_SGA$D4R1 + count_table_E8_SGA$D4R2 + count_table_E8_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=5000)

D4_average <- (count_table_E9_SGA$D4R1 + count_table_E9_SGA$D4R2 + count_table_E9_SGA$D4R3)/3
hist(D4_average, xlim=c(0,5000), breaks=100)


###############----------------------------------------------------###########################################################
###############-----------COUNTS_END-------------------------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------VCF_CLEANING_FOR_VEP_INPUT_START---------###########################################################
###############----------------------------------------------------###########################################################

#Define directories where VaLiAnT VCFs are stored
dir_vcf=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/data/vcfs/")
vcf_out=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/data/vep_inputs/")
#makes a list of count file names held in the dir 
temp_vcf = list.files(dir_vcf, pattern="*.vcf")
#makes a list of full paths to each file
file_all_vcf=paste0(dir_vcf,temp_vcf)
#reads the data from these paths into a list of dataframes
my_vcf_files = sapply(file_all_vcf, read.table, skip = 9, header =T, comment.char = "#", simplify=FALSE)
#names columns
colnames_vcf <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
my_vcf_files = lapply(my_vcf_files, setNames, colnames_vcf)
#gets list of directory names so one knows which targetons are which dataframe
vcf_names <-names(my_vcf_files)
#gets file names from full paths
vcf_names_base <-basename(vcf_names)
rm(vcf_names)

#########################################################################
#FUNCTIONS
#########################################################################

#Function to extract the oligo ID from the INFO field and place in the ID column
wrangle.vcf<-function(number) {
  tmp_vcf_data<-my_vcf_files[[number]]
  tmp_vcf_data[c('INFO1', 'INFO2', 'INFO3')] <- str_split_fixed(tmp_vcf_data$INFO, ';', 3)
  tmp_vcf_data$INFO2<-str_remove(tmp_vcf_data$INFO2, "SGE_OLIGO=")
  tmp_vcf_data[ ,c('INFO1', 'INFO3', 'ID')] <- list(NULL)
  tmp_vcf_data$ID = tmp_vcf_data$INFO2
  tmp_vcf_data$INFO2 = NULL
  tmp_vcf_data <- tmp_vcf_data %>% relocate(ID, .before = REF)
  colnames_vcf_out <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  setNames(tmp_vcf_data, colnames_vcf_out)
}

#run the above function over all VCFs
my_vcf_files = sapply (c(1:26), function(z) {wrangle.vcf(number=z)
}, simplify = FALSE, USE.NAMES=TRUE)
#name the dataframes based on the original VCF file 
names(my_vcf_files) <- vcf_names_base
#function to output each VCF as a separate dataframe
make.tables<-function(object_name) {
  table_output_vcf<-my_vcf_files[[object_name]]
  suffix<-names(my_vcf_files)[object_name]
  write.table(table_output_vcf, file=paste0(vcf_out,"VEP_INPUT_", suffix), sep = "\t", row.names = FALSE, quote = FALSE)
}
#run the above function over all VCFs to output appropriately named VEP INPUTS 
lapply (c(1:26), function(zy) {make.tables(object_name=zy)
})

###############----------------------------------------------------###########################################################
###############-----------VCF_CLEANING_FOR_VEP_INPUT_END-----------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------COUNT_FILTERING_INPUT_START--------------###########################################################
###############----------------------------------------------------###########################################################

#directoy of metafiles - ONLY A at the moment (expand 13 to 26 if want to do all)
dir_meta=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set/meta/")
#makes a list of count file names held in the dir 
temp_meta = list.files(dir_meta, pattern="*.csv")
#makes a list of full paths to each file
file_all_meta=paste0(dir_meta,temp_meta)
#reads the data from these paths into a list of dataframes
my_meta_files = sapply(file_all_meta, read.csv, header =T, comment.char = "#", simplify=FALSE)
#meta names, full paths
meta_names <-names(my_meta_files)
#gets file names from full paths
meta_names_base <-basename(meta_names)
rm(meta_names)
#name the dataframes based on the original VCF file 
names(my_meta_files) <- meta_names_base
#function to output each VCF as a separate dataframe
make.tables<-function(object_name) {
  table_output_meta<-my_meta_files[[object_name]]
  identifier<-names(my_meta_files)[object_name]
  names(table_output_meta)[names(table_output_meta) == 'oligo_name'] <- 'id'
  assign(identifier, table_output_meta, envir = .GlobalEnv)
}
#run the above function over all metas to output appropriately named metas into the global environment
lapply (c(1:26), function(zy) {make.tables(object_name=zy)
})

#merge meta and counts_SGA
E1_2_SGA_master<-merge(count_table_E1_2_SGA, chr17_58692607_58692852_plus_sgRNA_1_2_a_meta.csv, by="id", all.x=FALSE)
E2_1_SGA_master<-merge(count_table_E2_1_SGA, chr17_58694822_58695075_plus_sgRNA_2_1_a_meta.csv, by="id", all.x=FALSE)
E2_2_SGA_master<-merge(count_table_E2_2_SGA, chr17_58694985_58695224_plus_sgRNA_2_2_a_meta.csv, by="id", all.x=FALSE)
E3_1_SGA_master<-merge(count_table_E3_1_SGA, chr17_58696626_58696865_plus_sgRNA_3_1_a_meta.csv, by="id", all.x=FALSE)
E5_1_SGA_master<-merge(count_table_E5_1_SGA, chr17_58709712_58709961_plus_sgRNA_5_1_a_meta.csv, by="id", all.x=FALSE)
E6_SGA_master<-merge(count_table_E6_SGA, chr17_58720663_58720902_plus_sgRNA_6_a_meta.csv, by="id", all.x=FALSE) 
E7_SGA_master<-merge(count_table_E7_SGA, chr17_58723963_58724206_plus_sgRNA_7_a_meta.csv, by="id", all.x=FALSE)                     
E8_SGA_master<-merge(count_table_E8_SGA, chr17_58732410_58732651_plus_sgRNA_8_a_meta.csv, by="id", all.x=FALSE)                      
E9_SGA_master<-merge(count_table_E9_SGA, chr17_58734062_58734302_plus_sgRNA_9_a_meta.csv, by="id", all.x=FALSE)                     


#merge meta and counts_SGB

E3_2_SGB_master<-merge(count_table_E3_2_SGB, chr17_58696699_58696943_plus_sgRNA_3_2_b_meta.csv, by="id", all.x=FALSE)
E4_SGB_master<-merge(count_table_E4_SGB, chr17_58703145_58703389_plus_sgRNA_4_b_meta.csv, by="id", all.x=FALSE)
E5_2_SGB_master<-merge(count_table_E5_2_SGB, chr17_58709875_58710114_plus_sgRNA_5_2_b_meta.csv, by="id", all.x=FALSE)



#remove rows that have fewer than 10 counts total across all 9 samples - OPTIONAL REMOVE VARIANTS IN CONSTANT REGTIONS
filter.counts<-function(exon, sg){
  master_input<-paste0("E",exon,"_SG",sg,"_master")
  filtered_counts <- get(master_input)
  filtered_counts$rowsum <- rowSums(filtered_counts[, c("D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D14R1","D14R2","D14R3")])
  filtered_counts<- filtered_counts %>% filter(rowsum>=10)
  #optional_remove high counts
  #filtered_counts<- filtered_counts %>% filter(rowsum<=10000)
  #OPTIONAL - REMOVE VARIANTS PAM MUTATION CODONS
  #filtered_counts<- filtered_counts %>% filter(pam_mut_sgrna_id %in% "" | is.na(pam_mut_sgrna_id))
  #OPTIONAL - REMOVE VARIANTS IN CONSTANT REGIONS - CAN REMOVE THIS STEP IF ADAPTERS STRIPPED SAME AS SEQUENCING PRIMERS
  filtered_counts<- filtered_counts %>% filter(vcf_var_in_const==0 | is.na(vcf_var_in_const))
  tablename<-paste0(master_input,"_filtered_counts")
  assign(tablename, filtered_counts, envir = .GlobalEnv)
}


#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(y) {filter.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(y) {filter.counts(exon=y, sg="B")
})


###############----------------------------------------------------###########################################################
###############-----------COUNT_FILTERING_INPUT_END----------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------VEP_START--------------------------------###########################################################
###############----------------------------------------------------###########################################################

#mergeing vep outputs to create vep + meta
#Define directories where VEP OUTPUTS are stored
dir_vep=paste0("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/full_analysis_set/vep_outputs/")
#makes a list of count file names held in the dir 
temp_vep = list.files(dir_vep, pattern="*.tsv")
#makes a list of full paths to each file
file_all_vep=paste0(dir_vep,temp_vep)
#reads the data from these paths into a list of dataframes
my_vep_files = sapply(file_all_vep, read.table, skip = 92, header =T, comment.char = "", simplify=FALSE)
#meta names, full paths
vep_names <-names(my_vep_files)
#gets file names from full paths
vep_names_base <-basename(vep_names)
rm(vep_names)
#name the dataframes based on the original VEP file 
names(my_vep_files) <- vep_names_base
#function to output each VCF as a separate dataframe
make.vep.tables<-function(object_name) {
  table_output_vep<-my_vep_files[[object_name]]
  identifier<-names(my_vep_files)[object_name]
  names(table_output_vep)[names(table_output_vep) == 'X.Uploaded_variation'] <- 'id'
  table_output_vep <- table_output_vep %>% filter(Feature %in% "ENST00000337432.9")
  assign(paste0("VEP_OUTPUT_",identifier), table_output_vep, envir = .GlobalEnv)
}
#run the above function over all metas to output appropriately named metas into the global environment
lapply (c(1:26), function(zy) {make.vep.tables(object_name=zy)
})

#merge filtered master tables with VEP annotation tables and syn_filter
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGA 
E1_2_SGA_master_annotated<-merge(E1_2_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58692607_58692852_plus_sgRNA_1_2_a.tsv, by="id", all.x=TRUE)
E2_1_SGA_master_annotated<-merge(E2_1_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58694822_58695075_plus_sgRNA_2_1_a.tsv, by="id", all.x=TRUE)
E2_2_SGA_master_annotated<-merge(E2_2_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58694985_58695224_plus_sgRNA_2_2_a.tsv, by="id", all.x=TRUE)
E3_1_SGA_master_annotated<-merge(E3_1_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58696626_58696865_plus_sgRNA_3_1_a.tsv, by="id", all.x=TRUE)
E5_1_SGA_master_annotated<-merge(E5_1_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58709712_58709961_plus_sgRNA_5_1_a.tsv, by="id", all.x=TRUE)
E6_SGA_master_annotated<-merge(E6_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58720663_58720902_plus_sgRNA_6_a.tsv, by="id", all.x=TRUE)                                 
E7_SGA_master_annotated<-merge(E7_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58723963_58724206_plus_sgRNA_7_a.tsv, by="id", all.x=TRUE)                                
E8_SGA_master_annotated<-merge(E8_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58732410_58732651_plus_sgRNA_8_a.tsv, by="id", all.x=TRUE)                                
E9_SGA_master_annotated<-merge(E9_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58734062_58734302_plus_sgRNA_9_a.tsv, by="id", all.x=TRUE)                               
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGB
E3_2_SGB_master_annotated<-merge(E3_2_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58696699_58696943_plus_sgRNA_3_2_b.tsv, by="id", all.x=TRUE)
E4_SGB_master_annotated<-merge(E4_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58703145_58703389_plus_sgRNA_4_b.tsv, by="id", all.x=TRUE)
E5_2_SGB_master_annotated<-merge(E5_2_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr17_58709875_58710114_plus_sgRNA_5_2_b.tsv, by="id", all.x=TRUE)

###############----------------------------------------------------###########################################################
###############-----------VEP_END----------------------------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------DESEQ2_START-----------------------------###########################################################
###############----------------------------------------------------###########################################################

#filter the above tables to get normalization tables
normalization.counts<-function(exon, sg){
  master_input<-paste0("E",exon,"_SG",sg,"_master_annotated")
  filtered_counts <- get(master_input)
  filtered_counts <-filtered_counts[!duplicated(filtered_counts$mseq), ]
  filtered_counts<- filtered_counts %>% filter(Consequence %in% "synonymous_variant" | Consequence %in% "intron_variant")
  tablename<-paste0("E",exon,"_SG",sg,"_normalization_counts")
  assign(tablename, filtered_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(y) {normalization.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(y) {normalization.counts(exon=y, sg="B")
})


#make normalization matricies
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D14R1","D14R2","D14R3")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(y) {deseq.norm.input(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(y) {deseq.norm.input(exon=y, sg="B")
})


#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D14R1","D14R2","D14R3")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(y) {deseq.count.input(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(y) {deseq.count.input(exon=y, sg="B")
})


#read in experimental designs
annot_a<-read.csv("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/data/annot_a.csv")
annot_a_continuous<-read.csv("/Users/aw28/Documents/SGE_PROJECTS/RAD51C_rebecca/data/annot_a_continuous.csv")
annot_a<-as.matrix.data.frame(annot_a)

#########################################################################
#FUNCTIONS_DESEQ_START
#########################################################################

#OUTPUT FUNCTIONS 
#PCA SCREE PLOT
plotPCA.hk <- function (object, intgroup = "condition", ntop = 500, pc_1 = 1, pc_2 = 2, returnData = FALSE, scree = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc_1], PC2 = pca$x[, pc_2], group = group,
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc_1:pc_2]
    return(d)
  }
  if (scree) {
    xx<- barplot(round(percentVar, digits = 2)*100, names.arg=c(1:length(percentVar)),xlab="PC",ylab="% Variance",ylim=c(0,100), main="Scree Plot")
    text(x = xx, y = round(percentVar, digits = 4)*100, label = round(percentVar, digits = 4)*100, pos = 3, cex = 0.8, col = "black")
  }
  else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color="condition", shape="type", label = "sample")) + scale_shape_manual(values=seq(0,127)) + geom_point(size = 3) + xlab(paste0("PC",pc_1,": ", round(percentVar[pc_1] * 100, digits = 2), "% variance")) + ylab(paste0("PC",pc_2,": ", round(percentVar[pc_2] * 100, digits=2), "% variance")) + coord_fixed() + geom_text_repel(size=3)
  }
}

#needed for scatter pairs plots
ggpairs_ext <- function(data, mapping, pts=list(), smt=list(), ...){
  ggplot(data = data, mapping = mapping, ...) +
    do.call(geom_point, pts) +
    do.call(geom_smooth, smt)
}


# Function to estimate size factors (typically run for control oligos only)
estimate_control_size_factors <- function ( countData = countData, colData = colData, design = design, minRowSum = 10, ref = NULL ) {
  print( "Estimating control size factors...")
  dds <- DESeqDataSetFromMatrix( countData = countData, 
                                 colData = colData, 
                                 design = as.formula( design ) )
  
  dds <- dds[ rowSums( counts( dds ) ) > minRowSum, ]
  
  if ( ! is.null( ref ) ) {
    dds$condition <- relevel( dds$condition, ref = ref )
  }
  
  control_size_factors <- sizeFactors( estimateSizeFactors( dds ) )
  
  return( control_size_factors )
}

#function to add text to column names
appendDataFrameColumns<-function(df, prefix='', suffix='', sep='')
{
  colnames(df) <- paste(prefix, colnames(df), suffix, sep=sep)
  
  return(df)
}


#########################################################################
#OUTPUT FUNCTIONS END
#########################################################################


#the main deseq funtion - setwd() to where you want all the files outputted to
run.discrete.deseq<-function(exon_assay, sg_assay){
  x_normalization_counts = paste0("E",exon_assay,"_","SG",sg_assay,"_normalization_counts_MATRIX")
  x_filtered_counts = paste0("E",exon_assay,"_","SG",sg_assay,"_master_annotated_filtered_counts_MATRIX")
  x_normalization_counts=get(x_normalization_counts)
  x_filtered_counts =get(x_filtered_counts )
  # Get control size factors
  SF <- estimate_control_size_factors(  countData = x_normalization_counts,
                                        colData = annot_a,
                                        ref = "D4",
                                        design = "~ condition")
  #Prep DESEQ
  dds <- DESeqDataSetFromMatrix(countData = x_filtered_counts, colData = annot_a, design = ~condition)
  dds$condition <- factor(dds$condition, levels=c("D4", "D7","D14"))
  dds$condition <- relevel(dds$condition, ref = "D4")
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) <- SF
  #Run Deseq
  dds <- DESeq(dds)
  rld <- rlog(dds)
  res<-results(dds)
  res<-as.data.frame(res)
  z_score <- assay(rld) %>% as.matrix() %>% t() %>% scale() %>% t() %>% as.data.frame()
  colnames(z_score) <- paste0(colnames(z_score), "_z_score")
  z_score$mseq <- row.names(z_score)  
  row.names(z_score) <- NULL
  z_score <-  z_score %>% select(mseq, everything())
  table_wald <- degComps(dds, combs = "condition", contrast = list("condition_D7_vs_D4", "condition_D14_vs_D4"), alpha = 0.05, skip = FALSE, type = "apeglm", pairs = FALSE, fdr = "default")
  
  #Summary Table
  #production of summary table
  summary <- purrr::reduce(c(deg(table_wald[[1]], "shrunken") %>% appendDataFrameColumns(suffix="_D4_D7") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[2]], "shrunken") %>% appendDataFrameColumns(suffix="_D4_D14") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list()), left_join, by="mseq") %>% select(-c("baseMean_D4_D14")) %>% dplyr::rename(baseMean = "baseMean_D4_D7") %>% data.frame()
  #bind the z_scores to the summary table
  summary<-merge(summary, z_score, by="mseq", all.x=TRUE)
  
  #OUTPUT PLOTS AND TABLES - MAKE SURE TO RUN HK FUNCTIONS BEFORE THIS STEP
  #this is a way to get the exon and guide name from the dataframe
  #ident<-deparse(substitute(x_filtered_counts))
  #ident<-sub("_master_annotated_filtered_counts_MATRIX", "", ident)
  ident<-paste0("E",exon_assay,"_","SG",sg_assay)
  #SCREE
  pdf(paste0(ident,"_scree.pdf"))
  plotPCA.hk(rld,intgroup=c("condition", "type"), returnData=FALSE,pc_1=1, pc_2=2, scree=TRUE)
  dev.off()
  #DISPERSION
  pdf(paste0(ident,"_dispersion.pdf"))
  plotDispEsts(dds, ylim =c(1e-4,2e1))
  dev.off()
  #HEATMAP
  sampleDistMatrix <- as.matrix( dist( t( assay(rld) ) ) )
  pdf(paste0(ident,"_heatmap.pdf"))
  pheatmap(sampleDistMatrix, trace="none", col=colorRampPalette(rev(brewer.pal(9, "Blues")) )(255), adjRow = c(1,1))
  dev.off()
  remove(sampleDistMatrix)
  ###Scatter plot, checking replicate consistency
  pdf(paste0(ident,"_scatter_matrix.pdf"), width=8, height=8)
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D4")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D7")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D14")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  dev.off()
  
  #PCA PLOTS
  #PCA PLOTS - you are here save it 
  pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()+
    theme_classic()
  ggsave(paste0(ident,"_PCA.pdf"), pca, height=6, width=8)
  
  #Export rld
  write.table(as.data.frame(assay(rld)), sep="\t",file=paste0(ident,"_rld.txt"), col.names=NA)
  #Export normalized read count
  write.table(counts(dds,normalized=TRUE), sep="\t",file=paste0(ident,"_norm_count.txt"), col.names=NA)
  #Export size factor
  write.table(dds@assays@data@listData %>% as.data.frame(),sep="\t",file=paste0(ident,"_normalization_table.txt"), col.names=NA)
  #Export full table
  write.table(dds@rowRanges@elementMetadata@listData %>% as.data.frame() ,sep="\t",file=paste0(ident,"_disper_table.txt"), col.names=NA)
  
  #continuous DESEQ to add SGE RATE 
  # Get control size factors
  SF <- estimate_control_size_factors(  countData = x_normalization_counts,
                                        colData = annot_a_continuous,
                                        design = "~ condition")
  #Prep DESEQ
  dds <- DESeqDataSetFromMatrix(countData = x_filtered_counts, colData = annot_a_continuous, design = ~condition)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) <- SF
  #Run Deseq
  dds <- DESeq(dds)
  table_wald <- degComps(dds, combs = "condition", alpha = 0.05, skip = FALSE, type = "apeglm", pairs = FALSE, fdr = "default")
  rate <- as.data.frame(table_wald[[2]])%>%rownames_to_column(var="mseq")%>%appendDataFrameColumns(suffix="_continuous")
  colnames(`rate`)[colnames(`rate`) == "mseq_continuous"] <- "mseq"
  summary<-merge(summary, rate, by="mseq", all.x=TRUE)
  
  #Substract the LFC median of the synonymous_variant and inron_variant
  adj<- as.data.frame(x_normalization_counts)
  adj<-adj %>% mutate(mseq=rownames(adj)) %>% select(mseq, everything()) %>% remove_rownames() %>% select(1)
  adj <- adj %>% left_join(summary %>%select("mseq","log2FoldChange_D4_D7","log2FoldChange_D4_D14","log2FoldChange_continuous")) %>% summarise(median_D4_D7=median(log2FoldChange_D4_D7),median_D4_D14=median(log2FoldChange_D4_D14), median_continuous=median(log2FoldChange_continuous))
  median_scaled<- summary %>% mutate(median_D4_D7=adj$median_D4_D7) %>% mutate(median_D4_D14=adj$median_D4_D14) %>% mutate(median_continuous=adj$median_continuous) %>% mutate(adj_lfc_D4_D7=log2FoldChange_D4_D7-median_D4_D7) %>% mutate(adj_lfc_D4_D14=log2FoldChange_D4_D14-median_D4_D14) %>% mutate(adj_lfc_continuous=log2FoldChange_continuous-median_continuous) %>% mutate(adj_score_D4_D7=adj_lfc_D4_D7/lfcSE_D4_D7) %>% mutate(adj_score_D4_D14=adj_lfc_D4_D14/lfcSE_D4_D14) %>% mutate(adj_score_continuous=adj_lfc_continuous/lfcSE_continuous)                              
  median_scaled<-median_scaled %>% select("mseq", "adj_lfc_D4_D7", "adj_lfc_D4_D14", "adj_lfc_continuous", "adj_score_D4_D7", "adj_score_D4_D14", "adj_score_continuous")
  
  #recalculate p values based on z score produced from median scaling, standard error does not need to be recalculated after the scaling as this is a simple translation (ie moving the y axis up and down)
  median_scaled<-median_scaled %>% mutate(uncombined_two_tailed_p_D4_D7= pnorm(abs(adj_score_D4_D7),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D7 = p.adjust(uncombined_two_tailed_p_D4_D7, method = "BH")) %>% 
    mutate(uncombined_two_tailed_p_D4_D14= pnorm(abs(adj_score_D4_D14),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D14 = p.adjust(uncombined_two_tailed_p_D4_D14, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_continuous= pnorm(abs(adj_score_continuous),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_continuous = p.adjust(uncombined_two_tailed_p_continuous, method = "BH"))
  summary<-merge(summary, median_scaled, by="mseq", all.x=TRUE)
  #merge summary to annotation file - this is the final output file
  annotation_file <-paste0(ident,"_master_annotated")
  annotation_file <-get(annotation_file)
  OUT<-summary
  OUT$SG<-sg_assay
  OUT$IDENT<-ident
  OUT$EXON_GROUP<-exon_assay
  OUT<-merge(OUT, annotation_file[ , c("id","mseq","pam_mut_sgrna_id")],  by="mseq", all.x=TRUE)
  OUT$combined_name<-paste0(OUT$id,"_",OUT$EXON_GROUP)
  OUT<-OUT %>% relocate(SG, .after = mseq) %>% relocate(IDENT, .after = SG) %>% relocate(EXON_GROUP, .after = IDENT) %>% relocate(combined_name, .after = mseq) %>% relocate(pam_mut_sgrna_id, .after = EXON_GROUP) %>% relocate(id, .after = mseq)
  #make multiple outputs - one full like before and one input for the combination
  OUT_full_anot<-merge(OUT, annotation_file,  by="id", all.x=TRUE)
  names(OUT_full_anot)[names(OUT_full_anot) == 'pam_mut_sgrna_id.x'] <- 'pam_mut_sgrna_id'
  names(OUT_full_anot)[names(OUT_full_anot) == 'mseq.x'] <- 'mseq'
  OUT_full_anot$pam_mut_sgrna_id.y<-NULL
  OUT_full_anot$mseq.y<-NULL
  write.csv(OUT, file=(paste0("./",ident,"_OUT.csv")), row.names = FALSE)
  write.csv(OUT_full_anot, file=(paste0("./",ident,"_OUT_full_anot.csv")), row.names = FALSE)
  #name and output files to global environment
  #deseqname<-deparse(substitute(x_filtered_counts))
  deseqname<-paste0(ident,"_OUT")
  deseqname_full<-paste0(ident,"_OUT_full_anot")
  assign(deseqname, OUT, envir = .GlobalEnv)
  assign(deseqname_full, OUT_full_anot, envir = .GlobalEnv)
}


#########################################################################
#FUNCTIONS_DESEQ_END
#########################################################################

#RUN ALL OF THE ANALYSIS AND QC OUTPUT !!!!!!!!!!!!!!!!! ######START############################
#Function to apply to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(xx) {run.discrete.deseq(exon_assay=xx,sg_assay="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(xx) {run.discrete.deseq(exon_assay=xx,sg_assay="B")
})

###############----------------------------------------------------###########################################################
###############-----------DESEQ2_END-------------------------------###########################################################
###############----------------------------------------------------###########################################################




###############----------------------------------------------------###########################################################
###############-----------library tiling calculations-START -------###########################################################
###############----------------------------------------------------###########################################################


#combine and count 

#the above DESeq2 scripts use unique mseq per targeton - mseqs will be duplicated if not made unique as. This is because different 'id' are not unique for sequence 'mseq'.
#mutator function used by valiant is part of the 'id' nomenclature, so the same variant at the 'mseq' level is produced by varied means and will produce a duplicated mseq at the annotation stage
#the mseqs are expanded again at the final annotation stage in the deseq process and all ids will again be present in the below bindings

#tiled libraries
total_single_rad51c <- do.call("rbind.fill", list(E1_2_SGA_OUT,
                                                  E4_SGB_OUT,
                                                  E6_SGA_OUT,
                                                  E7_SGA_OUT,
                                                  E8_SGA_OUT,
                                                  E9_SGA_OUT))
total_single_rad51c$process<-"non_tiled"
#single libraries
total_tiled_rad51c <- do.call("rbind.fill", list(E2_1_SGA_OUT, 
                                                 E2_2_SGA_OUT,
                                                 E3_1_SGA_OUT, 
                                                 E3_2_SGB_OUT,
                                                 E5_1_SGA_OUT,
                                                 E5_2_SGB_OUT))
total_tiled_rad51c$process<-"tiled"                                                    
#can bind together libraries - the 'id' field produced by valiant in these dataframes are guide library agnostic (ie. the duplicate variant will will have the same id)
total_rad51c <- do.call("rbind.fill", list(total_single_rad51c,total_tiled_rad51c))

#this produces a data frame with all variants processed - not duplicated for 'id' but will be duplicated for mseq
total_unique_id_count<-total_rad51c[!duplicated(total_rad51c$id), ]
#this produces a data frame with all variants processed - not duplicated for 'mseq' - what the experimental flasks would have contained as collated unique libraries
total_unique_id_count<-total_rad51c[!duplicated(total_rad51c$mseq), ]

################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
#in this step we take only the libraries where there is a possibility of library duplication (ie. the same variant assayed with guide library a and library b)
#we group by the exon group in which the variant is found - so there will be duplicated ids as a whole
#fields of interest for each duplicated variant are put into a new field and separated by comma then key fields for calculation are ungrouped by comma into new fields
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
tiled_collapse <- total_rad51c  %>% 
  group_by(id) %>%
  dplyr::summarize(Variant_duplication= n(), 
                   process=paste(process, collapse = ','),
                   Variant_Sources = paste(SG, collapse = ','),
                   Targeton_Sources = paste(EXON_GROUP, collapse = ','),
                   PAM_status = paste(pam_mut_sgrna_id, collapse = ','),
                   mseq_combined = paste(mseq, collapse = ','),
                   baseMean = paste(baseMean, collapse = ','), 
                   raw_LFC_D4_D7 = paste(log2FoldChange_D4_D7, collapse = ','), 
                   raw_LFC_D4_D14 = paste(log2FoldChange_D4_D14, collapse = ','), 
                   raw_LFC_continuous = paste(log2FoldChange_continuous, collapse = ','),
                   
                   lfcSE_D4_D7 = paste(lfcSE_D4_D7, collapse = ','), 
                   lfcSE_D4_D14 = paste(lfcSE_D4_D14, collapse = ','), 
                   lfcSE_continuous = paste(lfcSE_continuous, collapse = ','), 
                   
                   adj_lfc_D4_D7 = paste(adj_lfc_D4_D7, collapse = ','),
                   adj_lfc_D4_D14 = paste(adj_lfc_D4_D14, collapse = ','),
                   adj_lfc_continuous = paste(adj_lfc_continuous, collapse = ',')) %>%
  
  ungroup() %>%
  separate(baseMean, sep=",", c("sga_baseMean","sgb_baseMean")) %>%
  separate(mseq_combined, sep=",", c("mseq_a","mseq_b")) %>%
  separate(adj_lfc_D4_D7, sep=",", c("sga_adj_lfc_D4_D7","sgb_adj_lfc_D4_D7")) %>%
  separate(adj_lfc_D4_D14, sep=",", c("sga_adj_lfc_D4_D14","sgb_adj_lfc_D4_D14")) %>%
  separate(adj_lfc_continuous, sep=",", c("sga_adj_lfc_continuous","sgb_adj_lfc_continuous")) %>%
  
  separate(raw_LFC_D4_D7, sep=",", c("sga_raw_LFC_D4_D7","sgb_raw_LFC_D4_D7")) %>%
  separate(raw_LFC_D4_D14, sep=",", c("sga_raw_LFC_D4_D14","sgb_raw_LFC_D4_D14")) %>%
  separate(raw_LFC_continuous, sep=",", c("sga_raw_LFC_continuous","sgb_raw_LFC_continuous")) %>%
  
  separate(lfcSE_D4_D7, sep=",", c("sga_lfcSE_D4_D7","sgb_lfcSE_D4_D7")) %>%
  separate(lfcSE_D4_D14, sep=",", c("sga_lfcSE_D4_D14","sgb_lfcSE_D4_D14"))%>%
  separate(lfcSE_continuous, sep=",", c("sga_lfcSE_continuous","sgb_lfcSE_continuous")) %>%
  
  #obtain values for instances where a variant is specific to one library - basemean
  mutate(sgb_baseMean=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_baseMean, TRUE ~ sgb_baseMean)) %>%
  mutate(sga_baseMean=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_baseMean)) %>%
  #obtain values for instances where a variant is specific to one library - mseq (the oligo used to generate the variant)
  mutate(mseq_b=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ mseq_a, TRUE ~ mseq_b)) %>%
  mutate(mseq_a=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ mseq_a)) %>%
  
  #obtain values for instances where a variant is specific to one library - median scaled LFCs (that is the median LFC of intron and synonymous variants for that targeton subtracted post-deseq2)
  #these values will not be used to calulate a weighted mean as there is only one variant instance, there will be no weight added due to error and the median scaled values will be progressed
  mutate(sgb_adj_lfc_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_D4_D7, TRUE ~ sgb_adj_lfc_D4_D7)) %>%
  mutate(sga_adj_lfc_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_D4_D7)) %>%
  mutate(sgb_adj_lfc_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_D4_D14, TRUE ~ sgb_adj_lfc_D4_D14)) %>%
  mutate(sga_adj_lfc_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_D4_D14)) %>%
  mutate(sgb_adj_lfc_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_continuous, TRUE ~ sgb_adj_lfc_continuous)) %>%
  mutate(sga_adj_lfc_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_continuous)) %>%
  
  #obtain values for instances where a variant is specific to one library - standard error from DESeq2 
  mutate(sgb_lfcSE_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_D4_D7, TRUE ~ sgb_lfcSE_D4_D7)) %>%
  mutate(sga_lfcSE_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_D4_D7)) %>%
  mutate(sgb_lfcSE_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_D4_D14, TRUE ~ sgb_lfcSE_D4_D14)) %>%
  mutate(sga_lfcSE_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_D4_D14)) %>%
  mutate(sgb_lfcSE_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_continuous, TRUE ~ sgb_lfcSE_continuous)) %>%
  mutate(sga_lfcSE_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_continuous)) %>%
  
  #obtain values for instances where a variant is specific to one library - raw log2fold changes (before median scaling), not used in calculation but to preserve provenance of value
  mutate(sgb_raw_LFC_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_D4_D7, TRUE ~ sgb_raw_LFC_D4_D7)) %>%
  mutate(sga_raw_LFC_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_D4_D7)) %>%
  mutate(sgb_raw_LFC_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_D4_D14, TRUE ~ sgb_raw_LFC_D4_D14)) %>%
  mutate(sga_raw_LFC_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_D4_D14)) %>%
  mutate(sgb_raw_LFC_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_continuous, TRUE ~ sgb_raw_LFC_continuous)) %>%
  mutate(sga_raw_LFC_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_continuous)) #%>%


#make key field contents nurmeric for further calculations
tiled_collapse[, c(9:28)] <- sapply(tiled_collapse[, c(9:28)], as.numeric)

#perform calculations to obtained weighted means -if a duplicated variant is in a pam codon, the weight for that library contribution will be 0 and the alternative library will contribute 100% to the combined value
tiled_collapse2 <-tiled_collapse %>%
  mutate(weight_a_D4_D7=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_D4_D7)^2)) %>% mutate(weight_b_D4_D7= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_D4_D7)^2)) %>%
  mutate(weight_a_D4_D14=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_D4_D14)^2)) %>% mutate(weight_b_D4_D14= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_D4_D14)^2)) %>%
  mutate(weight_a_continuous=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_continuous)^2)) %>% mutate(weight_b_continuous= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_continuous)^2)) %>%
  
  mutate(weight_a_D4_D7 = case_when(str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a") ~ 1/(sga_lfcSE_D4_D7)^2, TRUE ~ weight_a_D4_D7)) %>% mutate(weight_b_D4_D7 = case_when(str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a") ~ 1/(sgb_lfcSE_D4_D7)^2, TRUE ~ weight_b_D4_D7)) %>% 
  mutate(weight_a_D4_D14 = case_when(str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a") ~ 1/(sga_lfcSE_D4_D7)^2, TRUE ~ weight_a_D4_D14)) %>% mutate(weight_b_D4_D7 = case_when(str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a") ~ 1/(sgb_lfcSE_D4_D7)^2, TRUE ~ weight_b_D4_D7)) %>%
  mutate(weight_a_continuous = case_when(str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a") ~ 1/(sga_lfcSE_D4_D7)^2, TRUE ~ weight_a_continuous)) %>% mutate(weight_b_D4_D7 = case_when(str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a") ~ 1/(sgb_lfcSE_D4_D7)^2, TRUE ~ weight_b_D4_D7)) %>%
  
  mutate(sum_of_weight_D4_D7=rowSums(cbind(weight_a_D4_D7,weight_b_D4_D7), na.rm=TRUE))%>%
  mutate(sum_of_weight_D4_D14=rowSums(cbind(weight_a_D4_D14,weight_b_D4_D14), na.rm=TRUE))%>%
  mutate(sum_of_weight_continuous=rowSums(cbind(weight_a_continuous,weight_b_continuous), na.rm=TRUE))%>%
  
  mutate(SE_bind_D4_D7 = (sum_of_weight_D4_D7)^(-0.5)) %>%
  mutate(SE_bind_D4_D14 = (sum_of_weight_D4_D14)^(-0.5)) %>%
  mutate(SE_bind_continuous = (sum_of_weight_continuous)^(-0.5)) %>%
  
  mutate(weighted_sga_LFC_D4_D7= weight_a_D4_D7*sga_adj_lfc_D4_D7) %>% mutate(weighted_sgb_LFC_D4_D7= weight_b_D4_D7*sgb_adj_lfc_D4_D7) %>%
  mutate(weighted_sga_LFC_D4_D14= weight_a_D4_D14*sga_adj_lfc_D4_D14) %>% mutate(weighted_sgb_LFC_D4_D14= weight_b_D4_D14*sgb_adj_lfc_D4_D14) %>%
  mutate(weighted_sga_LFC_continuous= weight_a_continuous*sga_adj_lfc_continuous) %>% mutate(weighted_sgb_LFC_continuous= weight_b_continuous*sgb_adj_lfc_continuous) %>%
  
  mutate(sum_of_weighted_LFC_D4_D7=rowSums(cbind(weighted_sga_LFC_D4_D7,weighted_sgb_LFC_D4_D7), na.rm=TRUE)) %>% mutate(combined_LFC_D4_D7=sum_of_weighted_LFC_D4_D7/sum_of_weight_D4_D7) %>%
  mutate(sum_of_weighted_LFC_D4_D14=rowSums(cbind(weighted_sga_LFC_D4_D14,weighted_sgb_LFC_D4_D14), na.rm=TRUE)) %>% mutate(combined_LFC_D4_D14=sum_of_weighted_LFC_D4_D14/sum_of_weight_D4_D14) %>%
  mutate(sum_of_weighted_LFC_continuous=rowSums(cbind(weighted_sga_LFC_continuous,weighted_sgb_LFC_continuous), na.rm=TRUE)) %>% mutate(combined_LFC_continuous=sum_of_weighted_LFC_continuous/sum_of_weight_continuous) %>%
  
  # 5_2_B guide  does not overlap so variants observed only once in a pam codon need to be kept
  mutate(SE_bind_D4_D7 = case_when(Variant_duplication=="1" & PAM_status=="sgRNA_5_2_b" ~ sgb_lfcSE_D4_D7, TRUE ~ SE_bind_D4_D7)) %>%
  mutate(SE_bind_D4_D14 = case_when(Variant_duplication=="1" & PAM_status=="sgRNA_5_2_b" ~ sgb_lfcSE_D4_D14, TRUE ~ SE_bind_D4_D14)) %>%
  mutate(SE_bind_continuous = case_when(Variant_duplication=="1" & PAM_status=="sgRNA_5_2_b" ~ sgb_lfcSE_continuous, TRUE ~ SE_bind_continuous)) %>%
  
  mutate(combined_LFC_D4_D7 = case_when(Variant_duplication=="1" & PAM_status=="sgRNA_5_2_b" ~ sgb_adj_lfc_D4_D7, TRUE ~ combined_LFC_D4_D7)) %>%
  mutate(combined_LFC_D4_D14 = case_when(Variant_duplication=="1" & PAM_status=="sgRNA_5_2_b" ~ sgb_adj_lfc_D4_D14, TRUE ~ combined_LFC_D4_D14)) %>%
  mutate(combined_LFC_continuous = case_when(Variant_duplication=="1" & PAM_status=="sgRNA_5_2_b" ~ sgb_adj_lfc_continuous, TRUE ~ combined_LFC_continuous)) %>%
  
  # non-tiled data do not overlap so variants observed only once in a pam codon need to be kept
  mutate(SE_bind_D4_D7 = case_when(Variant_Sources=="A" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sga_lfcSE_D4_D7, TRUE ~ SE_bind_D4_D7)) %>%
  mutate(SE_bind_D4_D14 = case_when(Variant_Sources=="A" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sga_lfcSE_D4_D14, TRUE ~ SE_bind_D4_D14)) %>%
  mutate(SE_bind_continuous = case_when(Variant_Sources=="A" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sga_lfcSE_continuous, TRUE ~ SE_bind_continuous)) %>%
  
  mutate(SE_bind_D4_D7 = case_when(Variant_Sources=="B" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sgb_lfcSE_D4_D7, TRUE ~ SE_bind_D4_D7)) %>%
  mutate(SE_bind_D4_D14 = case_when(Variant_Sources=="B" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sgb_lfcSE_D4_D14, TRUE ~ SE_bind_D4_D14)) %>%
  mutate(SE_bind_continuous = case_when(Variant_Sources=="B" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sgb_lfcSE_continuous, TRUE ~ SE_bind_continuous)) %>%
  
  mutate(combined_LFC_D4_D7 = case_when(Variant_Sources=="A" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sga_adj_lfc_D4_D7, TRUE ~ combined_LFC_D4_D7)) %>%
  mutate(combined_LFC_D4_D14 = case_when(Variant_Sources=="A" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sga_adj_lfc_D4_D14, TRUE ~ combined_LFC_D4_D14)) %>%
  mutate(combined_LFC_continuous = case_when(Variant_Sources=="A" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sga_adj_lfc_continuous, TRUE ~ combined_LFC_continuous)) %>%
  
  mutate(combined_LFC_D4_D7 = case_when(Variant_Sources=="B" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sgb_adj_lfc_D4_D7, TRUE ~ combined_LFC_D4_D7)) %>%
  mutate(combined_LFC_D4_D14 = case_when(Variant_Sources=="B" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sgb_adj_lfc_D4_D14, TRUE ~ combined_LFC_D4_D14)) %>%
  mutate(combined_LFC_continuous = case_when(Variant_Sources=="B" & process=="non_tiled" & str_detect(PAM_status,"sgRNA_") ~ sgb_adj_lfc_continuous, TRUE ~ combined_LFC_continuous)) %>%
  
  
  #check_2<-tiled_collapse2 %>% filter(Targeton_Sources %in% "2_1" | Targeton_Sources %in% "2_2"| Targeton_Sources %in% "2_1,2_2")
  #check_2_step1<-tiled_collapse %>% filter(Targeton_Sources %in% "2_1" | Targeton_Sources %in% "2_2"| Targeton_Sources %in% "2_1,2_2")
  
  
  #remove variants that are only present in one targeton (even though covered by an overlapping tile) and are in a pam codon - ie. variants that are artifacts of fixed pam site edit
  filter(combined_LFC_D4_D7 != "NaN") %>%
  filter(combined_LFC_D4_D14 != "NaN") %>%
  filter(combined_LFC_continuous != "NaN")


###############----------------------------------------------------###########################################################
###############-----------library tiling calculations-END ---------###########################################################
###############----------------------------------------------------###########################################################


###############----------------------------------------------------###########################################################
###############-------filtering and BH FDR calculations-START------###########################################################
###############----------------------------------------------------###########################################################

#First step: we filter the inherited dataset (the output from the above processes) to obtain variants that are not NA in mseq_a 
#(some will be because they are from the unconbined source and from library b)
is_not_NA_mseq_a<- tiled_collapse2 %>% filter(!mseq_a %in% "NA") %>% filter(!is.na(mseq_a))
#group by these non-NA mseq_a sequences to make them unique
#we also create two new fields dup_mseqa, which is the number of instances that mseq_a is seen
#also dup_PAM_status_mseqa - this is similar to arrayed_mseq_a seen for the combination of guides, except here is captures and vectorizes each pam_status for each mseqa
#this is important because some mseq_a duplicates will not be listed as having a pam status and other identical sequences will be, particular to 1bp deletions of repetitive sequence
is_not_NA_mseq_a_unique <- is_not_NA_mseq_a  %>% 
  group_by(mseq_a) %>%
  dplyr::summarize(dup_mseqa= n(), 
                   dup_PAM_status_mseqa = paste(PAM_status, collapse = ','))
#join back to inherited dataframe to get the data
is_not_NA_mseq_a_unique <- inner_join(is_not_NA_mseq_a_unique, tiled_collapse2, by="mseq_a", all.x=FALSE)
#now it's safe to make unique for mseq_a
is_not_NA_mseq_a_unique <- is_not_NA_mseq_a_unique[!duplicated(is_not_NA_mseq_a_unique$mseq_a), ]
#now get the NA for mseqa
is_NA_mseq_a<- tiled_collapse2 %>% filter(mseq_a %in% "NA")
#now it's safe to make unique for mseq_b
is_NA_mseq_a <- is_NA_mseq_a[!duplicated(is_NA_mseq_a$mseq_b), ]
#fill the new fields appropriately 
is_NA_mseq_a$dup_mseqa<-NA
is_NA_mseq_a$dup_PAM_status_mseqa<-NA
#combine non-NA unique for mseq_a and NA mseq_b
total_mseq_a_unique<-rbind(is_not_NA_mseq_a_unique,is_NA_mseq_a)
#we cannot make unique for mseq_b and mseq_a together as some rows will be lost
#RENAME the same dataframe 
total_mseq_unique<-total_mseq_a_unique

#now that mseq is not duplicated can accurately calculate p-values and FDR 
total_mseq_unique<- total_mseq_unique %>%
  mutate(combined_Z_D4_D7=combined_LFC_D4_D7/SE_bind_D4_D7) %>% mutate(two_tailed_p_D4_D7= pnorm(abs(combined_Z_D4_D7),lower.tail = FALSE) *2) %>% mutate(BH_FDR_D4_D7 = p.adjust(two_tailed_p_D4_D7, method = "BH")) %>%
  mutate(combined_Z_D4_D14=combined_LFC_D4_D14/SE_bind_D4_D14) %>% mutate(two_tailed_p_D4_D14= pnorm(abs(combined_Z_D4_D14),lower.tail = FALSE) *2) %>% mutate(BH_FDR_D4_D14 = p.adjust(two_tailed_p_D4_D14, method = "BH")) %>%
  mutate(combined_Z_continuous=combined_LFC_continuous/SE_bind_continuous) %>% mutate(two_tailed_p_continuous= pnorm(abs(combined_Z_continuous),lower.tail = FALSE) *2) %>% mutate(BH_FDR_continuous = p.adjust(two_tailed_p_continuous, method = "BH"))

write.csv(total_mseq_unique,"./total_mseq_unique_processed_dataframe.csv", row.names = FALSE)

###############----------------------------------------------------###########################################################
###############-------filtering and BH FDR calculations-END--------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############------QC plots for library combinations-START-------###########################################################
###############----------------------------------------------------###########################################################

#simplify the PAM_status field so that each specific guide is not colored, but each guide per targeton (ie. not 'sgRNA_1_A' but 'A') 
library(ggpubr)
sga_sgb_combined_guide<-total_mseq_unique %>% mutate(PAM_status_simplified=case_when(str_detect(PAM_status,"_a")~ "A", str_detect(PAM_status,"_b")~ "B", str_detect(PAM_status,"")~ "-", str_detect(PAM_status,",")~ "-"))

#opens a new PDF in the working directory and pastes each plot comparison - only shows those targetons that have a posibility of duplication, some guides will only be present in library a or library b - filered by variant duplicatio field
pdf("correlated_adj_LFC_sga_sgb_combined2.pdf")
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="sgb_adj_lfc_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="combined_LFC_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="combined_LFC_D4_D7",y="sgb_adj_lfc_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="sgb_adj_lfc_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="combined_LFC_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="combined_LFC_D4_D14",y="sgb_adj_lfc_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sga_adj_lfc_continuous = as.numeric(sga_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="combined_LFC_continuous",y="sga_adj_lfc_continuous", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sgb_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sgb_adj_lfc_continuous = as.numeric(sgb_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="sgb_adj_lfc_continuous",y="combined_LFC_continuous", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
dev.off()


#same plot as above but coloured by exon according to VEP annotation
pdf("correlated_adj_LFC_sga_sgb_combined2_EXON.pdf")
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="sgb_adj_lfc_D4_D7", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n")) %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="combined_LFC_D4_D7", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n")) %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="combined_LFC_D4_D7",y="sgb_adj_lfc_D4_D7", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="sgb_adj_lfc_D4_D14", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="combined_LFC_D4_D14", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="combined_LFC_D4_D14",y="sgb_adj_lfc_D4_D14", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sga_adj_lfc_continuous = as.numeric(sga_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="combined_LFC_continuous",y="sga_adj_lfc_continuous", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sgb_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sgb_adj_lfc_continuous = as.numeric(sgb_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="sgb_adj_lfc_continuous",y="combined_LFC_continuous", color = "Targeton_Sources", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
dev.off()

###############----------------------------------------------------###########################################################
###############------QC plots for library combinations-END---------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############------QC plots for editing-START--------------------###########################################################
###############----------------------------------------------------###########################################################

#positional effect QC generation 
position.plas.plot<-function(exon,sg){
  locus_input <- paste0("E",exon,"_SG",sg,"_OUT_full_anot")
  locus_input <- get(locus_input)
  locus_input_plasmid<-locus_input
  locus_input_plasmid$D4_MEAN <-(locus_input_plasmid$D4R1 + locus_input_plasmid$D4R2 + locus_input_plasmid$D4R3)/3
  locus_input_plasmid$D4_DIV_PLAS<-locus_input_plasmid$D4_MEAN/locus_input_plasmid$PLASMID_MEAN
  
  locus_input_plasmid$D4R1_DIV_PLAS<-locus_input_plasmid$D4R1/locus_input_plasmid$PLASMID_MEAN
  locus_input_plasmid$D4R2_DIV_PLAS<-locus_input_plasmid$D4R2/locus_input_plasmid$PLASMID_MEAN
  locus_input_plasmid$D4R3_DIV_PLAS<-locus_input_plasmid$D4R3/locus_input_plasmid$PLASMID_MEAN
  
  pam_df<-locus_input_plasmid %>% filter(!is.na(pam_mut_sgrna_id))
  pam_df<-pam_df %>% filter(!pam_mut_sgrna_id %in% "")
  locus_input_plasmid$mut_type<-as.factor(locus_input_plasmid$mut_type)
  levels(locus_input_plasmid$mut_type)[levels(locus_input_plasmid$mut_type)==''] <- 'other'
  
  pos_plot_mean<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate Mean")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("Average Day_4 SGE Counts / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_mean)
  
  pos_plot_D4R1<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4R1_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4R1_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate 1")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("D4R1 SGE Counts / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_D4R1)
  
  pos_plot_D4R2<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4R2_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4R2_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate 2")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("D4R2 SGE Counts / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_D4R2)
  
  pos_plot_D4R3<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4R3_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4R3_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate 3")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("D4R3 / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_D4R3)
  
  p<-plot_grid(pos_plot_mean, pos_plot_D4R1, pos_plot_D4R2, pos_plot_D4R3, nrow = 1, ncol = 4, label_size = 12)
  
  pdf(paste0("E",exon,"_SG",sg,"_positional_effect_D4_DIV_PLASMID.pdf"), width=16, height=4, pointsize = 1)
  
  print(p)
  
  dev.off() 
}

#Function to apply to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(xxx) {position.plas.plot(exon=xxx,sg="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(xxx) {position.plas.plot(exon=xxx,sg="B")
})

###############----------------------------------------------------###########################################################
###############------QC plots for editing-END----------------------###########################################################
###############----------------------------------------------------###########################################################


###############----------------------------------------------------###########################################################
###############---Expansion of data frame to full annotation-START-###########################################################
###############----------------------------------------------------###########################################################

#annotation and re-expansion to full dataset
E1_2_SGA_META_VEP<-merge(chr17_58692607_58692852_plus_sgRNA_1_2_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58692607_58692852_plus_sgRNA_1_2_a.tsv, by="id", all.x=TRUE)
E2_1_SGA_META_VEP<-merge(chr17_58694822_58695075_plus_sgRNA_2_1_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58694822_58695075_plus_sgRNA_2_1_a.tsv, by="id", all.x=TRUE)
E2_2_SGA_META_VEP<-merge(chr17_58694985_58695224_plus_sgRNA_2_2_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58694985_58695224_plus_sgRNA_2_2_a.tsv, by="id", all.x=TRUE)
E3_1_SGA_META_VEP<-merge(chr17_58696626_58696865_plus_sgRNA_3_1_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58696626_58696865_plus_sgRNA_3_1_a.tsv, by="id", all.x=TRUE)
E5_1_SGA_META_VEP<-merge(chr17_58709712_58709961_plus_sgRNA_5_1_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58709712_58709961_plus_sgRNA_5_1_a.tsv, by="id", all.x=TRUE)
E6_SGA_META_VEP<-merge(chr17_58720663_58720902_plus_sgRNA_6_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58720663_58720902_plus_sgRNA_6_a.tsv, by="id", all.x=TRUE)
E7_SGA_META_VEP<-merge(chr17_58723963_58724206_plus_sgRNA_7_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58723963_58724206_plus_sgRNA_7_a.tsv, by="id", all.x=TRUE)
E8_SGA_META_VEP<-merge(chr17_58732410_58732651_plus_sgRNA_8_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58732410_58732651_plus_sgRNA_8_a.tsv, by="id", all.x=TRUE)
E9_SGA_META_VEP<-merge(chr17_58734062_58734302_plus_sgRNA_9_a_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58734062_58734302_plus_sgRNA_9_a.tsv, by="id", all.x=TRUE)

#annotation and re-expansion to full dataset
E3_2_SGB_META_VEP<-merge(chr17_58696699_58696943_plus_sgRNA_3_2_b_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58696699_58696943_plus_sgRNA_3_2_b.tsv, by="id", all.x=TRUE)
E4_SGB_META_VEP<-merge(chr17_58703145_58703389_plus_sgRNA_4_b_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58703145_58703389_plus_sgRNA_4_b.tsv, by="id", all.x=TRUE)
E5_2_SGB_META_VEP<-merge(chr17_58709875_58710114_plus_sgRNA_5_2_b_meta.csv,VEP_OUTPUT_VEP_OUTPUT_chr17_58709875_58710114_plus_sgRNA_5_2_b.tsv, by="id", all.x=TRUE)

#FUNCTION to add an additional field which is ID and SG combined - to label a variant annotation as either from SGA or SGB source - which is needed subsequently 
#label annotations with unique id
label <- function(exon, sg){
  input<-paste0("E",exon,"_SG",sg,"_META_VEP")
  input<-get(input)
  input$SG<-sg
  input$exon_group<-exon
  input$required_annotation_id<-paste0(input$id,"_",input$SG)
  countdfname<-paste0("E",exon,"_SG",sg,"_META_VEP")
  assign(countdfname, input, envir = .GlobalEnv)
}

#Function to apply to listed exons and designated sgRNA_A
lapply (c("1_2","2_1","2_2","3_1","5_1","6","7","8","9"), function(xx) {label(exon=xx,sg="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c("3_2","4","5_2"), function(xx) {label(exon=xx,sg="B")
})


total_META_VEP <- do.call("rbind.fill", list(E1_2_SGA_META_VEP,
                                             E2_1_SGA_META_VEP,
                                             E2_2_SGA_META_VEP,
                                             E3_1_SGA_META_VEP,
                                             E5_1_SGA_META_VEP,
                                             E6_SGA_META_VEP,
                                             E7_SGA_META_VEP,
                                             E8_SGA_META_VEP,
                                             E9_SGA_META_VEP,
                                             E3_2_SGB_META_VEP,
                                             E4_SGB_META_VEP,
                                             E5_2_SGB_META_VEP))
#remove any variant annotations that are variants in the constant regions 
total_META_VEP<- total_META_VEP %>% filter(vcf_var_in_const==0 | is.na(vcf_var_in_const))
total_META_VEP$pam_mut_sgrna_id<-NULL 
total_META_VEP$pam_mut_annot<-NULL 
total_META_VEP$exon_group<-NULL
total_META_VEP$ref_start<-NULL
total_META_VEP$ref_end<-NULL
total_META_VEP$pam_seq<-NULL
total_META_VEP$ref_seq<-NULL
total_META_VEP$mseq_no_adapt<-NULL
total_META_VEP$oligo_length<-NULL
total_META_VEP$vcf_var_in_const<-NULL


#take the inherited dataset (post combination of guides and combination of overlapping targetons)
#if the variant is observed once (ie is not from combined from datasets, or is only observed in lirbary a or library b) then the required source of annoration
#will be whichever library is the source of the data (either A or B)

total_mseq_unique_min<- total_mseq_unique[ ,c(1:10,40:42,50,52,54,55:63)]

annotation_processed_single<-total_mseq_unique_min%>% filter(Variant_duplication %in% "1") %>%
  mutate(required_annotation_source=case_when(str_detect(Variant_Sources,"A")~ "A", TRUE~ "B"))

#if the variant is derived from two libraries (A and B) then take metadata from source B unless the arrayed pam status contains a B guide, then take A
annotation_processed_duplicate<-total_mseq_unique_min %>% filter(Variant_duplication %in% "2") %>%
  mutate(required_annotation_source=case_when(str_detect(dup_PAM_status_mseqa,"_a") | str_detect(PAM_status,"_a")~ "B", TRUE~ "A")) %>% 
  mutate(required_annotation_source=case_when(str_detect(dup_PAM_status_mseqa,"sgRNA_2_1_a") | str_detect(PAM_status,"sgRNA_2_1_a,sgRNA_2_2_a")~ "A", TRUE~ required_annotation_source))

#BIND together these two dataframes
req_anno<-do.call("rbind.fill", list(annotation_processed_single, annotation_processed_duplicate))
#CREATE a new required annotation field in the dataframe
req_anno$required_annotation_id<-paste0(req_anno$id,"_",req_anno$required_annotation_source)
req_anno$mseq_a<-NULL
req_anno$mseq_b<-NULL

#create minimal fields needed to link mseq to all (mseq redundant ids)
minimal_annotation<-total_META_VEP[ ,c("required_annotation_id","mseq")]
#REQ_ANNO has one id per unique mseq - need to get all ids per unique mseq
#1.fish for the required mseq using the required_annotation_id
annotation_key_association<-merge(minimal_annotation, req_anno, by="required_annotation_id", all.x=FALSE)
annotation_key_association$required_annotation_id<-NULL
annotation_key_association$id<-NULL
#2.get all the ids associated with the appropriate mseqs by merging the data by mseq
annotation_full<-merge(total_META_VEP, annotation_key_association,  by="mseq", all.x=FALSE)
#check which ids are missing from the final dataset versus the metadata
check<-subset(total_unique_id_count, !(id %in% annotation_full$id))
#all of the variants that are missing are likely to be pam site edits that further modify the mseq to be something unique, that would not occur in wild-type mutation 
#most will also be duplciated due to merging with expanded metadata annotation
duplicated_id <- annotation_full  %>% 
  group_by(id) %>%
  filter(n() > 1)

####DUPLICATE ACCESSION PROCESSING - START 
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING

#There are some mseqs associated with an ID that we don't want because they are an meseq that is a duplicate of a mseq that is annotated as having a pam edit, but isn't itself annotated as a pam containing mseq
#this happens with 1del deletions, and in some cases 2del if the pam site is in a 2 del location 
#some ids are represented by multiple mseqs 
#these mseqs will have remained distinct in the annotation process
#for library a and library b the required annotation id will have pulled the correct mseq(s)
#some required annotation ids will have pulled mseqs that are duplications

##remove accessions that are duplicated for id where the annotation states that is not observed twice in the dataset - this is an artifact of merging mseqs and then pulling by mseqs again

##can delete these - start 
#annotation_full$variant_n_instances_between_libraries<-as.character(annotation_full$variant_n_instances_between_libraries)
#annotation_full$variant_n_instances_between_targetons<-as.character(annotation_full$variant_n_instances_between_targetons)

#filtered_annotation_full <- annotation_full[!(duplicated(annotation_full$id) & (annotation_full$variant_n_instances_between_libraries == 1)), ]
#filtered_annotation_full <- filtered_annotation_full[!(duplicated(filtered_annotation_full$id) & is.na(filtered_annotation_full$variant_n_instances_between_libraries)), ]
##can delete these - end

#create a column that counts the number of duplicates (1=not duplciated)
filtered_annotation_full<-annotation_full %>%
  group_by(id) %>% dplyr::mutate(duplicated_id=n())
#annotated this new dataframe with a column to state if the mseq contains pam site annotation, duplicates which don't contain this annotation are an artifact of merging
duplicated_set<- filtered_annotation_full %>% filter(duplicated_id >1) %>% mutate(is_guide=case_when(str_detect(dup_PAM_status_mseqa,"sgRNA") | str_detect(PAM_status,"sgRNA") ~ "Y", TRUE ~ "N"))
#check it in excel
write.csv(duplicated_set,"./duplicated_set.csv", row.names = FALSE)

#annotate the datafram with a removal indication column. Remove the row if it is duplicate for id and contains no guide annotation and is not exon17 (special case) or 
#is exon17 duplication without a pam status that indicates the row is from an array of values - this is a special case becuase the GG at boarder of overlap for 17_1 and 17_2, one G is in overlapping region and one is not, which leads to the same mseq being annotated with multiple ids
duplicated_set_2<- duplicated_set  %>% mutate(remove=case_when(PAM_status=="" ~ "delete", TRUE ~ "keep"))
#remove rows based on the removal indication column
duplicated_set_3<-duplicated_set_2 %>% filter(!remove %in% 'delete')
#check it in excel
write.csv(duplicated_set_3,"./duplicated_set_3.csv", row.names = FALSE)

#remove the columns that were generated to process the filtering  
duplicated_set_3$duplicated_id<-NULL
duplicated_set_3$remove<-NULL
duplicated_set_3$is_guide<-NULL

#MAKE UNIQUE and count the number of unique ids in the final list, should be the same as 'annotation_full' > 'variant_n_instances_between_targetons==1' > make unique
filtered_duplicates<-duplicated_set_3[!duplicated(duplicated_set_3$id), ]

###creation of final dataset - post duplicate removal
#remove those accessions that are duplicated from the main dataframe
annotation_full_minus_duplicates<-annotation_full[!annotation_full$id %in% filtered_duplicates$id,]
#add back in the filtered previously duplicated accesssions
final_annotated_dataset<-do.call("rbind.fill", list(annotation_full_minus_duplicates, filtered_duplicates))

#check that no duplicate ids remain
duplicated <- final_annotated_dataset  %>% 
  group_by(id) %>%
  filter(n() > 1)
#check which variants are not in the final dataset which were in the orignal post DEseq2 list
check<-subset(total_unique_id_count, !(id %in% final_annotated_dataset$id))
#are these accessions unique or are there duplicates - shouldn't be an any duplciates
check_unique<-check[!duplicated(check$id), ]

#this is your final annotated dataset: 'final_annotated_dataset', write it to file
write.csv(final_annotated_dataset,"./final_annotated_dataset.csv", row.names = FALSE)

####DUPLICATE ACCESSION PROCESSING - END
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING


###############----------------------------------------------------###########################################################
###############---Expansion of data frame to full annotation-END---###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############---Cleanup, classification and clustering START-----###########################################################
###############----------------------------------------------------###########################################################

#simplify the pam status field and create a flag field for variants that only come from one dataset and are pam variant codon
clean_final_annotated_dataset<- final_annotated_dataset%>% 
  mutate(pam_codon=case_when(str_detect(dup_PAM_status_mseqa,"sgRNA") | str_detect(PAM_status,"sgRNA")~ "Y", TRUE ~ "N")) %>%
  mutate(pam_flag=case_when(Variant_duplication=="1" & str_detect(dup_PAM_status_mseqa,"sgRNA") | Variant_duplication=="1" & str_detect(PAM_status,"sgRNA" ) ~ "Y", TRUE ~ "N"))

#select key columns
clean_final_annotated_dataset <- clean_final_annotated_dataset[ ,c(1:9,11:90,95:114,116,117)]

#reorder the dataframe
clean_final_annotated_dataset_reordered <-clean_final_annotated_dataset[, c(1:2,90:94,110,111,95:109,3:89)]

#rename columns and simplify some fields
clean_final_annotated_dataset_reordered <- clean_final_annotated_dataset_reordered %>% 
  dplyr::rename(variant_duplication = "Variant_duplication") %>% 
  dplyr::rename(variant_source = "Variant_Sources") %>%
  dplyr::rename(targeton_source = "Targeton_Sources") %>%
  dplyr::rename(processed_adj_lfc_D4_D7 = "combined_LFC_D4_D7") %>%
  dplyr::rename(processed_adj_lfc_D4_D14 = "combined_LFC_D4_D14") %>%
  dplyr::rename(processed_adj_lfc_continuous = "combined_LFC_continuous") %>%
  
  dplyr::rename(processed_z_D4_D7 = "combined_Z_D4_D7") %>%
  dplyr::rename(processed_z_D4_D14 = "combined_Z_D4_D14") %>%
  dplyr::rename(processed_z_continuous = "combined_Z_continuous") %>% 
  mutate(process = case_when(process=="tiled,tiled" ~ "tiled", TRUE ~ process))


#add a pos_ref_new field which could be useful - different id an give same variant, but are unique
clean_final_annotated_dataset_reordered$pos_ref_new <- paste0(clean_final_annotated_dataset_reordered$mut_position,"_",clean_final_annotated_dataset_reordered$ref,"_",clean_final_annotated_dataset_reordered$new)
#re-order some more
clean_final_annotated_dataset_reordered <- clean_final_annotated_dataset_reordered %>%
  relocate(pos_ref_new, .before = variant_duplication) %>%
  relocate(processed_adj_lfc_D4_D7, processed_adj_lfc_D4_D14, processed_adj_lfc_continuous, .after =pam_flag) %>%
  relocate(processed_z_D4_D7, processed_z_D4_D14, processed_z_continuous, .after = processed_adj_lfc_continuous) %>%
  relocate(two_tailed_p_D4_D7,two_tailed_p_D4_D14,two_tailed_p_continuous, .after = processed_z_continuous) %>%
  relocate(BH_FDR_D4_D7,BH_FDR_D4_D14,BH_FDR_continuous, .after = two_tailed_p_continuous)
#write.csv(clean_final_annotated_dataset_reordered ,"./clean_rad51c.csv", row.names = FALSE)


#MAKE dataset with full expaneded annotations
test<- clean_final_annotated_dataset_reordered
#set levels to reorder by mutator - so that when duplicate variants are later removed, the duplicate accession is preserved in order or preference
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")
#re-order by mutator according to set levels
test<-test %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)
#write the df to file - this is the full expanded annotation file
write.csv(test ,"./clean_rad51c.csv", row.names = FALSE)


######CLUSTERING and CLASSIFICATION
#clustering
cluster.frame<-test 

#stratify for variants that have FDR>0.01
unchanged<-cluster.frame %>% filter(BH_FDR_D4_D14>=0.01)
enriched<-cluster.frame %>% filter(BH_FDR_D4_D14<=0.01) %>% filter(processed_z_D4_D14>0) 
depleted<-cluster.frame %>% filter(BH_FDR_D4_D14<=0.01) %>% filter(processed_z_D4_D14<0) 

cluster.frame$mut_type<-as.factor(cluster.frame$mut_type)
levels(cluster.frame$mut_type)[levels(cluster.frame$mut_type)==''] <- 'other'
cluster.frame_no_pam<-cluster.frame %>% filter(!pam_flag %in% "Y")


library(mclust)
set.seed(0)
EEV <- Mclust(depleted%>% dplyr::select(c("processed_z_D4_D14","processed_z_D4_D7")), G=2, modelNames = "EEV")
EEV_plot<-plot(EEV, what="classification", xlim=c(-40,5), ylim=c(-40,5), xlab="D4D14", ylab="D4D7", symbols=c(3,0,4,2))
dev.off()

all_depleted <- cbind(depleted, EEV$z %>% as.data.frame() %>% `colnames<-` (c("C1_z","C2_z")),EEV$classification %>% as.data.frame() %>% `colnames<-` (c("Cluster")), EEV$uncertainty %>% as.data.frame() %>% `colnames<-` (c("uncertainty")))
enriched$Cluster<-3
unchanged$Cluster<-4

rates<- do.call("rbind.fill", list(all_depleted, enriched, unchanged))
rates$Cluster<-as.factor(rates$Cluster)

levels(rates$Cluster)[levels(rates$Cluster)=='1'] <- 'fast depleted'
levels(rates$Cluster)[levels(rates$Cluster)=='2'] <- 'slow depleted'
levels(rates$Cluster)[levels(rates$Cluster)=='3'] <- 'enriched'
levels(rates$Cluster)[levels(rates$Cluster)=='4'] <- 'unchanged'

rates$Cluster <- relevel(rates$Cluster, "unchanged", "fast depleted", "slow depleted", "enriched")

#move class next to id
x_final_frame_class <- rates %>% 
  mutate(chrom_pos_ref_alt = paste0("chr17_",pos_ref_new)) %>%
  relocate(chrom_pos_ref_alt, .after = id) %>% 
  relocate(Cluster, .after = chrom_pos_ref_alt)
x_final_frame_class$pos_ref_new <- NULL


#Consequence simplification
x_final_frame_class_v2<- x_final_frame_class  %>% mutate(slim_consequence=case_when(str_detect(Consequence,"stop_gained") ~ "stop_gained")) %>%
  mutate(slim_consequence=case_when(str_detect(Consequence,"UTR_variant") ~ "UTR", TRUE ~ slim_consequence)) %>%
  mutate(slim_consequence=case_when(str_detect(Consequence,"frameshift_variant") ~ "frameshift", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(mutator,"inframe") ~ "codon_deletion", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"intron_variant") ~ "intron", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"missense_variant") ~ "missense", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"synonymous_variant") ~ "synonymous", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"splice_acceptor_variant") ~ "splice_acceptor", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"splice_donor_variant") ~ "splice_donor", TRUE ~ slim_consequence)) %>%
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"start_lost") & mutator== "snvre" |
                                      str_detect(Consequence,"start_lost") & mutator== "snv" |
                                      str_detect(Consequence,"start_lost") & mutator== "custom" |
                                      str_detect(Consequence,"start_lost") & mutator== "ala"~ "start_lost", TRUE ~ slim_consequence)) %>% 
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"stop_lost") & mutator== "snvre" |
                                      str_detect(Consequence,"stop_lost") & mutator== "snv" |
                                      str_detect(Consequence,"stop_lost") & mutator== "custom" |
                                      str_detect(Consequence,"stop_lost") & mutator== "ala"~ "stop_lost", TRUE ~ slim_consequence)) %>% 
  
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"stop_retained_variant") & mutator=="snvre" |
                                      str_detect(Consequence,"stop_retained_variant") & mutator=="snv" |
                                      str_detect(Consequence,"stop_retained_variant") & mutator=="custom" ~ "synonymous", TRUE ~ slim_consequence)) %>% 
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"inframe_deletion") & mutator=="custom" ~ "clinical_inframe_deletion", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"inframe_insertion") & mutator=="custom" ~ "clinical_inframe_insertion", TRUE ~ slim_consequence)) 

#re-order some more
x_final_frame_class_v3 <- x_final_frame_class_v2 %>%
  relocate(slim_consequence, .before = pam_flag) %>%
  relocate(pam_flag, .after = pam_codon) %>%
  relocate(slim_consequence, .after = Cluster) %>%
  relocate(processed_z_D4_D7,processed_z_D4_D14,processed_z_continuous,BH_FDR_D4_D7,BH_FDR_D4_D14,BH_FDR_continuous, .before = slim_consequence) %>%
  relocate(HGVSc,HGVSp, .after = id)

x_final_frame_class_v5 <- x_final_frame_class_v3 %>% mutate(region=case_when(targeton_source=="4" ~ "4")) %>%
  mutate(region=case_when(targeton_source=="6" ~ "6", TRUE ~ region)) %>%
  mutate(region=case_when(targeton_source=="7" ~ "7", TRUE ~ region)) %>% 
  mutate(region=case_when(targeton_source=="8" ~ "8", TRUE ~ region)) %>% 
  mutate(region=case_when(targeton_source=="9" ~ "9", TRUE ~ region)) %>%
  mutate(region=case_when(str_detect(targeton_source,"1_") ~ "1", TRUE ~ region)) %>%
  mutate(region=case_when(str_detect(targeton_source,"2_") ~ "2", TRUE ~ region)) %>%
  mutate(region=case_when(str_detect(targeton_source,"3_") ~ "3", TRUE ~ region)) %>% 
  mutate(region=case_when(str_detect(targeton_source,"5_") ~ "5", TRUE ~ region))

#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")
#re-order by mutator according to set levels
final_expanded<-x_final_frame_class_v5 %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)
#write the df to file - this is the full expanded annotation file
write.csv(final_expanded ,"./sge_rad51c_full_annotations.csv", row.names = FALSE)


#organise expanded dataframe to factor by clinical observation
databases<-final_expanded %>% mutate(is_in_clinvar=case_when(ClinVar == "-" | ClinVar_CLNSIG== "not_provided" ~ "N", TRUE ~ "Y")) %>%
  mutate(is_in_gnomad=case_when(gnomAD == "-" ~ "N", TRUE ~ "Y")) %>%
  mutate(present_accession=case_when(is_in_clinvar == "N" ~ "Unobserved", TRUE ~ "ClinVar_only")) %>%
  mutate(present_accession=case_when(is_in_clinvar == "N" & is_in_gnomad =="Y" ~ "gnomAD_only", TRUE ~present_accession)) %>%
  mutate(present_accession=case_when(is_in_clinvar == "Y" & is_in_gnomad =="Y" ~ "ClinVar_and_gnomAD", TRUE ~present_accession)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Uncertain_significance" ~ "Uncertain_significance")) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="-" ~ "Unobserved", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Pathogenic/Likely_pathogenic" ~ "Pathogenic/Likely_pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Pathogenic" ~ "Pathogenic/Likely_pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Likely_pathogenic" ~ "Pathogenic/Likely_pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Likely_benign" ~ "Benign/Likely_benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Benign/Likely_benign" ~ "Benign/Likely_benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Benign" ~ "Benign/Likely_benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="Conflicting_interpretations_of_pathogenicity" ~ "Conflicting_interpretation", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(present_accession_filter=paste0(present_accession,"_",mseq))
databases <- databases[!duplicated(databases$present_accession_filter), ]    
databases_expanded <- databases

#organise expanded dataframe to factor by clinical observation
y_levels<-c("ClinVar_only","ClinVar_and_gnomAD","gnomAD_only","Unobserved")
#re-order by mutator according to set levels
databases_unique<-databases_expanded %>% mutate(present_accession =  factor(present_accession, levels = y_levels)) %>%
  arrange(present_accession)

#remove duplicated rows by HGVSc
databases_unique <- databases_unique[!duplicated(databases_unique$HGVSc), ]
#remove rows without a HGVSc - these are variants that revert PPEs back to wild-type
databases_unique <- databases_unique %>% filter(!HGVSc %in% "-")


#organise expanded dataframe to factor by clinical observation
databases_clinvar_and_unobs <-databases_unique
databases_clinvar_and_unobs$ClinVar_CLNSIG_slim<- factor(databases_clinvar_and_unobs$ClinVar_CLNSIG_slim, levels=c("Unobserved","Benign/Likely_benign","Uncertain_significance","Conflicting_interpretation","Pathogenic/Likely_pathogenic"))
databases_clinvar_and_unobs$ClinVar_CLNSIG_slim<-revalue(databases_clinvar_and_unobs$ClinVar_CLNSIG_slim, c("Unobserved" = "Unobserved",                                                                             
                                                                                                            "Benign/Likely_benign" ="Benign/Likely benign",
                                                                                                            "Pathogenic/Likely_pathogenic"= "Pathogenic/Likely pathogenic",
                                                                                                            "Uncertain_significance" = "Uncertain significance",
                                                                                                            "Conflicting_interpretation"="Conflicting interpretation"))
#remove duplicated rows by mseq
databases_clinvar_and_unobs <- databases_clinvar_and_unobs[!duplicated(databases_clinvar_and_unobs$mseq), ]
#remove duplicated rows by HGVSc to be sure that no duplications of variants remain
databases_clinvar_and_unobs <- databases_clinvar_and_unobs[!duplicated(databases_clinvar_and_unobs$HGVSc), ]


#re-organise the data frame for easier interpretation
dataframe_cleaned_up <- databases_clinvar_and_unobs %>% relocate(region, .before = EXON) %>%
  relocate(C1_z, .after = SE_bind_continuous) %>%
  relocate(C2_z,uncertainty, .after = C1_z) %>%
  relocate(is_in_clinvar,is_in_gnomad,present_accession,ClinVar_CLNSIG_slim,present_accession_filter, .after = dbSNP) %>%
  relocate(ClinVar_CLNSIG_slim, .after = ClinVar_CLNSIG) %>%
  relocate(vcf_alias,vcf_var_id, .after = present_accession) %>%
  relocate(Consequence, .after = slim_consequence) %>%
  relocate(region,EXON,INTRON, .after = Consequence) %>%
  relocate(cDNA_position,CDS_position,Protein_position,Amino_acids,Codons, .after = INTRON) %>%
  relocate(ref_chr,ref_strand,mut_position,ref,new,ref_aa,alt_aa,mutator, .after = Codons) %>%
  select(!present_accession_filter) %>%
  select(!Existing_variation) %>%
  select(!mut_type) %>%
  select(!MOTIF_NAME) %>%
  select(!MOTIF_POS) %>%
  select(!HIGH_INF_POS) %>%
  select(!MOTIF_SCORE_CHANGE) %>%
  select(!TRANSCRIPTION_FACTORS) %>%
  select(!CLIN_SIG) %>%
  select(!SOMATIC) %>%
  select(!PHENO) %>%
  select(!DISTANCE) %>%
  select(!STRAND) %>%
  select(!FLAGS) %>%
  select(!PICK) %>%
  select(!MANE_PLUS_CLINICAL) %>%
  select(!TREMBL) %>%
  select(!SOURCE) %>%
  select(!gnomAD_AF.1) %>%
  select(!Allele) %>%
  select(!PAM_status) %>%
  select(!variant_source) %>%
  select(!Location) %>%
  select(!Gene) %>%
  select(!Feature) %>%
  dplyr::rename(functional_classification="Cluster") %>% 
  dplyr::rename(z_score_D4_D7="processed_z_D4_D7") %>%
  dplyr::rename(z_score_D4_D14="processed_z_D4_D14") %>%
  dplyr::rename(z_score_continuous="processed_z_continuous") %>%
  dplyr::rename(GMM_fast_depleted_z="C1_z") %>%
  dplyr::rename(GMM_slow_depleted_z="C2_z") %>%
  dplyr::rename(GMM_uncertainty="uncertainty")%>%
  dplyr::rename(target_region_source="targeton_source")%>%
  arrange(mut_position)
#re-label target region (targeton) 1_2 as 1
dataframe_cleaned_up$target_region_source[which(dataframe_cleaned_up$target_region_source == '1_2')] <- "1"

#get the column names 
dataframe_cleaned_up_columns<-as.data.frame(colnames(dataframe_cleaned_up))
#write the df columns  to file - this is the full expanded annotation file
write.csv(dataframe_cleaned_up_columns ,"./dataframe_columns.csv", row.names = FALSE)  

#this is the final dataframe --------
#write the df to file - this is the full expanded annotation file
write.csv(dataframe_cleaned_up ,"./sge_rad51c.csv", row.names = FALSE)  


###############----------------------------------------------------###########################################################
###############---Cleanup, classification and clustering END-------###########################################################
###############----------------------------------------------------###########################################################


##END OF CODE



