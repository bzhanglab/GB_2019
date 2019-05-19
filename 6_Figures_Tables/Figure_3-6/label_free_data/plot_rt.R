library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)


load_data=function(f){
    a <- read_tsv(f)
    filename <- a$method
    combine_or_ind = ifelse(str_detect(filename,pattern = "^combine_"),"combine","individual")
    peptide_or_psm = ifelse(str_detect(filename,pattern = "peptide_level"),"peptide","PSM")
    se = detect_search_engine(filename)
    fdr_method = detect_fdr_method(filename)
    a$level = combine_or_ind
    a$peptide_or_psm = peptide_or_psm
    a$search_engine = se
    a$FDR_method = fdr_method
    return(a)
}

detect_search_engine=function(x){
    se = ifelse(str_detect(x,pattern = "_comet_"),"Comet",
                ifelse(str_detect(x,pattern = "_msgf_"),"MS-GF+","X!Tandem"))
    return(se)
}

detect_fdr_method=function(x){
    me <- ifelse(str_detect(x,pattern = "_global_fdr"),"Global FDR",
           ifelse(str_detect(x,pattern = "_separate_fdr"),"Separate FDR","Two-stage FDR"))
    return(me)
}

wilcox_test=function(x,y){
    pp <- wilcox.test(x~y)
    return(pp$p.value)
}

wilcox_test_multiple=function(x,y){
    x1 <- x[y=="Global FDR"]
    x2 <- x[y=="Separate FDR"]
    x3 <- x[y=="Two-stage FDR"]
    pp12 <- wilcox.test(x1,x2)
    pp13 <- wilcox.test(x1,x3)
    pp23 <- wilcox.test(x2,x3)
    return(data.frame(m12=pp12$p.value,m13=pp13$p.value,m23=pp23$p.value))
}

remove_redundant=function(ind_data){
    pep_pepquery <- ind_data %>% group_by(peptide,level,peptide_or_psm,search_engine,FDR_method) %>% 
        summarise(PepQuery=ifelse(any(pepquery==1),1,0)) %>%
        ungroup()
    ind_data$pepquery <- NULL
    ind_data <- inner_join(ind_data,pep_pepquery) %>% arrange(desc(score)) %>% 
        group_by(peptide,level,peptide_or_psm,search_engine,FDR_method) %>% filter(row_number()==1)
    #ind_data$error <- ifelse(ind_data$error>10,10,ind_data$error)
    ind_data$PepQuery <- ifelse(ind_data$PepQuery==1,"Pass","Fail")
    return(ind_data)
}

plot_boxplot=function(data,mae){
    gg <- data %>% ggplot(aes(y=error,x=FDR_method,colour=PepQuery)) +
        geom_boxplot(position = position_dodge(),width=0.75,outlier.size=0.2,fatten=0.8)+
        geom_text(data=mae,aes(x=FDR_method,y=c(11),label=sprintf("%.2f",mae)),position = position_dodge(width=0.75),size=2.3)+
        facet_grid(search_engine~.)+
        xlab("FDR estimation method")+
        ylab("Absolute RT error (minute)")+
        #ylim(0,12)+
        theme_bw()+
        theme(legend.position="top",legend.text=element_text(size=7),
              legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
        theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
        ggplot2::coord_cartesian(ylim=c(0, 13))
    return(gg)
}

## input parameter
input_dir <- "autort/pred_res/"

## load data
input_files <- list.files(path = input_dir,pattern = "test.csv",full.names = TRUE,include.dirs = TRUE,recursive = TRUE)
all_data <- bind_rows(lapply(input_files, load_data))
all_data$sample <- str_replace_all(all_data$sample,pattern = "^(\\d+).*$",replacement = "\\1")
all_data$error <- abs(all_data$y_pred - all_data$y)
all_data$method <- NULL

################################################################################
## RT based evaluation at experiment level. (no combination)
## peptide level
ind_data <- all_data %>% filter(level == "individual",peptide_or_psm=="peptide")
ind_data$PepQuery <- ifelse(ind_data$pepquery==1,"Pass","Fail")
#ind_data$error <- ifelse(ind_data$error>10,10,ind_data$error)
mae <- ind_data %>% group_by(FDR_method,search_engine,PepQuery) %>% summarise(mae=median(error))


gg <- ind_data %>% plot_boxplot(mae)
#stat_compare_means(comparisons = list(c(1,2),c(1,3)))

pdf("rt_individual_level_no_combine_peptide_level.pdf",width = 3,height = 4.5)
print(gg)
dev.off()

## RT based evaluation at experiment level. (no combination)
## level level
ind_data <- all_data %>% filter(level == "individual",peptide_or_psm=="PSM")
ind_data$PepQuery <- ifelse(ind_data$pepquery==1,"Pass","Fail")
#ind_data$error <- ifelse(ind_data$error>10,10,ind_data$error)
mae <- ind_data %>% group_by(FDR_method,search_engine,PepQuery) %>% summarise(mae=median(error))


gg <- ind_data %>% plot_boxplot(mae)
#stat_compare_means(comparisons = list(c(1,2),c(1,3)))
pdf("rt_individual_level_no_combine_PSM_level.pdf",width = 3,height = 4.5)
print(gg)
dev.off()


################################################################################
## peptide level combined
ind_data <- all_data %>% filter(level == "combine",peptide_or_psm=="peptide")
ind_data <- remove_redundant(ind_data)
mae <- ind_data %>% group_by(FDR_method,search_engine,PepQuery) %>% summarise(mae=median(error))

gg <- ind_data %>% plot_boxplot(mae)
    
pdf("rt_combine_peptide_level.pdf",width = 3,height = 4.5)
print(gg)
dev.off()


ind_data %>% group_by(search_engine,FDR_method) %>% summarise(pvalue=wilcox_test(error,PepQuery)) %>% write_tsv("combined_peptide_pvalue_before_after_pepquery.txt")
ind_data %>% filter(PepQuery == "Pass") %>% group_by(search_engine) %>% do(wilcox_test_multiple(.$error,.$FDR_method)) %>% write_tsv("combined_peptide_pvalue_after_pepquery.txt")

################################################################################
## PSM level combined
ind_data <- all_data %>% filter(level == "combine",peptide_or_psm=="PSM")
ind_data <- remove_redundant(ind_data)
mae <- ind_data %>% group_by(FDR_method,search_engine,PepQuery) %>% summarise(mae=median(error))

gg <- ind_data %>% plot_boxplot(mae)

pdf("rt_combine_PSM_level.pdf",width = 3,height = 4.5)
print(gg)
dev.off()


ind_data %>% group_by(search_engine,FDR_method) %>% summarise(pvalue=wilcox_test(error,PepQuery)) %>% write_tsv("combined_PSM_pvalue_before_after_pepquery.txt")
ind_data %>% filter(PepQuery == "Pass") %>% group_by(search_engine) %>% do(wilcox_test_multiple(.$error,.$FDR_method)) %>% write_tsv("combined_PSM_pvalue_after_pepquery.txt")

################################################################################
## peptide level individual
ind_data <- all_data %>% filter(level == "individual",peptide_or_psm=="peptide") %>% ungroup()
ind_data <- remove_redundant(ind_data)
mae <- ind_data %>% group_by(FDR_method,search_engine,PepQuery) %>% summarise(mae=median(error))

gg <- ind_data %>% plot_boxplot(mae)

pdf("rt_individual_peptide_level.pdf",width = 3,height = 4.5)
print(gg)
dev.off()


ind_data %>% group_by(search_engine,FDR_method) %>% summarise(pvalue=wilcox_test(error,PepQuery)) %>% write_tsv("individual_peptide_pvalue_before_after_pepquery.txt")
ind_data %>% filter(PepQuery == "Pass") %>% group_by(search_engine) %>% do(wilcox_test_multiple(.$error,.$FDR_method)) %>% write_tsv("individual_peptide_pvalue_after_pepquery.txt")

################################################################################
## PSM level individual
ind_data <- all_data %>% filter(level == "individual",peptide_or_psm=="PSM") %>% ungroup()
ind_data <- remove_redundant(ind_data)
mae <- ind_data %>% group_by(FDR_method,search_engine,PepQuery) %>% summarise(mae=median(error))

gg <- ind_data %>% plot_boxplot(mae)

pdf("rt_individual_PSM_level.pdf",width = 3,height = 4.5)
print(gg)
dev.off()


ind_data %>% group_by(search_engine,FDR_method) %>% summarise(pvalue=wilcox_test(error,PepQuery)) %>% write_tsv("individual_PSM_pvalue_before_after_pepquery.txt")
ind_data %>% filter(PepQuery == "Pass") %>% group_by(search_engine) %>% do(wilcox_test_multiple(.$error,.$FDR_method)) %>% write_tsv("individual_PSM_pvalue_after_pepquery.txt")


