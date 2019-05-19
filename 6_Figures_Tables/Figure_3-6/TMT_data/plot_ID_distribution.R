library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(VennDiagram)

load_data=function(f){
    a <- read_tsv(f)
    filename <- basename(f)
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
    se = "-"
    if(str_detect(x,pattern = "_comet_")){
        se = "Comet"
    }else if(str_detect(x,pattern = "_msgf_")){
        se = "MS-GF+"
    }else if(str_detect(x,pattern = "_tandem_")){
        se = "X!Tandem"    
    }else{
        stop(paste("Invalid search engine:",x,"\n",sep=""))
    }
    return(se)
}

detect_fdr_method=function(x){
    me = "-"
    if(str_detect(x,pattern = "_global_fdr")){
        me = "Global FDR"
    }else if(str_detect(x,pattern = "_separate_fdr")){
        me = "Separate FDR"
    }else if(str_detect(x,pattern = "_two_step")){
        me = "Two-stage FDR"    
    }else{
        stop(paste("Invalid FDR method:",x,"\n",sep=""))
    }
    return(me)
}

remove_redundant=function(ind_data){
    pep_pepquery <- ind_data %>% group_by(peptide,level,peptide_or_psm,search_engine,FDR_method) %>% 
        summarise(PepQuery=ifelse(any(pepquery==1),1,0)) %>%
        ungroup()
    ind_data$pepquery <- NULL
    ind_data <- inner_join(ind_data,pep_pepquery) 
    if("error" %in% names(ind_data)){
        ind_data <- ind_data %>% arrange(error) %>% group_by(peptide,level,peptide_or_psm,search_engine,FDR_method) %>% filter(row_number()==1)
        ind_data$error <- ifelse(ind_data$error>10,10,ind_data$error)
    }else{
        ind_data <- ind_data %>% select(peptide,level,peptide_or_psm,search_engine,FDR_method,PepQuery) %>% distinct()
    }
    ind_data$PepQuery <- ifelse(ind_data$PepQuery==1,"Pass","Fail")
    ind_data <- ind_data %>% ungroup()
    return(ind_data)
}

plot_upset=function(dat){
    x <- dat %>% filter(level=="individual",peptide_or_psm=="peptide")
    x$peptide <- str_replace_all(x$peptide,pattern = "I",replacement = "L")
    x$type <- paste(x$search_engine,x$FDR_method,sep=":")
    x <- x  %>% filter(pepquery==1) %>% select(type,peptide) %>% distinct()
    y <- list()
    a = sort(unique(x$type))
    for(tt in a){
        y[[tt]] <- x$peptide[x$type==tt]
    }
    
    pdf("individual_peptide_upset.pdf",width = 14,height = 6,onefile=FALSE)
    upset(fromList(y),nsets = length(a),nintersects = NA,order.by = "freq")
    dev.off()
    
    x <- dat %>% filter(level=="individual",peptide_or_psm=="PSM")
    x$peptide <- str_replace_all(x$peptide,pattern = "I",replacement = "L")
    x$type <- paste(x$search_engine,x$FDR_method,sep=":")
    x <- x  %>% filter(pepquery==1) %>% select(type,peptide) %>% distinct()
    y <- list()
    a = sort(unique(x$type))
    for(tt in a){
        y[[tt]] <- x$peptide[x$type==tt]
    }
    
    pdf("individual_PSM_upset.pdf",width = 14,height = 6,onefile=FALSE)
    upset(fromList(y),nsets = length(a),nintersects = NA,order.by = "freq")
    dev.off()
    
    x <- dat %>% filter(level=="combine",peptide_or_psm=="peptide")
    x$peptide <- str_replace_all(x$peptide,pattern = "I",replacement = "L")
    x$type <- paste(x$search_engine,x$FDR_method,sep=":")
    x <- x  %>% filter(pepquery==1) %>% select(type,peptide) %>% distinct()
    y <- list()
    a = sort(unique(x$type))
    for(tt in a){
        y[[tt]] <- x$peptide[x$type==tt]
    }
    
    pdf("combine_peptide_upset.pdf",width = 14,height = 6,onefile=FALSE)
    upset(fromList(y),nsets = length(a),nintersects = NA,order.by = "freq")
    dev.off()
    
    x <- dat %>% filter(level=="combine",peptide_or_psm=="PSM")
    x$peptide <- str_replace_all(x$peptide,pattern = "I",replacement = "L")
    x$type <- paste(x$search_engine,x$FDR_method,sep=":")
    x <- x  %>% filter(pepquery==1) %>% select(type,peptide) %>% distinct()
    y <- list()
    a = sort(unique(x$type))
    for(tt in a){
        y[[tt]] <- x$peptide[x$type==tt]
    }
    
    pdf("combine_PSM_upset.pdf",width = 14,height = 6,onefile=FALSE)
    upset(fromList(y),nsets = length(a),nintersects = NA,order.by = "freq")
    dev.off()
    
}

## one search engine, one figure
plotVenn=function(dat,level_type = "individual",peptide_or_psm_type="peptide"){
    
    dir.create("venn_search_engine",showWarnings = FALSE)
    x <- dat %>% filter(level==level_type,peptide_or_psm==peptide_or_psm_type)
    se_list <- sort(unique(dat$search_engine))
    x$peptide <- str_replace_all(x$peptide,pattern = "I",replacement = "L")
    #x$type <- paste(x$search_engine,x$FDR_method,sep=":")
    x <- x  %>% filter(pepquery==1) %>% select(peptide,FDR_method,search_engine) %>% distinct()
    for(se in se_list){
        y <- list()
        dd <- x %>% filter(search_engine==se)
        a = sort(unique(dd$FDR_method))
        for(tt in a){
            y[[tt]] <- dd$peptide[dd$FDR_method==tt]
        }
        
        figure <- paste("venn_search_engine/",se,level_type,peptide_or_psm_type,"-venn.pdf",sep="")
        cat(level_type," + ",peptide_or_psm_type," => ",figure,"\n");
        pdf(file = figure,width = 3,height = 3)
        
        vdat <- y
        vplot <- venn.diagram(vdat,filename = NULL,lty = rep("blank", 3), 
                              fill = c("red", "pink1","skyblue"), 
                              cat.col= c("red", "pink1","skyblue"), 
                              cex=0.8,
                              alpha = rep(0.5, 3),margin=0.05)
        grid.draw(vplot)
        dev.off()
    }
}

## same method, different search engines
plotVenn_method=function(dat,level_type = "individual",peptide_or_psm_type="peptide"){
    
    dir.create("venn_fdr_method",showWarnings = FALSE)
    
    x <- dat %>% filter(level==level_type,peptide_or_psm==peptide_or_psm_type)
    me_list <- sort(unique(dat$FDR_method))
    x$peptide <- str_replace_all(x$peptide,pattern = "I",replacement = "L")
    #x$type <- paste(x$search_engine,x$FDR_method,sep=":")
    x <- x  %>% filter(pepquery==1) %>% select(peptide,FDR_method,search_engine) %>% distinct()
    for(se in me_list){
        y <- list()
        dd <- x %>% filter(FDR_method==se)
        a = sort(unique(dd$search_engine))
        for(tt in a){
            y[[tt]] <- dd$peptide[dd$search_engine==tt]
        }
        
        figure <- paste("venn_fdr_method/",se,level_type,peptide_or_psm_type,"-venn.pdf",sep="")
        cat(level_type," + ",peptide_or_psm_type," => ",figure,"\n");
        pdf(file = figure,width = 3,height = 3)
        
        vdat <- y
        vplot <- venn.diagram(vdat,filename = NULL,lty = rep("blank", 3), 
                              fill = c("red", "pink1","skyblue"), 
                              cat.col= c("red", "pink1","skyblue"), 
                              cex=0.8,
                              alpha = rep(0.5, 3),margin=0.05)
        grid.draw(vplot)
        dev.off()
    }
}


## input parameter
input_dir <- "results_summary/prediction_v2/"
sample_or_experiment <- "TMT Experiment" # or "TMT Experiment" => TMT or iTRAQ
sample_or_experiment_rotate <- FALSE # or FALSE
process_sample_name = TRUE
prospective_colon_sample=FALSE

## load data
input_files <- list.files(path = input_dir,pattern = ".txt",full.names = TRUE,include.dirs = TRUE)
all_data <- bind_rows(lapply(input_files, load_data))
if(process_sample_name){
    all_data$sample <- str_replace_all(all_data$sample,pattern = "^(\\d+).*$",replacement = "\\1")
}

if(prospective_colon_sample==TRUE){
    sample_to_n <- data.frame(sample=all_data$sample %>% unique() %>% sort,stringsAsFactors = FALSE)
    sample_to_n$sample_n = 1:nrow(sample_to_n)
    all_data <- merge(all_data,sample_to_n,by="sample")
    all_data$sample <- all_data$sample_n
}


plotVenn(all_data,level_type = "individual",peptide_or_psm_type="peptide")
plotVenn(all_data,level_type = "individual",peptide_or_psm_type="PSM")
plotVenn(all_data,level_type = "combine",peptide_or_psm_type="PSM")
plotVenn(all_data,level_type = "combine",peptide_or_psm_type="peptide")

plotVenn_method(all_data,level_type = "individual",peptide_or_psm_type="peptide")
plotVenn_method(all_data,level_type = "individual",peptide_or_psm_type="PSM")
plotVenn_method(all_data,level_type = "combine",peptide_or_psm_type="PSM")
plotVenn_method(all_data,level_type = "combine",peptide_or_psm_type="peptide")

################################################################################
## Variant peptide identification distribution: individual level
ind_data <- all_data %>% filter(level == "individual")
dd <- ind_data %>% group_by(sample,peptide_or_psm,search_engine,FDR_method) %>% 
    summarise(n=length(unique(peptide)))
dd_pepquery <- ind_data %>% filter(pepquery == 1) %>% group_by(sample,peptide_or_psm,search_engine,FDR_method) %>% 
    summarise(n=length(unique(peptide)))
dd$PepQuery <- "No"
dd_pepquery$PepQuery <- "Yes"
max_y <- max(dd$n)*1.1
gg <- bind_rows(dd,dd_pepquery) %>% mutate(linetype=paste(PepQuery,FDR_method,sep="_")) %>% filter(peptide_or_psm=="peptide") %>% 
    ggplot(aes(x=sample,y=n,colour=FDR_method))+
    #geom_line(aes(group=FDR_method,linetype=PepQuery))+
    geom_line(aes(group=linetype,linetype=PepQuery))+
    geom_point(size=0.8)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab(sample_or_experiment)+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    ylim(0,max_y)+
    scale_color_manual(values=c("red", "blue", "green"))
if(sample_or_experiment_rotate){
    gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
pdf("Experiment_level_peptide_level_distribution.pdf",width = 6,height = 4.5)
print(gg)
dev.off()

gg <- bind_rows(dd,dd_pepquery) %>% mutate(linetype=paste(PepQuery,FDR_method,sep="_")) %>% filter(peptide_or_psm=="PSM") %>% 
    ggplot(aes(x=sample,y=n,colour=FDR_method))+
    #geom_line(aes(group=FDR_method,linetype=PepQuery))+
    geom_line(aes(group=linetype,linetype=PepQuery))+
    geom_point(size=0.8)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab(sample_or_experiment)+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    ylim(0,max_y)+
    scale_color_manual(values=c("red", "blue", "green"))
if(sample_or_experiment_rotate){
    gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
pdf("Experiment_level_PSM_level_distribution.pdf",width = 6,height = 4.5)
print(gg)
dev.off()


################################################################################
#### add side barplot
ind_data <- all_data %>% filter(level == "individual")
dd <- ind_data %>% group_by(sample,peptide_or_psm,search_engine,FDR_method,pepquery) %>% 
    summarise(n=length(unique(peptide)))
dat <- dd %>% spread(key = pepquery, value = n,fill = 0)
#dat$ratio <- dat$`0`/dat$`1`

## barplot at peptide level
x <- dat %>% filter(peptide_or_psm=="peptide") %>% group_by(search_engine,FDR_method) %>% 
    summarise(Fail=median(`0`),Pass=median(`1`)) %>%
    mutate(ratio=Fail/(Pass+Fail))
x1 <- x %>% gather(key = "PepQuery",value = "n",-search_engine,-FDR_method,-ratio) 
text_y <- 1.2*max(x1$n)
gg <- x1 %>% ggplot(aes(x=FDR_method,y=n,fill=PepQuery))+
    geom_bar(stat="identity",width = 0.5)+
    #geom_text(aes(x=FDR_method,y=))
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=2.5)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab("FDR estimation method")+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    geom_text(aes(x=FDR_method,y=text_y,label=sprintf("%.2f%%",100*ratio)),size=2.5)+
    theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
    ylim(0,text_y*1.1)

pdf("Experiment_level_peptide_level_bar.pdf",width = 2,height = 4.5)
print(gg)
dev.off()


## barplot at PSM level
x <- dat %>% filter(peptide_or_psm=="PSM") %>% group_by(search_engine,FDR_method) %>% summarise(Fail=median(`0`),Pass=median(`1`)) %>%
    mutate(ratio=Fail/(Pass+Fail))
x1 <- x %>% gather(key = "PepQuery",value = "n",-search_engine,-FDR_method,-ratio) 

text_y <- 1.2*max(x1$n)
gg <- x1 %>% ggplot(aes(x=FDR_method,y=n,fill=PepQuery))+
    geom_bar(stat="identity",width = 0.5)+
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=2.5)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab("FDR estimation method")+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    geom_text(aes(x=FDR_method,y=text_y,label=sprintf("%.2f%%",100*ratio)),size=2.5)+
    theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
    ylim(0,text_y*1.1)

pdf("Experiment_level_PSM_level_bar.pdf",width = 2,height = 4.5)
print(gg)
dev.off()


################################################################################
## combined individual
dat <- all_data %>% remove_redundant() %>%
    group_by(level,peptide_or_psm,search_engine,FDR_method,PepQuery) %>%
    summarise(n=n()) %>% 
    spread(key=PepQuery,value = n)

text_y <- 1.05*max(dat$Pass+dat$Fail)

dat <- dat %>% mutate(ratio=Fail/(Pass+Fail)) %>%
    gather(key="PepQuery",value="n",-level,-peptide_or_psm,-search_engine,-FDR_method,-ratio)



gg <- dat %>% filter(peptide_or_psm=="PSM",level=="combine") %>% ggplot(aes(x=FDR_method,y=n,fill=PepQuery))+
    geom_bar(stat="identity",width = 0.5)+
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=2.5)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab("FDR estimation method")+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    geom_text(aes(x=FDR_method,y=text_y,label=sprintf("%.2f%%",100*ratio)),size=2.5)+
    theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
    ylim(0,text_y*1.15)
pdf("combine_psm_bar.pdf",width = 2,height = 4.5)
print(gg)
dev.off()

gg <- dat %>% filter(peptide_or_psm=="peptide",level=="combine") %>% ggplot(aes(x=FDR_method,y=n,fill=PepQuery))+
    geom_bar(stat="identity",width = 0.5)+
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=2.5)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab("FDR estimation method")+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    geom_text(aes(x=FDR_method,y=text_y,label=sprintf("%.2f%%",100*ratio)),size=2.5)+
    theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
    ylim(0,text_y*1.1)

pdf("combine_peptide_bar.pdf",width = 2,height = 4.5)
print(gg)
dev.off()

gg <- dat %>% filter(peptide_or_psm=="PSM",level=="individual") %>% ggplot(aes(x=FDR_method,y=n,fill=PepQuery))+
    geom_bar(stat="identity",width = 0.5)+
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=2.5)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab("FDR estimation method")+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    geom_text(aes(x=FDR_method,y=text_y,label=sprintf("%.2f%%",100*ratio)),size=2.5)+
    theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
    ylim(0,text_y*1.1)

pdf("individual_psm_bar.pdf",width = 2,height = 4.5)
print(gg)
dev.off()

gg <- dat %>% filter(peptide_or_psm=="peptide",level=="individual") %>% ggplot(aes(x=FDR_method,y=n,fill=PepQuery))+
    geom_bar(stat="identity",width = 0.5)+
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),size=2.5)+
    facet_grid(search_engine~.)+
    ylab("# novel peptides")+
    xlab("FDR estimation method")+
    theme_bw()+
    theme(legend.position="top",legend.text=element_text(size=7),
          legend.title=element_text(size=7.5),legend.box.spacing = unit(0, "inch"))+
    geom_text(aes(x=FDR_method,y=text_y,label=sprintf("%.2f%%",100*ratio)),size=2.5)+
    theme(axis.text.x = element_text(size=9,angle = 45,hjust = 1))+
    ylim(0,text_y*1.1)
pdf("individual_peptide_bar.pdf",width = 2,height = 4.5)
print(gg)
dev.off()






