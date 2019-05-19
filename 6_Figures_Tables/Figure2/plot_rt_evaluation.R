
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)

# load AutoRT result
read_autort=function(x){
    dat <- read_tsv(x)
    file_name <- basename(dirname(x))
    dat$sample <- file_name
    return(dat)
}

read_gptime=function(x){
    dat <- read_delim(x,col_names = FALSE,delim = " ")
    dat <- dat[,1:3]
    names(dat) <- c("x","y","y_pred")
    file_name <- basename(dirname(x))
    dat$sample <- file_name
    return(dat)
}

read_train=function(x){
    dat <- read_tsv(x)
    file_name <- basename(x) %>% str_replace_all(pattern = "_train.tsv$",replacement = "")
    dat$sample <- file_name
    return(dat)
}

calc_metrics=function(x){
    x$ae <- abs(x$y_pred - x$y) 
    dat <- x %>% group_by(sample) %>%
        summarise(mae=median(ae))
    return(dat)
}

plot_mae=function(x,y,x_lab,y_lab){
    x$method <- x_lab
    y$method <- y_lab
    dat <- bind_rows(x,y)
    median_v <- dat %>% group_by(method) %>% summarise(mae=median(mae))
    print(median_v)
    gg <- ggplot(dat,aes(x=method,y=mae,colour=method))+
        geom_violin()+
        geom_boxplot(width=0.3)+
        geom_text(data=median_v,aes(x=method,y=mae,label=sprintf("%.2f",mae)),hjust=-0.35,size=3.5)+
        theme_bw()+
        theme(legend.position="none")+
        xlab("Method")+
        ylab("MAE (minute)")+
        ylim(0.3,4)
    return(gg)
}

plot_mae_size=function(x,y,z,x_lab,y_lab){
    x$Method <- x_lab
    y$Method <- y_lab
    x <- merge(x,z)
    y <- merge(y,z)
    dat <- bind_rows(x,y)
    gg <- ggplot(dat,aes(x=n,y=mae,colour=Method))+
        geom_point(alpha=0.6,size=0.8)+
        #theme(legend.position="none")+
        xlab("Size")+
        ylab("MAE (minute)")+
        geom_smooth()+
        theme_bw()+
        theme(legend.justification = c(1, 1), legend.position = c(1, 1),
              legend.background = element_rect(fill=alpha('white', 0)))+
        ylim(0.3,4)
    return(gg)
    
}


autort_test_dir <- "autort/test/"
gptime_test_dir <- "gptime/test/"
training_data_dir <- "training_testing_data/"

## AutoRT testing result
autort_test_files <- list.files(path = autort_test_dir,pattern = "test.csv",
                                recursive = TRUE,full.names = TRUE,
                                include.dirs = TRUE)

autort_data <- bind_rows(lapply(autort_test_files, read_autort))

## GPTime testing result
gptime_test_files <- list.files(path = gptime_test_dir,pattern = "test.tsv",
                                recursive = TRUE,full.names = TRUE,
                                include.dirs = TRUE)

gptime_data <- bind_rows(lapply(gptime_test_files, read_gptime))


## load training data
training_data_files <- list.files(path = training_data_dir,pattern = "_train.tsv",
                                recursive = TRUE,full.names = TRUE,
                                include.dirs = TRUE)

training_data <- bind_rows(lapply(training_data_files, read_train))
training_data <- training_data %>% group_by(sample) %>% summarise(n=n())


## metrics
au <- calc_metrics(autort_data)
gp <- calc_metrics(gptime_data)

pdf("mae.pdf",width = 2,height = 3)
print(plot_mae(au,gp,"AutoRT","GPTime"))
dev.off()

## 
pdf("mae_size.pdf",width = 3,height = 3)
print(plot_mae_size(au,gp,training_data,"AutoRT","GPTime"))
dev.off()

merge(au,training_data) -> m;lm(m$mae~m$n);cor.test(m$mae,m$n,method = "spe")
merge(gp,training_data) -> m;lm(m$mae~m$n);cor.test(m$mae,m$n,method = "spe")


