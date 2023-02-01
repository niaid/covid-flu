library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggridges)
library(dplyr)

parameter.scatter.plot <- function(df,param1,param2,rug=F) {
  df$param1 <- df[,param1]
  df$param2 <- df[,param2]
  g <- ggpubr::ggscatter(df,x="param1",y="param2",rug=rug,alpha=0.6,color="sex",size=1.5,shape="group") + 
    stat_cor(method="spearman",size=3,cor.coef.name = "rho") + theme_bw() + scale_color_nejm()  + 
    theme(aspect.ratio = 0.9,panel.spacing = unit(1,"lines")) + xlab(param1) + ylab(param2)
  return(g)
}

# scatter plot with x-axis showing TSD (covid.diagnosis.start.date.to.sample.drawn) and y-axis representing the parameter of interest
tsd.scatter.plot <- function(df,param.colname,additional.facet.term = NA) {
  colnames(df) <- make.names(colnames(df))
  param.colname <- make.names(param.colname)
  # average of HCs
  hc.df <- subset(df,group == "HC")
  if (nrow(hc.df) > 0) {
    if (is.na(additional.facet.term)) {
      hc.median <- aggregate(hc.df[,param.colname],list(hc.df$sex),median,na.rm=T)
      colnames(hc.median) <- c("sex","value")
    } else {
      hc.median <- aggregate(hc.df[,param.colname],list(hc.df$sex,hc.df[,additional.facet.term]),median,na.rm=T)
      colnames(hc.median) <- c("sex",additional.facet.term,"value")
    }
  } 
  g <- ggplot(df,aes_string("covid.diagnosis.start.date.to.sample.drawn",param.colname,group="sex"))
  
  # add median line of healthy
  if (nrow(hc.df) > 0)
    g <- g + geom_hline(data=hc.median,aes(yintercept=value,color=sex),linetype="dashed",alpha=0.4,size=0.8) 
  
  g <-  g + 
    geom_point(aes(fill=sex),pch=21,alpha=0.5,size=2) + 
    stat_smooth(method="lm",formula = y ~ x,aes(color=sex,fill=sex),alpha=0.3) + 
    stat_cor(aes(color=sex),method = "spearman",cor.coef.name = "rho",size=3) +
    scale_color_nejm() + scale_fill_nejm() + 
    theme_bw() + xlab("Time since diagnosis (Days)") 
  return(g)
}

geneset.ridge.plot <- function(geneset.fc,highlight.genes=NULL) {
  if (nrow(geneset.fc) > 0) {
    # calculate summary statistics
    mean.abs.effect <- as.data.frame(geneset.fc %>% group_by(visit,sex) %>% summarise(n=length(z.std),mean=mean(z.std),
                                                                                      median=median(z.std)))
    mean.abs.effect$visit <- factor(mean.abs.effect$visit,c("D28","D1","Baseline"))
    #print(mean.abs.effect)
    mean.abs.effect$group.label <- paste0(mean.abs.effect$sex," (",mean.abs.effect$n," genes)")
    
    # ridge plots
    geneset.fc <- merge(geneset.fc,mean.abs.effect,by=c("sex","visit"))
    g <- ggplot(geneset.fc,aes(y=visit,x=z.std)) +
      geom_vline(xintercept = 0) + 
      geom_density_ridges(aes(fill=visit),quantile_lines=F,quantiles=2,scale=1,alpha=0.75,color="grey15",
                          jittered_points=T,point_size=1,point_alpha = 0.6,point_shape="|",
                          position = position_points_jitter(width=0,height=0)) +
      geom_segment(data = mean.abs.effect,aes(x=median,xend=median,y=as.numeric(visit),yend=as.numeric(visit) + 0.9),
                   color="red",linetype="dashed") + xlab("Effect Size") +
      scale_fill_brewer(palette = 4) + scale_color_d3() + theme_ridges() + scale_y_discrete(expand = c(0.01, 0)) +
      facet_grid(~group.label,scales = "free") + 
      theme(title = element_text(size=8),legend.position = "bottom",legend.box = "vertical",
            strip.background = element_rect(fill=NA),
            panel.spacing = unit(2,"lines")) 
    if (!is.null(highlight.genes))
      g <- g + geom_point(data=subset(geneset.fc,gene.symbol %in% highlight.genes),aes(color=gene.symbol),size=2,alpha=.95)
    return(g)
  }
}

overall.group.visit.change.boxplot <- function(df,param.colname,pval.df=NULL) {
  visit.days <- c("D0","D1","D7","D28")
  df$overall.group.visit <- factor(paste0(df$group,"\n",gsub("Day ","D",df$visit)),levels = c(paste0("HC\n",visit.days),paste0("COVR\n",visit.days)))
  
  g <- ggplot(df,aes_string("overall.group.visit",param.colname)) +
    geom_boxplot(aes(fill=sex),color = "black", alpha = 0.4, width = 0.75, outlier.shape = NA) +
    geom_line(aes(group=alt.subject.id,color=sex),alpha=0.75) +
    geom_point(aes(color=sex),size=0.5) +
    scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6")) +
    scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6")) +
    scale_y_continuous(expand = c(0,.2)) +
    facet_grid(~sex) + theme_bw(base_size = 14)
  if (is.null(pval.df)) {
    g <- g + stat_compare_means(comparisons = list(c(1,2),c(3,4)), method = "wilcox.test", size = 3,paired = T)
  } else {
    value.max <- max(df[,param.colname])
    value.min <- min(df[,param.colname])
    g <- g + stat_pvalue_manual(pval.df,size=3,y.position = (value.max + (value.max-value.min)*0.1))
  }
  return(g)
}

overall.sex.group.boxplot <- function(df,param.colname) {
  df$value <- df[,param.colname]
  df$overall.group.sex <- factor(paste0(df$group,"\n",df$sex),
                                 levels = c("HC\nFemale","COVR\nFemale","HC\nMale","COVR\nMale"))
  df$visit <- gsub("Day ","D",df$visit)
  
  g <- ggplot(df,aes(overall.group.sex,value)) + 
    geom_boxplot(aes(fill=sex,color=sex), alpha = 0.4, width = 0.3, outlier.shape = NA)+
    geom_point(aes(fill=sex), shape = 21, color = "white", size = 2, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
    scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6")) +
    scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6")) +
    stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,4),c(3,4)), method = "wilcox.test", size = 2.5)
  return(g)
}

overall.group.sex.boxplot <- function(df,param.colname) {
  df$value <- df[,param.colname]
  df$overall.group.sex <- factor(paste0(df$group,"\n",df$sex),
                                 levels = c("HC\nFemale","COVR\nFemale","HC\nMale","COVR\nMale"))
  g <- ggplot(df,aes(overall.group.sex,value)) + 
    geom_boxplot(aes(color=group),outlier.shape = NA,width=0.6) + 
    geom_jitter(aes(fill=sex),height = 0,width=0.1,alpha=0.75,pch=21,color="white",size=2) +
    scale_fill_nejm() + scale_color_jama() + xlab("") +
    theme_pubclean() + stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,4),c(3,4)),hide.ns = F,label = "p.format",tip.length = 0) +
    theme(legend.position = "None")
  return(g)
}

overall.group.sex.barplot <- function(df,param.colname,to.plot.colname=NA) {
  df$value <- df[,param.colname]
  df <- subset(df,!is.na(value))
  df$overall.group.sex <- factor(paste0(df$group,"\n",df$sex),
                                 levels = c("HC\nFemale","COVR\nFemale","HC\nMale","COVR\nMale"))
  df$visit <- gsub("Day ","D",df$visit)
  
  # calculate mean and std error for each bar
  df.stat <- df %>% group_by(overall.group.sex,group,sex) %>% summarise(mean=mean(value),sem=plotrix::std.error(value))
  if (!is.na(to.plot.colname)) {
    df$value <- df[,to.plot.colname]
    df <- subset(df,!is.na(value))
  }
  
  g <- ggplot(df,aes(overall.group.sex,value)) + 
    geom_hline(yintercept = 0) + 
    geom_jitter(aes(fill=sex),alpha=0.6,height=0,width=0.15,pch=21) +
    geom_bar(data=df.stat,aes(y=mean,color=group),stat="identity",width=0.4,fill=NA) + 
    geom_errorbar(data=df.stat,aes(y=mean,ymin=mean-sem,ymax=mean+sem,color=group),width=.2,position=position_dodge(.9),size=0.8) + 
    stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,4),c(3,4)),size=3) + theme_pubr() +
    ylab(param.colname) + xlab("") + scale_y_continuous(expand = expansion(mult=c(0.05,.12))) +
    scale_color_manual(values=c(setNames(pal_nejm()(2),c("Female","Male")),setNames(pal_jama()(2),c("COVR","HC")))) + 
    scale_fill_manual(values=c(setNames(pal_nejm()(2),c("Female","Male")),setNames(pal_jama()(2),c("COVR","HC")))) +
    theme(axis.line.x = element_line(color=NA),legend.position = "none") 
  return(g)
}









