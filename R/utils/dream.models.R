loadLibraries <- function() {
  library(variancePartition)
  library(emmeans)
  library(ggpubr)
  library(ggsci)
  library(dplyr)
  library(lmerTest)
}

# add min non-zero value of column divided by d (d=2 by default) to each column with one or more 0 values
# to allow for log transformation
add.dynamic.offset <- function(data,d=2) {
  new.data <- sapply(1:ncol(data),function(x){
    if(sum(data[,x] == 0,na.rm = T) > 0){
      col.values <- data[,x]
      data[,x] + min(col.values[col.values > 0],na.rm = T)/d
    }else{
      data[,x]
    }})
  dimnames(new.data) <- dimnames(data)
  return(new.data)
}

# run geneset enrichment analysis separately for each term in the input
# default to use signed -log10(P.Value) to order the gene list
# return both intra-ontology and cross-ontology FDR (adj.P.Val and p.adjust, respectively) using the BH method
run.gsea <- function(model.de.res,gmt.df,sort.by.column=c("P.Value","logFC","z.std"),minSize=10,maxSize=300,plot.genesets = c()) {
  set.seed(786)
  overall.gsea.res <- data.frame()
  selected <- match.arg(sort.by.column)
  cat("Order genes by",selected,"\n")
  for (i in unique(model.de.res$term)) {
    tmp <- subset(model.de.res,term == i)
    gene.list <- tmp[,selected]
    if(selected == "P.Value") 
      gene.list <- -log10(tmp[,selected])*sign(tmp[,"logFC"])
    names(gene.list) <- tmp$gene.symbol
    
    gsea.obj <- GSEA(sort(gene.list,decreasing = T),pvalueCutoff = 1,nPerm = 10000,TERM2GENE = gmt.df,
                                   minGSSize = minSize,maxGSSize = maxSize,seed = T)
    for (to.plot in plot.genesets) {
      print(enrichplot::gseaplot2(gsea.obj,to.plot,title = paste0(i," - ",to.plot),pvalue_table = T))
    }
    gsea.res <- as.data.frame(gsea.obj)
    gsea.res$ont <- gmt.df[match(gsea.res$ID,gmt.df$ont),"source"]
    gsea.res <- gsea.res %>% group_by(ont) %>% mutate(adj.P.Val= p.adjust(pvalue,method="BH"))
    overall.gsea.res <- rbind(overall.gsea.res,data.frame(term=i,gsea.res))
  } 
  colnames(overall.gsea.res)[c(2,7)] <- c("parameter","P.Value")
  return(overall.gsea.res[,-3])
}

# compare COVR vs. HC and COVR males vs. females at baseline (pre-vaccination) timepoints, i.e. day -7 and/or day 0
# include subjects of all ages for the baseline
# default model: ~ 0 + overall.group.sex + age + race + (1|alt.subject.id)
# contrasts -
# COVR.v.HC: (COVR.m + COVR.f)/2 - (HC.m + HC.f)/2
# COVR.male.v.female: (COVR.m - COVR.f) - (HC.m - HC.f)
# Note: exprObj must be already normally distributed and can be a matrix of expression data 
# or the output returned by either voom or voomWithDreamWeights
baseline.condition.model <- function(exprObj,metadata,random.effects="alt.subject.id",additional.fixed.effects=NA,bayes=T) {
  weights <- NULL
  if (class(exprObj) == "EList") {
    weights = exprObj$weights
    exprObj <- exprObj$E
  }
  stopifnot(colnames(exprObj) == rownames(metadata)) 
  stopifnot(sum(is.na(exprObj)) == 0) 
  rownames(metadata) <- paste0("Sample",1:nrow(metadata))
  colnames(exprObj) <- rownames(metadata)
  if (!is.null(weights)) colnames(weights) <- rownames(metadata)
  include.samples <- subset(metadata,alt.subject.id != "Control")
  include.samples$overall.group.sex <- interaction(include.samples$group,include.samples$sex)
  
  # build model
  model <- "~ 0 + overall.group.sex + age + race"
  if (!is.na(additional.fixed.effects)) {
    model <- paste0(model," + ",paste0(additional.fixed.effects,collapse = " + "))
  } 
  random.ef <- NA
  if (!is.na(random.effects) & random.effects != "") {
    random.ef <- paste0("(1|",paste0(random.effects,collapse = ") + (1|"),")")
    model <- paste0(model," + ",random.ef)
  }
  cat("Model: ",model,"\n")
  model <- as.formula(model)
  
  # build contrast
  COVR.m <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                         c(paste0("overall.group.sexCOVR.Male"),paste0("overall.group.sexHC.Male")))
  COVR.f <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                         c(paste0("overall.group.sexCOVR.Female"),paste0("overall.group.sexHC.Female")))
  COVR.male.v.female <- COVR.m - COVR.f
  COVR.v.HC <- (COVR.m + COVR.f)/2
  L <- cbind(COVR.v.HC,COVR.male.v.female,COVR.m,COVR.f)
  
  # remove parameters that are constant
  param.var <- apply(exprObj[,rownames(include.samples)],1,var)
  if (sum(param.var == 0) > 0)
    cat("Removing parameters with 0 variance:",names(which(param.var == 0)),"\n")
  
  # run dream
  fit.baseline.all <- dream(exprObj[param.var > 0,rownames(include.samples)],model,include.samples,L = L,quiet=T,
                            weightsMatrix = weights[param.var > 0,rownames(include.samples)],suppressWarnings = T)
  # if there are no random factors, limma was used and need to do eBayes
  if (is.na(random.ef) & bayes == T)
    fit.baseline.all <- eBayes(fit.baseline.all,trend=T,robust=T)
  fit.baseline.all$model <- model
  return(fit.baseline.all)
}

# determine which parameters are associated with time since diagnosis in the COVR group and whether the trajectories are 
# different between males and females
# include subjects of all ages for the baseline
# default model: ~ 0 + sex + sex:TSD + age + race + (1|alt.subject.id)
# contrasts -
# TSD: (sexFemale:TSD + sexMale:TSD)/2
# TSD.male.v.female: sexMale:TSD - sexFemale:TSD
# Note: exprObj must be already normally distributed and can be a matrix of expression data 
# or the output returned by either voom or voomWithDreamWeights
baseline.tsd.association.model <- function(exprObj,metadata,random.effects="alt.subject.id",additional.fixed.effects=NA,bayes=T) {
  weights <- NULL
  if (class(exprObj) == "EList") {
    weights = exprObj$weights
    exprObj <- exprObj$E
  }
  stopifnot(colnames(exprObj) == rownames(metadata)) 
  rownames(metadata) <- paste0("Sample",1:nrow(metadata))
  colnames(exprObj) <- rownames(metadata)
  if (!is.null(weights)) colnames(weights) <- rownames(metadata)
  include.samples <- subset(metadata,!is.na(covid.diagnosis.start.date.to.sample.drawn))
  include.samples$TSD <- scale(include.samples$covid.diagnosis.start.date.to.sample.drawn)
  stopifnot(sum(is.na(exprObj[,rownames(include.samples)])) == 0) 
  
  # build model
  model <- "~ sex + sex:TSD + age + race"
  if (!is.na(additional.fixed.effects)) {
    model <- paste0(model," + ",paste0(additional.fixed.effects,collapse = " + "))
  } 
  random.ef <- NA
  if (!is.na(random.effects) & random.effects != "") {
    random.ef <- paste0("(1|",paste0(random.effects,collapse = ") + (1|"),")")
    model <- paste0(model," + ",random.ef)
  }
  cat("Model: ",model,"\n")
  model <- as.formula(model)
  
  # # build contrast
  TSD.male.v.female <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                         c(paste0("sexMale:TSD"),paste0("sexFemale:TSD")))
  TSD <- TSD.male.v.female
  TSD[TSD != 0] <- 0.5
  L <- cbind(TSD,TSD.male.v.female)
  
  # remove parameters that are constant
  param.var <- apply(exprObj[,rownames(include.samples)],1,var)
  if (sum(param.var == 0) > 0)
    cat("Removing parameters with 0 variance:",names(which(param.var == 0)),"\n")
  
  # run dream
  fit.baseline.tsd <- dream(exprObj[param.var > 0,rownames(include.samples)],
                            model,include.samples,L=L,weightsMatrix = weights[param.var > 0,rownames(include.samples)],suppressWarnings = T)
  # if there are no random factors, limma was used and need to do eBayes
  if (is.na(random.ef) & bayes == T)
    fit.baseline.tsd <- eBayes(fit.baseline.tsd,trend=T,robust=T)
  fit.baseline.tsd$model <- model
  
  return(fit.baseline.tsd)
}

# Timepoint model comparing COVR vs. HC and COVR male vs. female at each timepoint and changes of these contrasts relative to baseline
# Additional contrasts include changes within each group relative to their own baseline and COVR-M/F relative to the baseline of their HC counterpart
# Exclude subjects > 65 years of age
# default model: ~ 0 + visit.overall.group.sex + age + race + flu.vax.count.10yr + (1|alt.subject.id)
# contrasts -
# individual timepont (d1, d7, and d28, etc.) contrast to baseline (day -7 and day 0) for each group:sex 
# timepoint - baseline contrast for COVR vs. HC for male and female separately
# timepoint contrast for COVR male vs. female using both healthy and baseline as reference
timepoint.condition.model <- function(exprObj,metadata,baseline.timepoints,include.vax.history=T,
                                      additional.fixed.effects=NA,additional.random.effects=NA) {
  weights <- NULL
  if (class(exprObj) == "EList") {
    weights = exprObj$weights
    exprObj <- exprObj$E
  }
  stopifnot(colnames(exprObj) == rownames(metadata)) 
  stopifnot(sum(is.na(exprObj)) == 0) 
  rownames(metadata) <- paste0("Sample",1:nrow(metadata))
  colnames(exprObj) <- rownames(metadata)
  if (!is.null(weights)) colnames(weights) <- rownames(metadata)
  include.samples <- subset(metadata,age < 65 & alt.subject.id != "Control" & !is.na(flu.vax.count.10yr))
  include.samples$visit.group <- as.character(include.samples$visit)
  include.samples[which(include.samples$visit.group %in% baseline.timepoints),"visit.group"] <- "Baseline" 
  include.samples$visit.group <- gsub("Day ","D",include.samples$visit.group)
  include.samples$visit.overall.group.sex <- make.names(interaction(include.samples$visit.group,include.samples$group,include.samples$sex))
  
  # build model
  fixed.ef <- "~ 0 + visit.overall.group.sex + age + race"
  if (include.vax.history) fixed.ef <- paste0(fixed.ef," + flu.vax.count.10yr")
  if (!is.na(additional.fixed.effects))
    fixed.ef <- paste0(fixed.ef," + ",paste0(additional.fixed.effects,collapse = " + "))
  random.ef <- "(1|alt.subject.id)"
  if (!is.na(additional.random.effects))
    random.ef <- paste0(random.ef," + (1|",paste0(additional.random.effects,collapse = "+ (1|"),")")
  model <- as.formula(paste0(fixed.ef," + ",random.ef))
  cat("Model: ",paste0(fixed.ef," + ",random.ef),"\n")
  
  timepoint.contrast <- list()
  baseline.m.diff <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                      c(paste0("visit.overall.group.sex","Baseline",".COVR.Male"),
                        paste0("visit.overall.group.sex","Baseline",".HC.Male")))
  baseline.f.diff <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                      c(paste0("visit.overall.group.sex","Baseline",".COVR.Female"),
                        paste0("visit.overall.group.sex","Baseline",".HC.Female")))
  timepoint.contrast[["Baseline.COVR.m"]] <- baseline.m.diff
  timepoint.contrast[["Baseline.COVR.f"]] <- baseline.f.diff
  timepoint.contrast[["Baseline.sex.diff"]] <- baseline.m.diff - baseline.f.diff
  timepoint.contrast[["Baseline.COVR.v.HC"]] <- (baseline.m.diff + baseline.f.diff)/2
  to.contrast <- setdiff(make.names(unique(include.samples$visit.group)),"Baseline")
  for (i in to.contrast) {
    timepoint.m.diff <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                           c(paste0("visit.overall.group.sex",i,".COVR.Male"),paste0("visit.overall.group.sex",i,".HC.Male")))
    timepoint.v.baseline.HC.m <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                                    c(paste0("visit.overall.group.sex",i,".COVR.Male"),paste0("visit.overall.group.sex","Baseline",".HC.Male")))
    timepoint.f.diff <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                           c(paste0("visit.overall.group.sex",i,".COVR.Female"),paste0("visit.overall.group.sex",i,".HC.Female")))
    timepoint.v.baseline.HC.f <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                                    c(paste0("visit.overall.group.sex",i,".COVR.Female"),paste0("visit.overall.group.sex","Baseline",".HC.Female")))
    L3 <- (timepoint.m.diff - timepoint.f.diff) 
    L3.sub.baseline <- L3 - (baseline.m.diff - baseline.f.diff)
    L4 <- (timepoint.m.diff + timepoint.f.diff)/2 
    L4.sub.baseline <- L4 - (baseline.m.diff + baseline.f.diff)/2
    timepoint.contrast[[paste0(i,".COVR.m")]] <- timepoint.m.diff
    timepoint.contrast[[paste0(i,".COVR.v.baseline.HC.m")]] <- timepoint.v.baseline.HC.m
    timepoint.contrast[[paste0(i,".COVR.m.sub.baseline")]] <- timepoint.m.diff - timepoint.contrast[["Baseline.COVR.m"]]
    timepoint.contrast[[paste0(i,".COVR.f")]] <- timepoint.f.diff
    timepoint.contrast[[paste0(i,".COVR.v.baseline.HC.f")]] <- timepoint.v.baseline.HC.f
    timepoint.contrast[[paste0(i,".sex.diff.v.baseline.HC")]] <- timepoint.v.baseline.HC.m - timepoint.v.baseline.HC.f
    timepoint.contrast[[paste0(i,".COVR.f.sub.baseline")]] <- timepoint.f.diff - timepoint.contrast[["Baseline.COVR.f"]]
    timepoint.contrast[[paste0(i,".sex.diff")]] <- L3
    timepoint.contrast[[paste0(i,".sex.diff.sub.baseline")]] <- L3.sub.baseline
    timepoint.contrast[[paste0(i,".COVR.v.HC")]] <- L4
    timepoint.contrast[[paste0(i,".COVR.v.HC.sub.baseline")]] <- L4.sub.baseline
  }
  # timepoint vs. baseline contrast within individual groups
  for (i in to.contrast) {
    for (g in c("COVR.Male","COVR.Female","HC.Male","HC.Female")) {
      timepoint.v.baseline <- getContrast(exprObj[,rownames(include.samples)], model, include.samples,
                                    c(paste0("visit.overall.group.sex",i,".",g),paste0("visit.overall.group.sex","Baseline",".",g)))
      timepoint.contrast[[paste0(i,".v.baseline.",g)]] <- timepoint.v.baseline
    }
  }
  timepoint.contrast <- do.call("cbind",timepoint.contrast)
  
  # remove parameters that are constant
  param.var <- apply(exprObj[,rownames(include.samples),drop=F],1,var)
  if (sum(param.var == 0) > 0)
    cat("Removing parameters with 0 variance:",names(which(param.var == 0)),"\n")
  
  # run dream
  fit.timepoint <- dream(exprObj[param.var > 0,rownames(include.samples),drop=F],model,include.samples,L=timepoint.contrast,
                         weightsMatrix = weights[param.var > 0,rownames(include.samples),drop=F],suppressWarnings = T)
  fit.timepoint$model <- model
  return(fit.timepoint)
}

# assume that the values in exprObj is already in log space, calculte FC by subtracting 
# baseline mean values from all values (e.g. D7 - D0)
calculate.fold.change <- function(exprObj,metadata,baseline.timepoints,sample.suffix="") {
  stopifnot(colnames(exprObj) == rownames(metadata)) 
  orig.sample.id <- rownames(metadata)
  colnames(exprObj) <- paste0("S",colnames(exprObj))
  rownames(metadata) <- paste0("S",rownames(metadata))
  fc.values <- list()
  baseline.mean <- list()
  for (subj in unique(metadata$alt.subject.id)) {
    baseline.mean[[subj]] <- rowMeans(exprObj[,which(metadata$alt.subject.id == subj & metadata$visit %in% baseline.timepoints),drop=F],na.rm = T)
    fc.values[[subj]] <- exprObj[,which(metadata$alt.subject.id == subj)] - baseline.mean[[subj]] 
  }
  baseline.mean <- do.call(cbind,baseline.mean)
  fc.values <- do.call(cbind,fc.values)
  fc.values <- fc.values[,rownames(metadata)]
  colnames(fc.values) <- paste0(orig.sample.id,sample.suffix) 
  return(list(fc=fc.values,baseline=baseline.mean))
}

# extract model coefficients for specified terms
model.term.statistics <- function(model.fit,terms) {
  de.results <- data.frame()
  for (i in terms) {
    param.stat <- topTable(model.fit,coef = i,number = Inf)
    de.results <- rbind(de.results,data.frame(term=i,parameter=rownames(param.stat),param.stat))
  }
  rownames(de.results) <- paste0("res",1:nrow(de.results))
  return(de.results[order(de.results$adj.P.Val),])
}






