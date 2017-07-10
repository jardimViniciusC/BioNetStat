library(pathview)

####################################################################
selecting.list <- function(var.diff.list,threshold,thr.value){
  if(dim(var.diff.list)[1]==0) stop("The differential table have 0 rows")
  if(is.null(threshold)) tab<-var.diff.list[,c(1,5:(dim(var.diff.list)[2]))]
    else{
      if(threshold=="pvalue")tab<-var.diff.list[var.diff.list$Nominal.p.value<=thr.value,c(1,5:(dim(var.diff.list)[2]))]
      if(threshold=="qvalue")tab<-var.diff.list[var.diff.list$Q.value<=thr.value,c(1,5:(dim(var.diff.list)[2]))]
      if(dim(tab)[1]==0) stop("Threshold is too low, the differential table remains with 0 rows")
    }
  rownames(tab)<-tab[,1]
  return(tab[,-1])
}

####################################################################
var.list<-function(expr,labels){
  mat<-cbind(labels,expr)
  media<-aggregate(mat,by = list(labels),FUN = median)[-1]
  tab<-t(scale(media[,-1]))#gene.id
  colnames(tab)<-media[,1]
  return(tab)
}

####################################################################
plot.names<-function(var.diff.list,threshold,thr.value){
  if(dim(var.diff.list)[1]==0) stop("The differential table have 0 rows")
    if(is.null(threshold)) names<-var.diff.list[,1]
    else{
      if(threshold=="pvalue")names<-var.diff.list[var.diff.list$Nominal.p.value<=thr.value,1]
      if(threshold=="qvalue")names<-var.diff.list[var.diff.list$Q.value<=thr.value,1]
      if(length(names)==0) stop("Threshold is too restrict, the differential table remains with 0 rows")
    }
  return(names)
}

####################################################################
centralityPathPlot<- function(gene.data=NULL, cpd.data=NULL, threshold=NULL, thr.value=0.05, species , out.suffix, pathway.id, kegg.native=T, file.name="path"){
  
  if(!is.null(gene.data)) {
    tab.gene<-selecting.list(gene.data,threshold=threshold,thr.value=thr.value)
    max.gene<-max(tab.gene)
  }
    else{
      tab.gene<-NULL
      max.gene<-1
    }
  if(!is.null(cpd.data)){
    tab.cpd<-selecting.list(cpd.data,threshold=threshold,thr.value=thr.value)
    max.cpd<-max(tab.cpd)
  }
    else{
      tab.cpd<-NULL
      max.cpd<-1
    }
  if(is.null(cpd.data) & is.null(gene.data)) stop("You have to insert a matrix data")
  pv.out <- pathview(gene.data=tab.gene,cpd.data = tab.cpd, pathway.id = pathway.id,
                     species = species, out.suffix = file.name, kegg.native = kegg.native,
                     limit = list(gene = max.gene, cdp = max.cpd), bins = list(gene = 15,cpd = 15), 
                     both.dirs= list(gene = F,cpd = F), mid =list(gene = "white", cpd = "white"),high = list(gene = "red",cpd = "red"))
                     
                     # limit = limit, bins = bins, both.dirs= both.dirs,
                     # mid =mid,high = high)
  return(pv.out)
}

####################################################################
pathPlot<- function(gene.data=NULL, cpd.data=NULL, labels, varr.diff.list=NULL, threshold=NULL, thr.value=0.05, species , pathway.id, kegg.native=T, file.name="path"){
    
    if(!is.null(gene.data)) {
      tab.gene<-var.list(gene.data,labels=labels)
      names.gene<-plot.names(gene.data,threshold=threshold,thr.value=thr.value)
      # max.gene<-max(tab.gene)
    }
    else{
      tab.gene<-NULL
      names.gene<-NULL
      # max.gene<-1
    }
    if(!is.null(cpd.data)){
      tab.cpd<-var.list(cpd.data,labels=labels)
      names.cpd<-plot.names(var.diff.list=a,threshold=threshold,thr.value=thr.value)
      # max.cpd<-max(tab.cpd)
    }
    else{
      tab.cpd<-NULL
      names.cpd<-NULL
      # max.cpd<-1
    }
  if(is.null(cpd.data) & is.null(gene.data)) stop("You have to insert a matrix data")
  pv.out <- pathview(gene.data = tab.gene[names.gene,], cpd.data = tab.cpd[names.cpd,], pathway.id = pathway.id,
                     species = species, out.suffix = file.name, kegg.native = kegg.native)#"05200"
}