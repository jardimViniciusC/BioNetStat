library(pathview)

####################################################################
selecting.list <- function(var.diff.list,threshold,thr.value){
  if(dim(var.diff.list)[1]==0) stop("The differential table have 0 rows")
  if(is.null(threshold)) tab<-var.diff.list[,c(1,5:(dim(var.diff.list)[2]))]
    else{
      if(threshold=="pvalue")tab<-var.diff.list[var.diff.list[,3]<=thr.value,c(1,5:(dim(var.diff.list)[2]))]
      if(threshold=="qvalue")tab<-var.diff.list[var.diff.list[,4]<=thr.value,c(1,5:(dim(var.diff.list)[2]))]
      if(dim(tab)[1]==0) stop("Threshold is too low, the differential table remains with 0 rows")
    }
  rownames(tab)<-tab[,1]
  return(tab[,-1])
}

####################################################################
var.list<-function(expr,labels,FUN){
  mat<-cbind(labels,expr)
  media<-aggregate(mat,by = list(labels),FUN = FUN)[-1]
  tab<-t(scale(media[,-1]))#gene.id
  colnames(tab)<-media[,1]
  return(tab)
}

####################################################################
plot.names<-function(var.diff.list,threshold,thr.value){
  if(dim(var.diff.list)[1]==0) stop("The differential table have 0 rows")
    if(is.null(threshold)) names<-var.diff.list[,1]
    else{
      if(threshold=="pvalue")names<-var.diff.list[var.diff.list[,3]<=thr.value,1]
      if(threshold=="qvalue")names<-var.diff.list[var.diff.list[,4]<=thr.value,1]
      if(length(names)==0) stop("Threshold is too restrict, the differential table remains with 0 rows")
    }
  return(names)
}

####################################################################
#' Structural measures of vertices view in metabolic pathways
#' @description Vertices centralities or clustering coefficient view in KEGG metabolic pathways.
#' centralityPathPlot and pathplot are functions based on pathview function of Pathview package. Pathview is a tool set for pathway based data integration and visualization. It maps and renders user data on relevant pathway graphs. All users need is to supply their gene or compound data and specify the target pathway. Pathview automatically downloads the pathway graph data, parses the data file, maps user data to the pathway, and render pathway graph with the mapped data. Pathview generates both native KEGG view and Graphviz views for pathways. keggview.native and keggview.graph are the two viewer functions, and pathview is the main function providing a unified interface to downloader, parser, mapper and viewer functions.
#' @param gene.data: an output dataframe from diffNetAnalysis function. Data frame structure has genes as rows and statistical test, Nominal p-value, Q-value (p-value FDR adjust for multiple tests)and networks measures, for each network, as columns. Row names should be gene IDs. Here gene ID is a generic concepts, including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs. KEGG ortholog IDs are also treated as gene IDs as to handle metagenomic data. Check details for mappable ID types. Default gene.data=NULL.
#' numeric, character, continuous
#' @param cpd.data: the same as gene.data, excpet named with IDs mappable to KEGG compound IDs. Over 20 types of IDs included in CHEMBL database can be used here. Check details for mappable ID types. Default cpd.data=NULL. Note that gene.data and cpd.data can't be NULL simultaneously.
#' @param threshold: a character indicating which column has to be used to filter which genes or coumponds will be drawn in metabolic map. The options are "pvalue" or "qvalue" to filter by Nominal p-value or Q-value (p-value FDR adjust for multiple tests), respectively. The default threshold=NULL, do not filter any row of data frame.
#' @param thr.values: a numeric value indicating the upper threshold value to filter data frame rows.
#' @param species:character, either the kegg code, scientific name or the common name of the target species. This applies to both pathway and gene.data or cpd.data. When KEGG ortholog pathway is considered, species="ko". Default species="hsa", it is equivalent to use either "Homo sapiens" (scientific name) or "human" (common name).
#' @param pathway.id:character vector, the KEGG pathway ID(s), usually 5 digit, may also include the 3 letter KEGG species code.
#' @param kegg.nattive:logical, whether to render pathway graph as native KEGG graph (.png) or using graphviz layout engine (.pdf). Default kegg.native=TRUE.
#' @param file.name:character, the suffix to be added after the pathway name as part of the output graph file. Sample names or column names of the gene.data or cpd.data are also added when there are multiple samples. Default out.suffix="pathview".
#' @details This function uses pathview to visualize the vertex structural measures in metabolic maps. Pathview maps and renders user data on relevant pathway graphs. Pathview is a stand alone program for pathway based data integration and visualization. It also seamlessly integrates with pathway and functional analysis tools for large-scale and fully automated analysis. Pathview provides strong support for data Integration. It works with: 1) essentially all types of biological data mappable to pathways, 2) over 10 types of gene or protein IDs, and 20 types of compound or metabolite IDs, 3) pathways for over 2000 species as well as KEGG orthology, 4) varoius data attributes and formats, i.e. continuous/discrete data, matrices/vectors, single/multiple samples etc. To see mappable external gene/protein IDs do: data(gene.idtype.list), to see mappable external compound related IDs do: data(rn.list); names(rn.list). Pathview generates both native KEGG view and Graphviz views for pathways. Currently only KEGG pathways are implemented. Hopefully, pathways from Reactome, NCI and other databases will be supported in the future.
#' @return From viersion 1.9.3, pathview can accept either a single pathway or multiple pathway ids. The result returned by pathview function is a named list corresponding to the input pathway ids. Each element (for each pathway itself is a named list, with 2 elements ("plot.data.gene", "plot.data.cpd"). Both elements are data.frame or NULL depends on the corresponding input data gene.data and cpd.data. These data.frames record the plot data for mapped gene or compound nodes: rows are mapped genes/compounds, columns are:
#' kegg.names
#' standard KEGG IDs/Names for mapped nodes. It's Entrez Gene ID or KEGG Compound Accessions.
#' labels
#' Node labels to be used when needed.
#' all.mapped
#' All molecule (gene or compound) IDs mapped to this node.
#' type
#' node type, currently 4 types are supported: "gene","enzyme", "compound" and "ortholog".
#' x
#' x coordinate in the original KEGG pathway graph.
#' y
#' y coordinate in the original KEGG pathway graph.
#' width
#' node width in the original KEGG pathway graph.
#' height
#' node height in the original KEGG pathway graph.
#' other columns
#' columns of the mapped gene/compound data and corresponding pseudo-color codes for individual vertex measures
#' The results returned by keggview.native and codekeggview.graph are both a list of graph plotting parameters. These are not intended to be used externally.
#' @references This function is an adaptation of Luo, W. and Brouwer, C., Pathview: an R/Bioconductor package for pathway based data integration and visualization. Bioinformatics, 2013, 29(14): 1830-1831, doi: 10.1093/bioinformatics/btt285
#' @export
centralityPathPlot<- function(gene.data=NULL, cpd.data=NULL, threshold=NULL, thr.value=0.05, species , pathway.id, kegg.native=T, file.name="path",
                              limit = list(gene = max.gene, cdp = max.cpd), bins = list(gene = 15,cpd = 15), both.dirs= list(gene = F,cpd = F),
                              mid =list(gene = "white", cpd = "white"),high = list(gene = "red",cpd = "red")){

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
  if(is.null(cpd.data) & is.null(gene.data)) stop("You have to insert data in gene.data or cpd.data")
  pv.out <- pathview(gene.data=tab.gene,cpd.data = tab.cpd, pathway.id = pathway.id,
                     species = species, out.suffix = file.name, kegg.native = kegg.native,
                     limit = limit, bins = bins, both.dirs= both.dirs, mid =mid,high = high)
  return(pv.out)
}

####################################################################
#' Variable values view in metabolic pathways
#' @description Variable values view in KEGG metabolic pathways
#' @param gene.data: either vector (single sample) or a matrix-like data (multiple sample). Vector should be numeric with gene IDs as names or it may also be character of gene IDs. Character vector is treated as discrete or count data. Matrix-like data structure has genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a generic concepts, including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs. KEGG ortholog IDs are also treated as gene IDs as to handle metagenomic data. Check details for mappable ID types. Default gene.data=NULL.
#' numeric, character, continuous
#' @param cpd.data: the same as gene.data, excpet named with IDs mappable to KEGG compound IDs. Over 20 types of IDs included in CHEMBL database can be used here. Check details for mappable ID types. Default cpd.data=NULL. Note that gene.data and cpd.data can't be NULL simultaneously.
#' @param labels:a vector of -1s, 0s, and 1s associating each sample with a phenotype. The value 0 corresponds to the first phenotype class of interest, 1 to the second phenotype class of interest, and -1 to the other classes, if there are more than two classes in the gene expression data.
#' @param varr.diff.list: an output dataframe from diffNetAnalysis function. Data frame structure has genes as rows and statistical test, Nominal p-value, Q-value (p-value FDR adjust for multiple tests)and networks measures, for each network, as columns. Row names should be gene IDs. Here gene ID is a generic concepts, including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs.
#' @param threshold: a character indicating which column of "varr.diff.list" has to be used to filter which genes or coumponds will be drawn in metabolic map. The options are "pvalue" or "qvalue" to filter by Nominal p-value or Q-value (p-value FDR adjust for multiple tests), respectively. The default threshold=NULL, do not filter any row of data frame.
#' @param thr.values: a numeric value indicating the upper threshold value to filter data frame rows.
#' @param FUN: a function to define what value will be used in metabolic map.
#' @param species:character, either the kegg code, scientific name or the common name of the target species. This applies to both pathway and gene.data or cpd.data. When KEGG ortholog pathway is considered, species="ko". Default species="hsa", it is equivalent to use either "Homo sapiens" (scientific name) or "human" (common name).
#' @param pathway.id:character vector, the KEGG pathway ID(s), usually 5 digit, may also include the 3 letter KEGG species code.
#' @param kegg.nattive:logical, whether to render pathway graph as native KEGG graph (.png) or using graphviz layout engine (.pdf). Default kegg.native=TRUE.
#' @param file.name:character, the suffix to be added after the pathway name as part of the output graph file. Sample names or column names of the gene.data or cpd.data are also added when there are multiple samples. Default out.suffix="pathview".
#' @details Pathview maps and renders user data on relevant pathway graphs. Pathview is a stand alone program for pathway based data integration and visualization. It also seamlessly integrates with pathway and functional analysis tools for large-scale and fully automated analysis. Pathview provides strong support for data Integration. It works with: 1) essentially all types of biological data mappable to pathways, 2) over 10 types of gene or protein IDs, and 20 types of compound or metabolite IDs, 3) pathways for over 2000 species as well as KEGG orthology, 4) varoius data attributes and formats, i.e. continuous/discrete data, matrices/vectors, single/multiple samples etc. To see mappable external gene/protein IDs do: data(gene.idtype.list), to see mappable external compound related IDs do: data(rn.list); names(rn.list). Pathview generates both native KEGG view and Graphviz views for pathways. Currently only KEGG pathways are implemented. Hopefully, pathways from Reactome, NCI and other databases will be supported in the future.
#' @return From viersion 1.9.3, pathview can accept either a single pathway or multiple pathway ids. The result returned by pathview function is a named list corresponding to the input pathway ids. Each element (for each pathway itself is a named list, with 2 elements ("plot.data.gene", "plot.data.cpd"). Both elements are data.frame or NULL depends on the corresponding input data gene.data and cpd.data. These data.frames record the plot data for mapped gene or compound nodes: rows are mapped genes/compounds, columns are:
#' kegg.names
#' standard KEGG IDs/Names for mapped nodes. It's Entrez Gene ID or KEGG Compound Accessions.
#' labels
#' Node labels to be used when needed.
#' all.mapped
#' All molecule (gene or compound) IDs mapped to this node.
#' type
#' node type, currently 4 types are supported: "gene","enzyme", "compound" and "ortholog".
#' x
#' x coordinate in the original KEGG pathway graph.
#' y
#' y coordinate in the original KEGG pathway graph.
#' width
#' node width in the original KEGG pathway graph.
#' height
#' node height in the original KEGG pathway graph.
#' other columns
#' columns of the mapped gene/compound data and corresponding pseudo-color codes for individual samples
#' The results returned by keggview.native and codekeggview.graph are both a list of graph plotting parameters. These are not intended to be used externally.
#' @references Luo, W. and Brouwer, C., Pathview: an R/Bioconductor package for pathway based data integration and visualization. Bioinformatics, 2013, 29(14): 1830-1831, doi: 10.1093/bioinformatics/btt285
#' @export
pathPlot<- function(gene.data=NULL, cpd.data=NULL, labels, varr.diff.list=NULL, threshold=NULL, thr.value=0.05, FUN=median,species , pathway.id, kegg.native=T, file.name="path"){

    if(!is.null(gene.data)) {
      tab.gene<-var.list(gene.data,labels=labels,FUN=FUN)
      names.gene<-plot.names(gene.data,threshold=threshold,thr.value=thr.value)
      # max.gene<-max(tab.gene)
    }
    else{
      tab.gene<-NULL
      names.gene<-NULL
      # max.gene<-1
    }
    if(!is.null(cpd.data)){
      tab.cpd<-var.list(cpd.data,labels=labels,FUN=FUN)
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
