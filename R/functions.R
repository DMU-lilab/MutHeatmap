require(data.table)
require(RColorBrewer)

#' Mutation heatmap with multiple mutation types
#'
#' Plot a heat map with multiple mutation types within a gene
#'
#' The gene heat map that we usually see is a grid of genes, a type of mutation,
#' but in fact the same gene often has multiple mutation types in the same patient,
#' so the traditional heat map drawing tool can't satisfy us. For the problem, this
#'  function uses image to draw the preliminary heat map, and then uses points to
#' add the second mutation in the form of a square. The third mutation is added in
#' turn. At the same time, the square position is slightly moved and accompanied by
#' a slight reduction in size to achieve a better display effect, and up to four
#' mutations can be represented on one heat map grid.
#'
#' Before run MutHeatmap, one should known:
#'
#'    1. If you want to draw a heat map with a histogram, because the image size is too large, use the pdf
#'function and give it a large enough width and length;
#'
#'    2. The default is the type of mutation annotated by annovar;
#'
#'    3. Because there are too many factors affecting the alignment of the heat map and the histogram during
#'drawing, it is difficult to adjust the corresponding mar, mex, and oma parameters to achieve better
#'results. Therefore, it is recommended to quickly draw a rough estimate and then use inkscape or
#'adobe for layout alignment.
#'
#'   4. If the point of the mutation type in the heat map is too small, the width and length of the pdf file
#'should be reduced.
#'
#' @param vr The data frame containing the variance data, there are three columns,
#'  the column names should be: Gene, Type, Patient
#' @param pal Character vectors, Need to be customized according to the number of
#'  mutation types in the data frame, need more than the mutation type one color, it
#'  is used as the background color, and the background color is placed in the
#'  first place
#' @param type Character vectors. The corresponding mutation type, type have to equal
#'  to or more than all types appearing in the data; by default, the four types of copy
#'  number plus neutral plus annovar mutations. Default: "DEL","LOSS","NEUTRAL","GAIN",
#'  "AMPL","nonsynonymous SNV","synonymous SNV","intronic","stopgain","nonframeshift
#'  deletion", "splicing", "frameshift deletion","UTR3","frameshift insertion","UTR5"
#'  mutations annovar; the length of type is less  1 than pal
#' @param order_gene A logical scalar. Default T, sorted according to the number of
#'  patients gene mutations
#' @param order_patient A logical scalar. Default T, according to the patient The number of mutated
#'  genes is sorted
#' @param hist_plot A logical scalar.Default T, plus the corresponding histogram above and to the right
#' @param legend_dist Numeric. Default 0.4, adjust the distance between the legends, generally need to
#'  adjust
#' @param col_text_cex Numeric. Adjust the size of the patient name, default. 1
#' @param row_text_cex Numeric. Gene name adjust the size of the default. 1
#' @param xlab_adj Numeric. Patient name and adjustment FIG distance between
#' @param sub_gene Numeric. Character vectors. Drawing only selected portions of genes
#'  and the gene is present in the data needed, the default NULL
#' @param heatmap_mar Numeric. Mar parameter, adjusting the length of the left and right edges
#'  of the front and rear heat FIG default C ( 5,17,1,2)
#' @param Heatmap_oma Numeric. Oma parameter, adjust the outer edge length of the heat map before
#'  and after, the default c (0.2, 0.2, 0.2, 0.2)
#' @param mex Numeric. Adjust the heat map's mex parameter, used to describe the coordinates of
#'  the drawing edge, the default 0.5
#' @param legend_mar Numeric. Legend Mar parameter, adjust the position of the legend, default
#'  c (1,0,4,1)
#' @param order_omit Character vectors. The type of variation that is ignored when sorting, these mutation
#'  types will also be filtered in the histogram, the default c ("NEUTRAL"), If there is
#'   no mutation type "NEUTRAL", you can also keep the default parameters
#' @param heatmap_height Numeric. The height of the heat map. For the case of a large number of
#'  genes, based on the length of the pdf, you can set the heatmap to cover the height of
#'   the canvas. Grid
#' @param heatmap_width Numeric. The width of the heat map. For a very large number of patients,
#'  you can set the length of the heat map
#' @param Annotation_col Data.frame. Information annotating the patient,  the line name
#'  is the patient, the column name is different pathological information, the data must be
#'   a factor, rownames should be patients. The example is as follows (in the code example I briefly SmokingInfo is Smoking,
#'    OldInfo is Old)
#'        SmokingInfo   OldInfo
#' P1     Smoking   Old
#' P2  NonSmoking Unold
#' P3  NonSmoking Unold
#' P4  NonSmoking   Old
#' P5     Smoking Unold
#' P6  NonSmoking Unold
#' @param Annotation_colors Character vectors list. Match different color to pathological information, examples are
#' as follows, need to match the column name and factor in annotation_col
#' $SmokingInfo
#' Smoking NonSmoking
#' "#FDB462"  "#80B1D3"
#' $OldInfo
#' Old     Unold
#' "#FB8072" "#8DD3C7"
#' @param Anno_height Sets the height of the annotation, which is automatically adjusted by default.
#'
#' @examples
#' #without hist plot
#' pdf("heatmap_cnv_mut.pdf", height=12, width = 12)
#' MutHeatmap(vr, heatmap_mar = c(17,17,1,2),hist_plot = F, legend_dist=0.1, xlab_adj = 1.2, order_patient = T, order_gene = T)
#' dev.off()
#'
#' #with hist plot
#' pdf("heatmap_hist_cnv_mut.pdf", height=12, width = 12)
#' MutHeatmap(vr, heatmap_mar = c(17,7,1,2),hist_plot = T, legend_dist=0.3, xlab_adj = 1.2, order_patient = T, order_gene = T)
#' dev.off()
#'
#' #only a few gene
#' pdf("Assoc_CN1.pdf", height=2,width = 14)
#' # have to adjust legend by manual
#' MutHeatmap(vr, heatmap_mar = c(7,17,1,2), sub_gene = c("CDKN2A", "GNAQ", "NOTCH1", "RB1", "SMAD4", "ABL1"),hist_plot = F,legend_dist=0.2, xlab_adj = 0.9)
#' dev.off()
#'
#' #with annotation and hist
#' annotation_col <- data.frame(Smoking = factor(sample(c("Smoking", "NonSmoking"), length(unique(vr$Patient)), replace = T)), Old=factor(sample(c("Old", "Unold"), length(unique(vr$Patient)), replace = T)))
#' rownames(annotation_col) <- unique(vr$Patient)
#' annotation_colors <- list(Smoking =c(Smoking = "#FDB462", NonSmoking = "#80B1D3"), Old=c(Old = "#FB8072", Unold = "#8DD3C7"))
#' pdf("heatmap_hist_cnv_mut.pdf", height=15, width = 12)
#' MutHeatmap(vr, heatmap_mar = c(15,10,1,2), legend_mar = c(1,0,1,1), hist_plot = T, legend_dist=0.2, xlab_adj = 1.2, order_patient = T, order_gene = T, annotation_col=annotation_col, annotation_colors=annotation_colors)
#' dev.off()
#'
#' #with annotation and without hist
#' pdf("heatmap_cnv_mut.pdf", height=12, width = 12)
#' MutHeatmap(vr, heatmap_mar = c(17,10,1,2),hist_plot = F, legend_dist=0.1, xlab_adj = 1.2, order_patient = T, order_gene = T, annotation_col=annotation_col, annotation_colors=annotation_colors)
#' dev.off()
#'
#' #Start with the first step(there is not test data in MutHeatmap package)
#' # read variant data
#' dt <- read.csv("~/Downloads/selected_genetic_2.csv", sep = ";")
#' setDT(dt)
#' setnames(dt, c("Patient", "Gene", "Type")) # change column name
#' # read clinical data
#' ClinicalInfo <- read.csv("~/Downloads/selected_tma_clinic_2.csv", sep = ";")
#' annotation_col <- ClinicalInfo[,-1]
#' rownames(annotation_col) <- ClinicalInfo[,1]
#' annotation_colors <- list(Age = c("<50" = "#FDB462", ">50" = "#80B1D3"), Smoking=c("yes" = "#FB8072", "no" = "#8DD3C7"),PDL1_status = c("pos" = "#3288BD", "neg" = "#5E4FA2") )
#' # plot heatmap
#' pdf("example_4.pdf", height=15, width = 12)
#' MutHeatmap(dt, pal = c("#F2F2F2",brewer.pal(8, "Paired")) ,type =
#'             c("missense", "nonsense", "frameshiftDeletion", "frameshiftInsertion", "nonframeshiftDeletion", "Amplification", "Deletion", "Fusion"),order_gene = T, order_patient = T, hist_plot =
#'             T, legend_dist = 0.4, col_text_cex = 1, row_text_cex = 1, sub_gene= NULL,heatmap_mar = c(20,17,1,2), heatmap_oma=c(0.2,0.2,0.2,0.2),heatmap_mex=0.5, legend_mar = c(1,0,4,1),xlab_adj=0.2,
#'           order_omit=c("NEUTRAL"), annotation_col = annotation_col, annotation_colors = annotation_colors, heatmap_height = 5,
#'           heatmap_width = 6, anno_height = 0.5)
#' dev.off()
#'
#' @export MutHeatmap
#' @import data.table
#' @import RColorBrewer


MutHeatmap <- function(vr, pal = c("#F2F2F2",colorRampPalette(c("blue", "white", "red"))(5)[c(1,2)],"#F2F2F2",colorRampPalette(c("blue", "white", "red"))(5)[c(4,5)],brewer.pal(n = 8, name ="Accent")[c(1,4,6,8,2,3,5,7)],"#E31A1C","#6A3D9A"),type = c("DEL","LOSS","NEUTRAL","GAIN","AMPL","nonsynonymous SNV","synonymous SNV","intronic","stopgain","nonframeshift deletion","splicing", "frameshift deletion","UTR3","frameshift insertion","UTR5"),order_gene = T, order_patient = T, hist_plot = T, legend_dist = 0.4, col_text_cex = 1, row_text_cex = 1, sub_gene= NULL,heatmap_mar = c(5,17,1,2), heatmap_oma=c(0.2,0.2,0.2,0.2),heatmap_mex=0.5, legend_mar = c(1,0,4,1),xlab_adj=1, order_omit=c("NEUTRAL"), annotation_col=NULL, annotation_colors = NULL, heatmap_height = 3, heatmap_width = 3, anno_height=NULL)
{
  if(length(pal) > 18){stop("Error! max number of mutation type is 17!")}
  if((length(pal) - length(type)) !=1 ){stop("Error! Pal must be one longer than type, because first one pal is col for no mutation")}
  if(sum(unique(vr$Type) %in%  type) != length(unique(vr$Type))){stop("Error! Some mutations in data don't have corresponding color")}
  if(!is.null(annotation_col) & sum(colnames(annotation_col) %in% names(annotation_colors)) != length(colnames(annotation_col))){stop("Error! The annotation_colors are not match to annotation_col")}
  if(!is.null(sub_gene)){
    pal_dt <- data.table(pal, type=c("NoMut",type))
    vr <- vr[Gene %in% sub_gene,]
    type <- pal_dt[type %in% unique(vr$Type),type]
    pal <- c(pal[1],pal_dt[type, on="type"][,pal])
  }else{
    pal_dt <- data.table(pal, type=c("NoMut",type))
    type <- pal_dt[type %in% unique(vr$Type),type]
    pal <- c(pal[1],pal_dt[type, on="type"][,pal])
  }
  dt <- unique(vr[,.(Gene,Type,Patient)])
  dt$Type <- factor(dt$Type, levels = type)
  if(order_gene){gene <- dt[!Type %in% order_omit,.(N=length(unique(Patient))),by=Gene][order(N),Gene]}else{gene <- unique(dt[!Type %in% order_omit, Gene])}
  dt$Gene <- factor(dt$Gene, levels = gene)
  if(order_patient){patient <- data.table(table(vr[!Type %in% order_omit,]$Patient))[order(-N),V1]}else{patient <- unique(dt[!Type %in% order_omit, Patient])}
  dt$Patient <- factor(dt$Patient, levels = c(patient, setdiff(unique(dt$Patient),patient)))

  setkey(dt, "Type")

  n <- length(unique(dt$Type))

  dt$Gene_Patients <- paste(dt$Gene, dt$Patient)
  dt_inf <- dt[,.N,by=.(Gene, Patient)]
  max_mut_num <- max(dt_inf$N)
  dt[,Mut_num:=seq_len(.N),by=.(Patient,Gene)]
  #main plot
  dt1 <- copy(dt)
  dt1[Mut_num !=1, Type:=NA]
  dc <- data.frame(dcast(dt1, Patient ~ Gene, value.var = "Type", fun.aggregate = function(x)(x[!is.na(x)][1])))
  rownames(dc)<- dc[,1]
  data_matrix<-data.matrix(dc[,-1])
  data_matrix[is.na(data_matrix)] <- 0
  pal=pal
  breaks<-seq(-1,17,1)
  if(!hist_plot & is.null(annotation_col)){
    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(8,2), heights=c(1,1))
    par(mar=heatmap_mar, oma=heatmap_oma, mex=heatmap_mex)
  }else if(hist_plot & is.null(annotation_col)){
    layout(matrix(c(2,4,1,3),2,2,byrow=TRUE), widths=c(heatmap_width,1),
           heights=c(1,heatmap_height), TRUE)
    par(mar=heatmap_mar)
  }else if(hist_plot & !is.null(annotation_col)){
    if(is.null(anno_height)){anno_height <- 0.02 * ncol(annotation_col)}
    layout(matrix(data=c(3,5,2,5,1,4), nrow=3, ncol=2, byrow=TRUE), widths=c(heatmap_width,1), heights=c(1, anno_height, heatmap_height))
    par(mar=heatmap_mar, oma=heatmap_oma, mex=heatmap_mex)
  }else if(!hist_plot & !is.null(annotation_col)){
    if(is.null(anno_height)){anno_height <- 0.02 * ncol(annotation_col)}
    layout(matrix(data=c(2,4,1,3), nrow=2, ncol=2, byrow=TRUE), widths=c(8,2), heights=c(anno_height,1))
    par(mar=heatmap_mar, oma=heatmap_oma, mex=heatmap_mex)
  }




  image(x=1:nrow(data_matrix),y=1:ncol(data_matrix),
        z=data_matrix,xlab="",ylab="",
        breaks = breaks,
        col=pal[1:18],axes=FALSE)

  #sub plot
  add_plot <- function(dt, i){
    dt1 <- copy(dt)
    dt1[Mut_num != i, Type:=NA]
    dc <- data.frame(dcast(dt1, Patient ~ Gene, value.var = "Type", fun.aggregate = function(x){ifelse(length(x) >1,x[!is.na(x)][1],factor(NA))}))
    rownames(dc)<- dc[,1]
    data_matrix <- data.matrix(dc[,-1])
    xy <- which(data_matrix !=0, arr.ind = T)
    for(r in 1:nrow(xy)){x <- xy[r,]; points(x[1]-0.6+i*0.25, x[2],pch=15, cex=1.2 - i*0.08, col=pal[data_matrix[x[1],x[2]]+1])}
  }

  ploti <- 2:max_mut_num
  for(i in ploti){add_plot(dt, i)}

  text(x=1:nrow(data_matrix)+0.1, y=par("usr")[1] - xlab_adj,
       srt = 90, adj = 0.5, labels = rownames(data_matrix),
       xpd = TRUE, cex=col_text_cex)
  axis(2,at=1:ncol(data_matrix),labels=colnames(data_matrix),
       col="white",las=1, cex.lab=0.1, cex.axis=row_text_cex)
  abline(h=c(1:ncol(data_matrix))+0.5,v=c(1:nrow(data_matrix))+0.5,
         col="white",lwd=2,xpd=F)

  #add annotation plot
  if(!is.null(annotation_col)){
    for (i in 1:ncol(annotation_col)){
      if(!is.factor(annotation_col[,i])){
        warning(paste0(names(annotation_col)[i], " is not factor, now transforming to factor!"))
        annotation_col[,i] <- as.character(annotation_col[,i])
        if(length(which(is.na(annotation_col[,i]))) > 0){
          warning(names(annotation_col)[i], " has NA value, now transforming to factor named NA!")
          annotation_col[,i][is.na(annotation_col[,i])] <- "NA"
        }
        annotation_col[,i] <- factor(annotation_col[,i])
      }
    }
    change_factor <- function(x){k <- 1:length(unique(x)); names(k) <- as.character(unique(x)); k[as.character(x)]} #change infomation to numeric
    colname_tmp <- colnames(annotation_col)
    #annotation_col_mt <- as.matrix(apply(annotation_col, 2, change_factor))
    annotation_col_mt <- NULL;for (i in 1:ncol(annotation_col)){annotation_col_mt <- cbind(annotation_col_mt, (change_factor(annotation_col[,i])))}
    rownames(annotation_col_mt) <- rownames(annotation_col)
    colnames(annotation_col_mt) <- colname_tmp
    tryCatch(annotation_col_mt <- annotation_col_mt[rownames(data_matrix),],error = function(e){stop("Error! All the patients in the data should have corresponding information in annotation_col")})
    ## change infomation numric to unique number
    cumsum <- 0
    if(!is.null(dim(annotation_col_mt))){#if more than one column, cummulate info numeric
      for(i in 1:ncol(as.data.frame(annotation_col_mt))){
        annotation_col_mt[,i] <- cumsum + annotation_col_mt[,i]
        cumsum <- max(annotation_col_mt[,i])
      }
    }
    ## get color according to infomation
    get_color <- function(anno){return(annotation_colors[[anno]][levels(annotation_col[,anno])])}
    palAnn <- NULL
    if(is.null(dim(annotation_col_mt))){
      rowname_tmp <- rownames(annotation_col_mt)
      annotation_col_mt <- as.matrix(annotation_col_mt, nrow=1)
      rownames(annotation_col_mt) <- rowname_tmp
      colnames(annotation_col_mt) <- colname_tmp
    }
    for(anno in colnames(annotation_col_mt)){
      palAnn <- c(palAnn, get_color(anno))
    }
    par(mar=c(0,heatmap_mar[2], 0,  heatmap_mar[4]))
    image(x=1:nrow(annotation_col_mt), y=1:ncol(annotation_col_mt), z= annotation_col_mt, col=palAnn, xlab="",ylab="",
          axes=FALSE)
    axis(2,at=1:ncol(annotation_col_mt),labels=colnames(annotation_col_mt),
         col="white",las=1, cex.lab=0.1, cex.axis=row_text_cex)
    abline(h=c(1:ncol(annotation_col_mt))+0.5,v=c(1:nrow(annotation_col_mt))+0.5,
           col="white",lwd=2,xpd=F)
  }
  if(hist_plot){
    #hist
    par(mar=c(0,2+0.5,3,heatmap_mar[4]-0.9))
    patient_dt <- dt[,.N,by=.(Patient,Type)]
    mt <- data.frame(dcast(patient_dt, Type ~ Patient, value.var = "N"))
    data_matrix <- data.matrix(mt[,-1])
    rownames(data_matrix) <- mt[,1]
    tryCatch(data_matrix <- data_matrix[setdiff(type, order_omit), patient], error = function(e){stop("Error! Type argument or your patient name format(include "-" and so on )")})
    data_matrix[is.na(data_matrix)] <- 0
    omit_idx <- NULL
    for(i in order_omit){omit_idx <- c(omit_idx,1+which(type == i))}
    barplot(data_matrix, col=pal[-c(1,omit_idx)],space=0,border = "white",axes=T,xlab="",ann=F, xaxt="n")

    par(mar=c( heatmap_mar[1]-2 , 0.8, heatmap_mar[3]+2.2, 3),las=1)
    gene_dt <- dt[,.N,by=.(Gene,Type)]
    mt <- data.frame(dcast(gene_dt, Type ~ Gene, value.var = "N"))
    data_matrix <- data.matrix(mt[,-1])
    rownames(data_matrix) <- mt[,1]
    gene <- gsub("ATM,", "ATM.", gene)
    tryCatch(data_matrix <- data_matrix[setdiff(type, order_omit), gene], error = function(e){stop("Error! Type argument or check your gene name format(please not include "-" and so on)")})
    data_matrix[is.na(data_matrix)] <- 0
    barplot(data_matrix, col=pal[-c(1,omit_idx)],space=0,border = "white",axes=T,xlab="", ann=F, horiz = T, yaxt="n")

  }

  #add legend
  par(mar=legend_mar)
  plot(3, 8,  axes=F, ann=F, type="n")
  if(is.null(annotation_col)){
    ploti <- 1:length(type)
  }else{
    #add annotation legend
    ploti <- 1:(length(type) + max(annotation_col_mt))
    pal <- c("NULL", palAnn, pal[-1])
    anno_label <- NULL
    for (anno in colnames(annotation_col)){
      anno_label <- c(anno_label, levels(annotation_col[[anno]]))
    }
    type <- c(anno_label,type)
  }
  if(!hist_plot){
    tmp <- for(i in ploti){points(2, 10+(length(type)-i)*legend_dist, pch=15, cex=2, col=pal[i+1])}
    tmp <- for(i in ploti){text(3, 10+(length(type)-i)*legend_dist, labels = type[i],pch=15, cex=1, col="black")}
  }
  if(hist_plot){
    tmp <- for(i in ploti){points(2, 5+(length(type)-i)*legend_dist, pch=15, cex=0.9, col=pal[i+1])}
    tmp <- for(i in ploti){text(2.8, 5+(length(type)-i)*legend_dist, labels = type[i],pch=15, cex=0.9, col="black")}
  }

}
