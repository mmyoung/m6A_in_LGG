
##This script contains funtions: call_featurecounts, call_DESeq2, anno_DESeq2

## Usage: Rscript featurecounts.R out_dir bam_dir spices labels
## Spices must be one of "human" or "mouse"
## labels would be the group label of all samples by the order in hisat_results folder

library(stringr)
base_dir <- "/data/nym/20210339_LGG/"
output_dir <- "/data/nym/20210339_LGG/all_call_peak/mergeBed_res/"
call_featurecounts <- function(out_dir,bam_dir,spices){
	out_dir <- out_dir
	anno_file <- ifelse(spices=="human","/data/nym/reference/reference_hg38/Annotate_85/Homo_sapiens.GRCh38.85.gtf","/data/nym/reference/reference_mouse_mm10/Annotation/Mus_musculus.GRCm38.86.gtf")
	setwd(out_dir)
	bam_list <- list.files("/data/nym/20210339_LGG/hisat2_res",pattern="Q20_sorted_nonrRNA.bam$", full.names=TRUE)
	library(Rsubread)
	fc<-featureCounts(files=bam_list,annot.ext=anno_file,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",isPairedEnd=TRUE,nthreads=4)
	write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],fc$counts,stringsAsFactors=FALSE),file=paste0(out_dir,"trans_counts.txt"),quote=FALSE,sep="\t",row.names=FALSE)
}

call_featurecounts(out_dir=output_dir,paste0(base_dir,"hisat_result"),"human")

## usage: call_DESeq2(input_dir, labels, class)
## return: res
call_DESeq2 <- function(input_dir, labels, class){
	library(DESeq2)
	library(dplyr)
	labels <- labels
	class <- class
	setwd(input_dir)
	raw_data <- read.table(paste0(input_dir,"counts.txt"),header=T,row.names = 1) %>%
			select(contains(c("Length",labels)))
	expr_data <- raw_data[,c(2:dim(raw_data)[2])]
	colnames(expr_data) <- labels
	coldata <- data.frame(row.names=labels, condition=factor(class))
	dds <- DESeqDataSetFromMatrix(countData = expr_data,colData = coldata,design = ~ condition)
	dds <- DESeq(dds)
	res <- results(dds,contrast=c("condition","treated","control"))
	### make sample correlation plot
	rld<- rlogTransformation(dds, blind=TRUE)
	distsRL <- dist(t(assay(rld)))
	mat<- as.matrix(distsRL)
	rownames(mat) <- colnames(mat) <- colnames(expr_data)
	library('RColorBrewer')
	library('gplots')
	hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
	hc <- hclust(distsRL)
	pdf(paste0(input_dir,"hclust_heatmap.pdf"),width=5,height=4)
	heatmap.2(mat, Rowv=as.dendrogram(hc),
        	symm=TRUE, trace='none',
        	col = rev(hmcol), margin=c(9, 13))
	dev.off()
	return(list(res, raw_data))
}

#call_DESeq_return <- call_DESeq2(output_dir,labels,c("control","control","treated","treated"))
#DESeq2_res <- call_DESeq_return[1]
#raw_data <- call_DESeq_return[2][[1]]

## usage: annot_DESeq2(DESeq2_res, raw_data, labels, spices)
annot_DESeq2 <- function(DESeq2_res,raw_data, output_dir,labels, spices){
	DESeq2_res <- as.data.frame(DESeq2_res)
	raw_data <- raw_data
	labels <- labels
	if(spices == "human"){
		anno_file <- read.table("/data/nym/reference/reference_hg38/Annotate_85/GRch38.85_name_type_pos",stringsAsFactors=F, header=F)
	}
	else{
		anno_file <- read.table("/data/nym/reference/reference_mouse_mm10/Annotation/GRCm38.86_name_type_pos",stringsAsFactors=F, header=F)
	}
	colnames(anno_file) <- c("chr","start","end","gene_id","gene_name","gene_type")
	colnames(raw_data) <- c("length",labels)
	for(i in labels){
		j=paste(i,"_rpkm",sep = "")
		raw_data[,j] <- raw_data[,i]*10**9/(colSums(raw_data)[i]*raw_data$length)
	}
	anno_res <- merge(raw_data, anno_file, by.x=0,by.y="gene_id")
	anno_res_2 <- merge(anno_res, DESeq2_res, by.x="Row.names",by.y = 0)
	write.table(anno_res_2,file = paste0(output_dir,"DESeq2_anno.txt"),sep = "\t", row.names = F,quote = F)
}

#annot_DESeq2(DESeq2_res, raw_data, output_dir, labels, "human")

