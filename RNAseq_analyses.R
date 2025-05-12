## Code for conducting differential expression analyses on bulk RNAseq data using DESeq2 and limma ##
## This code has been developed by Jon Sanchez-Valle (BSC) and Andrea Marti Sarrias (UB) ##
## QC, metadata and counts are downloaded from GREIN, counts are also downloaded from ARCHS4 for comparative analyses ##

## Load the needed libraries ##
library("data.table")
library("DESeq2")
library("lsa")
library("dplyr")
library("PCAtools")
library("sva")
library("edgeR")
library("EnhancedVolcano")
library("VennDiagram")
library("SummarizedExperiment")
library("devtools")
library("calibrate")
library("dendextend")

## Provide and argument
args = commandArgs(trailingOnly=TRUE)

## Preliminary analysis on the quality of the sampels ##
if(args[1]=="check_quality_of_samples"){
  if("Aligned_reads"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/Aligned_reads")}
  ## List the tables donwloaded from GREIN ##
  qcs<-list.files("GREIN/QC/")
  ## Extract the number of samples with at least 70% of mapped reads ##
  qbysample<-c() ; qcsmedios<-c()
  for(a in 1:length(qcs)){
    qctab<-fread(paste("GREIN/QC/",qcs[a],sep=""),stringsAsFactors = F,sep="\t")
    qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
    qbysample<-qctab[,1:3]
    qbysample[[1]]<-gsub("_transcripts_quant","",qbysample[[1]])
    write.table(qbysample,paste("GREIN/Aligned_reads/",qcs[a],sep=""),quote=F,sep="\t",row.names = F)
    qcsmedios<-rbind(qcsmedios,c(length(which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)),length(qctab$`% Aligned`),
                                 round((length(which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70))/length(qctab$`% Aligned`))*100,3),mean(as.numeric(gsub("%","",qctab$`% Aligned`)))))
  }
  rownames(qcsmedios)<-gsub(".txt","",qcs)
  colnames(qcsmedios)<-c("morethan70","total_samples","percentage","mean_percentage_aligned_reads")
  ## Select only those studies with at least 10 samples with a percentage of aligned reads larger than 70% ## 
  workingfiles<-qcsmedios[which(qcsmedios[,1]>=10),]
  ## Check if we have metadata and counts for all the selected studies ##
  write.table(qcsmedios,"GREIN/Summary_by_study.txt",quote=F,sep="\t")
  write.table(workingfiles,"GREIN/Summary_by_study_working_datasets.txt",quote=F,sep="\t")
}


## Create the Summarized Experiment file with the samples with enough quality (at least 70% of mapped reads), giving the same name to the key columns of the metadata ##
if(args[1]=="create_summarized_experiment"){
  if("Processed_metadata"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/Processed_metadata")}
  if("SummarizedExperiments"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/SummarizedExperiments")}
  if("SummarizedExperiments_CaseControl"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/SummarizedExperiments_CaseControl")}
  metadatas<-list.files("GREIN/Metadata/")
  counts<-list.files("GREIN/Counts/")
  #### AD_brain_GSE95587 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[1],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[1]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[1],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,11,13,14,15,16,17,18)]
  colnames(metadata)<-c("sample","organism","tissue","sex","age","category","characteristic1","characteristic2","characteristic3")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("diagnosis: Alzheimer's disease","case",metadata$category)
  metadata$category<-gsub("diagnosis: control","control",metadata$category)
  metadata$sex<-gsub("Sex: F","female",metadata$sex) ; metadata$sex<-gsub("Sex: M","male",metadata$sex)
  metadata$age<-gsub("age at death: ","",metadata$age)
  metadata$tissue<-gsub("tissue: ","",metadata$tissue)
  metadata<-metadata[,c(1,6,4,5,3,2,7:10)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/AD_brain_GSE95587.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/AD_brain_GSE95587.rds")
  
  #### COVID19_cl_GSE147507 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[2],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[2]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[2],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,13,11,12,14)]
  colnames(metadata)<-c("sample","organism","category","tissue","characteristic1","characteristic2")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("treatment: SARS-CoV-2 infected.+","COVID",metadata$category)
  metadata$category<-gsub("treatment: RSV.+","RSV",metadata$category)
  metadata$category<-gsub("treatment: HPIV3.+","HPIV3",metadata$category)
  metadata$category<-gsub("treatment: Mock treatment","control",metadata$category)
  metadata$category<-gsub("treatment: ","",metadata$category)
  metadata$tissue<-gsub("cell line: ","",metadata$tissue)
  metadata$characteristic1<-gsub("cell line: ","",metadata$characteristic1)
  metadata$characteristic2<-gsub("time point: ","",metadata$characteristic2)
  metadata<-metadata[,c(1,3,4,5,6,2,7)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## consider those infected by COVID as cases
  metadata$category<-gsub("^COVID$","case",metadata$category)
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/COVID19_cl_GSE147507.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/COVID19_cl_GSE147507.rds")
  
  #### COVID19_dopaminergic_GSE174745 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[3],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[3]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[3],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,11,47,48)]
  colnames(metadata)<-c("sample","organism","tissue","treatment","category")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$tissue<-gsub("cell type: ","",metadata$tissue)
  metadata<-metadata[,c(1,5:2,6)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## Focus on non-treated samples
  metadata$category<-gsub("SARS-Cov-2 infected 48 hours","case",metadata$category)
  metadata$category<-gsub("no","control",metadata$category)
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/COVID19_dopaminergic_GSE174745.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/COVID19_dopaminergic_GSE174745.rds")
  
  #### COVID19_hipsc_GSE179923 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[4],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[4]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[4],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,9,10,12)]
  colnames(metadata)<-c("sample","tissue","organism","category")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("treatment: SARS-CoV-2","case",metadata$category)
  metadata$category<-gsub("treatment: Mock","control",metadata$category)
  metadata$category<-gsub("treatment: ","",metadata$category)
  metadata<-metadata[,c(1,4,2,3,5)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/COVID19_hipsc_GSE179923.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/COVID19_hipsc_GSE179923.rds")
  
  #### COVID19_pbmc_GSE152418 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[5],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[5]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[5],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,9,12,11,13,14)]
  colnames(metadata)<-c("sample","organism","tissue","sex","days_post_symptom_onset","category","severity")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("disease state: COVID-19","case",metadata$category)
  metadata$category<-gsub("disease state: Healthy","control",metadata$category)
  metadata$sex<-gsub("gender: F","female",metadata$sex) ; metadata$sex<-gsub("gender: M","male",metadata$sex)
  metadata$days_post_symptom_onset<-gsub("days_post_symptom_onset: ","",metadata$days_post_symptom_onset)
  metadata$severity<-gsub("severity: ","",metadata$severity)
  metadata<-metadata[,c(1,6,4,5,3,2,7,8)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/COVID19_pbmc_GSE152418.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/COVID19_pbmc_GSE152418.rds")
  
  #### COVID19_pbmc_GSE251849 - has Long COVID, Long COVID brain fog, and COVID convalescent in the detailed_category ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[6],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[6]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[6],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,9,12,13)]
  colnames(metadata)<-c("sample","organism","tissue","sex","category")
  metadata$category<-gsub("group: ","",metadata$category)
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("Healthy","control",metadata$category)
  metadata$category[which(metadata$category!="control")]<-"case"
  metadata$sex<-gsub("Sex: female","female",metadata$sex) ; metadata$sex<-gsub("Sex: male","male",metadata$sex)
  metadata<-metadata[,c(1,5,4,3,2,6)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/COVID19_pbmc_GSE251849.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/COVID19_pbmc_GSE251849.rds")
  
  #### MS_brain_GSE123496 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[7],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[7]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[7],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,9,12)]
  colnames(metadata)<-c("sample","organism","tissue","category")
  metadata$category<-gsub("disease state: ","",metadata$category)
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("MS","case",metadata$category)
  metadata$category<-gsub("healthy control","control",metadata$category)
  metadata<-metadata[,c(1,4,3,2,5)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/MS_brain_GSE123496.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/MS_brain_GSE123496.rds")
  
  #### PD_brain_GSE136666 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[8],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[8]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[8],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,9,11,40)]
  colnames(metadata)<-c("sample","organism","tissue","sex","category")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category<-gsub("Parkinson’s disease","case",metadata$category)
  metadata$category<-gsub("Control","control",metadata$category)
  metadata$sex<-gsub("gender: f","female",metadata$sex) ; metadata$sex<-gsub("gender: m","male",metadata$sex)
  metadata$tissue<-gsub("tissue: ","",metadata$tissue)
  metadata<-metadata[,c(1,5,4,3,2,6)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/PD_brain_GSE136666.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/PD_brain_GSE136666.rds")
  
  #### PD_brain_GSE216281 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[9],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[9]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[9],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,9,48,44,45,46)]
  colnames(metadata)<-c("sample","organism","tissue","sex","age","category","characteristic1")
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category[which(metadata$category>0)]<-"case" ; metadata$category[which(metadata$category==0)]<-"control"
  metadata$sex<-gsub("Female","female",metadata$sex) ; metadata$sex<-gsub("Male","male",metadata$sex)
  metadata<-metadata[,c(1,6,4,5,3,2,7,8)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/PD_brain_GSE216281.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/PD_brain_GSE216281.rds")
  
  #### PD_brain_GSE68719 ####
  ## Select the samples with at least a 70% of aligned reads ##
  qctab<-fread(paste("GREIN/QC/",rownames(workingfiles)[10],".txt",sep=""),stringsAsFactors = F,sep="\t")
  qctab<-qctab[grep("transcripts_quant",qctab$`Sample Name`)]
  thesamples<-unique(gsub(".+_","",gsub("_trans.+","",qctab$`Sample Name`[which(as.numeric(gsub("%","",qctab$`% Aligned`))>=70)])))
  ## Select only the counts for the samples with enough quality and put them as numbers and not characters ##
  rawcounts<-read.csv2(paste("GREIN/Counts/",intersect(counts[grep(gsub(".txt","",rownames(workingfiles)[10]),counts)],counts[grep("Raw",counts)]),sep=""),stringsAsFactors = F,sep=",")
  therawcounts<-c()
  for(a in 1:length(thesamples)){
    if(length(which(colnames(rawcounts)==thesamples[a]))==0){print(paste("Be careful with:",a))}
    therawcounts<-cbind(therawcounts,as.numeric(rawcounts[,thesamples[a]]))
  }
  rownames(therawcounts)<-rawcounts$X ; colnames(therawcounts)<-thesamples
  ## put the metadata in the proper order ##
  metadata<-fread(paste("GREIN/Metadata/",metadatas[grep(rownames(workingfiles)[10],metadatas)],sep=""),stringsAsFactors = F)
  metadata<-metadata[,c(2,10,9,12,15,3,13)]
  colnames(metadata)<-c("sample","organism","tissue","sex","age","category","characteristic1")
  metadata$category<-gsub(" [reanalysis of GSM1580881]","",metadata$category)
  metadata<-cbind(metadata,metadata$category) ; colnames(metadata)[length(metadata[1,])]<-"detailed_category"
  metadata$category[grep("^C_",metadata$category)]<-"control"
  metadata$category[grep("^P_",metadata$category)]<-"case"
  metadata$sex<-gsub("gender: female","female",metadata$sex) ; metadata$sex<-gsub("gender: male","male",metadata$sex)
  metadata$age<-gsub("age at death: ","",metadata$age)
  metadata<-metadata[,c(1,6,4,5,3,2,7,8)]
  cuales<-c() ; for(a in 1:length(thesamples)){cuales<-c(cuales,which(metadata$sample==thesamples[a]))}
  metadata<-metadata[cuales,]
  ## create the summarized experiment ##
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments/PD_brain_GSE68719.rds")
  ## Only with the cases and controls in the category column ##
  cuales<-c(which(metadata$category=="case"),which(metadata$category=="control"))
  metadata<-metadata[cuales,] ; therawcounts<-therawcounts[,cuales]
  dim(metadata) ; dim(therawcounts)
  sumexp<-SummarizedExperiment(assays=list(counts=therawcounts),colData=metadata)
  saveRDS(sumexp,file="GREIN/SummarizedExperiments_CaseControl/PD_brain_GSE68719.rds")
}


## Remove lowly expressed genes and prepare the data for conducting differential expression analysis ##
if(args[1]=="organize_data_as_needed"){
  if("processed_summarized_experiments"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/processed_summarized_experiments")}
  # if("Normalized_counts"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/Normalized_counts")}
  ficheros<-list.files("GREIN/SummarizedExperiments_CaseControl/")
  for(zeta in 1:length(ficheros)){
    # zeta<-1
    sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",ficheros[zeta],sep=""))
    ## Create the DGEList object ##
    ## @@ @@ @@ @@ @ @@ @@ @@ @@ ##
    dgeobj<-DGEList(counts = assays(sumexp)$counts, genes = as.data.frame(mcols(sumexp)), group = sumexp$category)
    assays(sumexp)$logCPM<-cpm(dgeobj, log=TRUE, prior.count=0.5)
    # Remove lowly expressed genes
    ezabatu<-rowSums(assays(sumexp)$logCPM <= 1) >= min(table(sumexp$category))
    fsumexp<-sumexp[!ezabatu, ]
    fdgeobj<-dgeobj[!ezabatu, ]
    # Calculate normalizing factors
    nfdgeobj<-calcNormFactors(fdgeobj)
    saveRDS(nfdgeobj, paste("GREIN/processed_summarized_experiments/dgeobj_",ficheros[zeta],sep=""))
    
    ## Create the DESeqDataSet object ##
    ## @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ ##
    ## here we are including the contrast matrix so it will be different for each dataset ##
    ## as we want to take into consideration different variables that are not present in all the studies ##
    if(zeta == 1){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + sex, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 2){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + tissue, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 3){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + treatment, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 4){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 5){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + sex, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 6){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + sex, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 7){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + tissue, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 8){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + sex + tissue, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 9){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category  + sex, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    if(zeta == 10){
      deseqobj<-DESeqDataSetFromMatrix(countData=round(nfdgeobj$counts), colData=colData(sumexp), design = ~ category, tidy = F) 
      deseqobj<-DESeq(deseqobj)
      saveRDS(deseqobj, paste("GREIN/processed_summarized_experiments/deseqobj_",ficheros[zeta],sep=""))
    }
    print(paste(round((zeta/length(ficheros))*100,2),"%",sep=""))
  }
}

## Conduct differential expression analyses ##
if(args[1]=="differential_expression_analyses"){
  workingfiles<-fread("GREIN/Summary_by_study_working_datasets.txt",stringsAsFactors = F,sep="\t")
  workingfiles<-cbind(workingfiles,0,0,0,0,0,0,0,0,0,0)
  colnames(workingfiles)[c(1,6:15)]<-c("dataset","case_control","sex","tissue","severity","age","adjusting_for","limma_up","limma_down","deseq2_up","deseq2_down")
  if("Differential_expression_profiles"%in%list.files("GREIN/")==FALSE){dir.create("GREIN/Differential_expression_profiles")}
  ficheros<-list.files("GREIN/processed_summarized_experiments/")
  nfdgeobjfiles<-ficheros[grep("dgeobj_",ficheros)] ;  deseqobjfiles<-ficheros[grep("deseqobj_",ficheros)]
  
  for(zeta in 1:length(deseqobjfiles)){
    ## differential expression analysis with deseq2 ##
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    # zeta<-1
    deseqobj<-readRDS(paste("GREIN/processed_summarized_experiments/",deseqobjfiles[zeta],sep=""))
    results<-results(object = deseqobj, contrast = c("category","case", "control"), alpha = 0.05, pAdjustMethod = "BH") 
    summary(results)
    deseqres<-as.data.frame(results)
    write.table(deseqres,paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",deseqobjfiles[zeta]),sep=""),quote=F,sep="\t")
    workingfiles$deseq2_up[zeta]<-length(intersect(which(deseqres$padj<=0.05),which(deseqres$log2FoldChange>0)))
    workingfiles$deseq2_down[zeta]<-length(intersect(which(deseqres$padj<=0.05),which(deseqres$log2FoldChange<0)))
    ## differential expression analysis with limma ##
    ## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ ##
    if(zeta==1){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + sex, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"yes" ; workingfiles$severity[zeta]<-"yes" ; workingfiles$age[zeta]<-"yes" ; workingfiles$adjusting_for[zeta]="sex" ; workingfiles$tissue[zeta]<-1
    }
    if(zeta==2){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + tissue, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"no" ; workingfiles$severity[zeta]<-"no" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="tissue" ; workingfiles$tissue[zeta]<-3
    }
    if(zeta==3){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + treatment, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"no" ; workingfiles$severity[zeta]<-"no" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="treatment" ; workingfiles$tissue[zeta]<-1
    }
    if(zeta==4){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"no" ; workingfiles$severity[zeta]<-"no" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="no" ; workingfiles$tissue[zeta]<-1
    }
    if(zeta==5){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + sex, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"yes" ; workingfiles$severity[zeta]<-"yes" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="sex" ; workingfiles$tissue[zeta]<-1
    }
    if(zeta==6){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + sex, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"yes" ; workingfiles$severity[zeta]<-"yes" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="sex" ; workingfiles$tissue[zeta]<-1
    }
    if(zeta==7){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + tissue, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"no" ; workingfiles$severity[zeta]<-"no" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="tissue" ; workingfiles$tissue[zeta]<-5
    }
    if(zeta==8){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + tissue + sex, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"yes" ; workingfiles$severity[zeta]<-"no" ; workingfiles$age[zeta]<-"no" ; workingfiles$adjusting_for[zeta]="tissue + sex" ; workingfiles$tissue[zeta]<-2
    }
    if(zeta==9){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category + sex, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"yes" ; workingfiles$severity[zeta]<-"yes" ; workingfiles$age[zeta]<-"yes" ; workingfiles$adjusting_for[zeta]="sex" ; workingfiles$tissue[zeta]<-1
    }
    if(zeta==10){
      sumexp<-readRDS(paste("GREIN/SummarizedExperiments_CaseControl/",gsub("dgeobj_","",nfdgeobjfiles[zeta]),sep=""))
      nfdgeobj<-readRDS(paste("GREIN/processed_summarized_experiments/",nfdgeobjfiles[zeta],sep=""))
      sumexp$category<-tryCatch({
        # Ensure 'category' is a factor before releveling
        if (!is.factor(sumexp$category)) {
          sumexp$category<-factor(sumexp$category)
        }
        relevel(sumexp$category, ref = "control")
      }, error = function(e) {
        cat("Error at iteration", i, ":", e$message, "\n")
        sumexp$category  # Return the original variable if error
      })
      ## Create the design matrix ##
      designmat<-model.matrix(~sumexp$category, colData(sumexp))
      # Mean-variance relationship
      elist<-voom(nfdgeobj, designmat, plot = TRUE)
      fit<-lmFit(elist, designmat)
      efit<-eBayes(fit)
      results<-decideTests(efit, p.value = 0.05)
      summary(results)
      limmares<-topTable(efit, coef = 2, n = Inf)
      write.table(limmares,file=paste("GREIN/Differential_expression_profiles/",gsub(".rds",".csv",nfdgeobjfiles[zeta]),sep=""),quote=F,sep="\t")
      
      workingfiles$case_control[zeta]<-paste(length(which(colData(sumexp)$category=="case")),"/",length(which(colData(sumexp)$category=="control")),sep="")
      workingfiles$limma_up[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC>0)))
      workingfiles$limma_down[zeta]<-length(intersect(which(limmares$adj.P.Val<=0.05),which(limmares$logFC<0)))
      workingfiles$sex[zeta]<-"males" ; workingfiles$severity[zeta]<-"no" ; workingfiles$age[zeta]<-"yes" ; workingfiles$adjusting_for[zeta]="no" ; workingfiles$tissue[zeta]<-1
    }
  }
  write.table(workingfiles,"GREIN/Summary_number_sDEGs.txt",quote=F,sep="\t",row.names=F)
}








