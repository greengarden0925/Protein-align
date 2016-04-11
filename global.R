source("http://bioconductor.org/biocLite.R")
mpackages=c("shiny","seqinr","Biostrings","tcltk","gtools","plotrix","googleVis","shinydashboard")
CheckandInstallPackages=function(mpackages){
  for (p in mpackages){
    lapply(p, require, character.only=T)
    if (!lapply(p,exists)[[1]]){ 
      if (!lapply(p, require, character.only=T)[[1]]){ 
        cat(p,"此package不存在下載中","\n") 
        lapply(p,biocLite)
        lapply(p, require, character.only=T)
      }
    }
  }
}
CheckandInstallPackages(mpackages)

convert_string_upper=function(seq){
  seq_string=c2s(seq)
  seq_string_toupper=toupper(seq_string)
  return(seq_string_toupper)
}

alignPossibleMaxScore=function(seq1,seq2){
  seq1no=nchar(seq1) 
  seq2no=nchar(seq2) 
  if (seq1no>=seq2no) {seq=seq2} else {seq=seq1} #  
  data(BLOSUM50) 
  maxScore =pairwiseAlignment(seq,seq,substitutionMatrix=BLOSUM50,gapOpening=-2,gapExtension=-8,scoreOnly=TRUE,type="local")
  return(maxScore)
}


printAlignNo=function(seq1_list,seq2_list){
  seq1_seqs=getSequence(seq1_list) #另存序列為seq1_seqs
  seq2_seqs=getSequence(seq2_list) #另存序列為seq2_seqs
  seq1_seqs_unique=seq1_seqs[!duplicated(seq1_seqs)] 
  seq2_seqs_unique=seq2_seqs[!duplicated(seq2_seqs)]
  seq1_no=length(seq1_seqs_unique)
  seq2_no=length(seq2_seqs_unique)
  
  outcome=list()
  outcome[[1]]=seq1_no
  outcome[[2]]=seq2_no
  outcome[[3]]=seq1_no*seq2_no
  return(outcome)
}

PieChartData=function(alignoutcome, cutpoint){
  subscores=subset(alignoutcome,alignscore>=cutpoint)
  subscores=na.omit(subscores)
  seq1_tab=as.vector(subscores$seq1Name)
  seq2_tab=as.vector(subscores$seq2Name)
  seq1data_pie=data.frame(table(seq1_tab))
  seq2data_pie=data.frame(table(seq2_tab))
  Pies=list()
  Pies[[1]]= seq1data_pie
  Pies[[2]]= seq2data_pie
  return(Pies)
}
