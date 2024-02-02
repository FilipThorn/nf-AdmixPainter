#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

max<-as.numeric(args[1])
chr<-as.character(args[2])

if(max > 40000000){
  fixed<-9
}else{
  if(max > 20000000){
    fixed<-7
  }else{
    fixed<-5
  }
}

window=as.integer((max/fixed))

i=1
j=window

tot<-data.frame()
while( (j+10)  < max ){
 
    string<-paste0(chr, ":", format(i,scientific = F), "-", format(j,scientific = F) )
    string
    tot<-rbind(tot, string)

    tot

    i = i+window
    j = j+window
}

end<-paste0(chr,":",format(i,scientific = F),"-", format(max,scientific = F) )
tot<-rbind(tot,end)

colnames(tot)<-c("region")

path<-paste0(chr,"_windows.txt")

write.table(tot, file=path, row.names=F, quote = F, col.names = F)
