logs<-read.table("./logs.txt")

max<-subset(logs, logs$V1 == max(logs[1]))[1,]

write.table(max, file="./best.val", row.names=F, quote = F, col.names = F)
