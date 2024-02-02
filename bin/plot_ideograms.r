#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(paste0("Loadning packages and dependensies"))
library("tidyverse")
library("ggnewscale")
library("cowplot")
print(paste0("Finnished packages and dependensies"))

#list of paths to the csv files
list_tmp=args[1]
list<-read.table(list_tmp)

chromosome<-read.table(args[2], header=T)
colnames(chromosome)<-c("chrom","max")
order<-chromosome[order(chromosome$max,  decreasing = F),]$chrom

chrSize<-data.frame()
for(j in 1:length(chromosome$chrom)){
  chr<-chromosome[chromosome$chrom == chromosome$chrom[j],]
  if(chr$chrom == "chrZ" || chr$chrom == "chrW" || chr$chrom == "chrY"){
    chr$size<-"Zsex"
    chrSize<-rbind(chrSize,chr)
    }else{
    if(chr$max > 40000000){
      chr$size<-"macro"
      chrSize<-rbind(chrSize,chr)
    }else{
      if(chr$max > 20000000){
        chr$size<-"medium"
        chrSize<-rbind(chrSize,chr)
      }else{
        if(chr$max < 6000000){
        chr$size<-"ultramicro"
        chrSize<-rbind(chrSize,chr) 
        }else{
        chr$size<-"micro"
        chrSize<-rbind(chrSize,chr)
        }
      }
    }
  }
}

ordering<-arrange(chrSize, interaction(desc(max),size))

df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=TRUE))
colnames(df)<-"csv"

for(x in 1:length(df$csv)){
  df$chrom[x]<-unlist(str_split(df$csv[x], "_"))[1]
}

df_ord<-arrange(merge(ordering, df, by="chrom"), interaction(desc(max),size))

meanAdmix<-data.frame()
for( i in 1:length(df_ord$csv)){
  listtmp<-data.frame()
  list2<-read.csv(df_ord$csv[i])
  for( u in 1:length(unique(list2$region))){
    sub<-unique(list2[list2$position == sort(unique(list2$position))[u], ])
    sub$index<-u
    #if(unique(sub$divergent) == "no"){
    #  sub$meanAdmix<-NA
    #}
    if(u == 1){
      sub1<-sub
      sub2<-sub 
      sub3<-sub
      sub4<-sub
      sub5<-sub
      sub6<-sub
      range<-unlist(str_split(sub$region, ":"))[2]
      sub3$range<-as.numeric(unlist(str_split(range, "-"))[2]) - (as.numeric(unlist(str_split(range, "-"))[2]) * 0.75)
      sub4$range<-as.numeric(unlist(str_split(range, "-"))[2]) - (as.numeric(unlist(str_split(range, "-"))[2]) * 0.75)
      sub3$x<-i
      sub4$x<-i+0.5
      sub5$range<-as.numeric(unlist(str_split(range, "-"))[2])
      sub6$range<-as.numeric(unlist(str_split(range, "-"))[2])
      sub5$x<-i+0.5
      sub6$x<-i
      #rounded endcaps
      sub1$range<-1
      sub2$range<-1
      sub1$x<-i+0.15
      sub2$x<-i+0.35
      sub_comb<-rbind(sub1,sub2,sub4,sub5,sub6, sub3)
      sub_comb$ILSpoint<-i+0.25
      listtmp<-rbind(listtmp,sub_comb)
    } else{
        if(u == max(length(unique(list2$region)))){
        sub1<-sub
        sub2<-sub
        sub3<-sub
        sub4<-sub
        sub5<-sub
        sub6<-sub
        range<-unlist(str_split(sub$region, ":"))[2]
        sub1$range<-as.numeric(unlist(str_split(range, "-"))[1])
        sub2$range<-as.numeric(unlist(str_split(range, "-"))[1])
        sub1$x<-i
        sub2$x<-i+0.5
        sub3$range<-as.numeric(unlist(str_split(range, "-"))[2]) - ((as.numeric(unlist(str_split(range, "-"))[2]) - as.numeric(unlist(str_split(range, "-"))[1]) )* 0.25)
        sub4$range<-as.numeric(unlist(str_split(range, "-"))[2]) - ((as.numeric(unlist(str_split(range, "-"))[2]) - as.numeric(unlist(str_split(range, "-"))[1]) )* 0.25)
        sub3$x<-i
        sub4$x<-i+0.5
        #rounded endcaps
        sub5$range<-as.numeric(unlist(str_split(range, "-"))[2])
        sub6$range<-as.numeric(unlist(str_split(range, "-"))[2])
        sub5$x<-i+0.35
        sub6$x<-i+0.15
        sub_comb<-rbind(sub6,sub5,sub4,sub2,sub1,sub3)
        sub_comb$ILSpoint<-i+0.25
        listtmp<-rbind(listtmp,sub_comb)
        }else{
        sub1<-sub
        sub2<-sub
        sub3<-sub
        sub4<-sub
        range<-unlist(str_split(sub$region, ":"))[2]
        sub1$range<-as.numeric(unlist(str_split(range, "-"))[1])
        sub2$range<-as.numeric(unlist(str_split(range, "-"))[1])
        sub1$x<-i
        sub2$x<-i+0.5
        sub3$range<-as.numeric(unlist(str_split(range, "-"))[2])
        sub4$range<-as.numeric(unlist(str_split(range, "-"))[2])
        sub3$x<-i+0.5
        sub4$x<-i
        sub_comb<-rbind(sub1,sub2,sub3,sub4)
        sub_comb$ILSpoint<-i+0.25
        listtmp<-rbind(listtmp,sub_comb)
        }
    }
  }
  meanAdmix<-rbind(meanAdmix,listtmp)
}

meanAdmix<-merge(meanAdmix, ordering, by="chrom")
meanAdmix$meanAdmix<-as.numeric(meanAdmix$meanAdmix)
meanAdmix$medianAdmix<-as.numeric(meanAdmix$medianAdmix)

##make shapes for chromosomes 
#take out regions higher than the mean admixture
#mark rest as NA
#Regions should be protected from ILS 
#Draw polygonal shapes for the regions as chromosomes 

meanAll<-arrange(meanAdmix,interaction(desc(max),size))

all<-ggplot(meanAll, aes(x = x, y = range)) +
  geom_polygon(data=meanAll, aes(fill = medianAdmix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
  theme_void() +
  scale_fill_gradient(high="red", low="white", limits=c(0,1))+
  guides(fill="none")
  #geom_point(data=meanAll, aes(x = ILSpoint, y = position, shape=divergent), size =0.1, col="blue")+
  #scale_shape_manual(values=c(8, NA)) +
  #guides(shape="none") 



 all_ILS<- ggplot(meanAll, aes(x = x, y = range)) +
  geom_polygon(data=meanAll, aes(fill = medianAdmix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
  theme_bw()+theme_void()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_fill_gradient(high="red", low="white", limits=c(0,1))+
  guides(fill="none") + 
  geom_point(data=meanAll, aes(x = ILSpoint, y = position, shape=divergent), size =0.5, col="blue")+
  scale_shape_manual(values=c(no=4, yes=NA)) +
  guides(shape="none") + scale_x_continuous(breaks = 1.25:(length(unique(meanAll$chrom)) + 0.25),
                                            labels = unique(meanAll$chrom))

plot<-plot_grid(all, all_ILS, ncol=1, nrow=2)

ggsave("Ideograms.pdf",
       plot,
       device = "pdf",
       scale = 1,
       width = 6,
       height = 6,
       units = "in")

#separate into macro, medium ,  micro , ultramicro
mac<-meanAdmix[meanAdmix$size %in% c("macro"),]
length(unique(mac$chrom))

med<-meanAdmix[meanAdmix$size %in%  c("medium"),]

mic<-meanAdmix[meanAdmix$size %in%  c("micro"),]

umic<-meanAdmix[meanAdmix$size %in%  c("ultramicro"),]

med$x<-med$x-length(unique(mac$chrom))
med$ILSpoint<-med$ILSpoint-length(unique(mac$chrom))
me<-arrange(med,interaction(desc(max),size))

umic$x<-umic$x-length(unique(mac$chrom))-length(unique(med$chrom))-length(unique(mic$chrom))
umic$ILSpoint<-umic$ILSpoint-length(unique(mac$chrom))-length(unique(med$chrom))-length(unique(mic$chrom))
umi<-arrange(umic,interaction(desc(max),size))


mic$x<-mic$x-length(unique(mac$chrom))-length(unique(med$chrom))
mic$ILSpoint<-mic$ILSpoint-length(unique(mac$chrom))-length(unique(med$chrom))
mi<-arrange(mic,interaction(desc(max),size))

sex<-meanAdmix[meanAdmix$size %in% c( "Zsex"),]
sex$x<-sex$x-length(unique(med$chrom))-length(unique(mic$chrom))-length(unique(umic$chrom))
sex$ILSpoint<-sex$ILSpoint-length(unique(med$chrom))-length(unique(mic$chrom))-length(unique(umic$chrom))

mas<-arrange(rbind(mac, sex),interaction(desc(max),size))

##shape guide 
shape_guide<-c("yes"=NA,
                "no"=4)


#macro and sex chromosomes
pmas<-ggplot(mas, aes(x = x, y = range)) +
  geom_polygon(data=mas, aes(fill = medianAdmix,  group = interaction(chrom, index)), 
               col="black", lwd=0.1) +
  theme_bw() +
  scale_fill_gradient(high="red", low="white",
                      name="Admix proportion",
                      breaks=c(0.25,0.5, 0.75),
                      limits=c(0,1)) +
  geom_point(data=mas, aes(x = ILSpoint, y = position, shape=divergent), 
             size =0.5, col="blue")+
  scale_shape_manual(values=shape_guide) +
  guides(shape="none")  + scale_x_continuous(breaks = 1.25:(length(unique(mas$chrom)) + 0.25),
                                             labels = unique(mas$chrom), limits = c(0.5,
                                            max(c(length(unique(mas$chrom)),
                                                  length(unique(me$chrom)),
                                                  length(unique(mi$chrom)),
                                                  length(unique(umi$chrom))))+0.5))+
  xlab("") + ylab("Scaffold size") + ggtitle("Macro and Z chromosomes") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.x = element_blank(),
        title = element_text(size=8))+
  guides(fill="none") 

#intermediate chromosomes
pme<-ggplot(me, aes(x = x, y = range)) +
  geom_polygon(data=me, aes(fill = medianAdmix,  group = interaction(chrom, index)), 
               col="black", lwd=0.1) +
  theme_bw() +
  scale_fill_gradient(high="red", low="white", 
                      name="Admix proportion", 
                      breaks=c(0.25,0.5, 0.75), 
                      limits=c(0,1)) +
  geom_point(data=me, aes(x = ILSpoint, y = position, shape=divergent), 
             size =0.5, col="blue")+
  scale_shape_manual(values=shape_guide) +
  guides(shape="none") +
  scale_x_continuous(breaks = 1.25:(length(unique(me$chrom)) + 0.25),
                     labels = unique(me$chrom), limits = c(0.5,
                     max(c(length(unique(mas$chrom)),
                           length(unique(me$chrom)),
                           length(unique(mi$chrom)),
                           length(unique(umi$chrom))))+0.5))+
  xlab("") + ylab("")+ ggtitle("Intermediate chromosomes") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.x = element_blank(),
        title = element_text(size=8)) +
  guides(fill="none")

#micro chromosomes
pmi<-ggplot(mi, aes(x = x, y = range)) +
  geom_polygon(data=mi, aes(fill = medianAdmix,  group = interaction(chrom, index)), 
               col="black", lwd=0.1) +
  theme_bw() +
  scale_fill_gradient(high="red", low="white", 
                      name="Admix proportion", 
                      breaks=c(0.25,0.5, 0.75), 
                      limits=c(0,1)) +
  geom_point(data=mi, aes(x = ILSpoint, y = position, shape=divergent), 
             size =0.5, col="blue")+
  scale_shape_manual(values=shape_guide) +
  guides(shape="none") +
  scale_x_continuous(breaks = 1.25:(length(unique(mi$chrom)) + 0.25),
                     labels = unique(mi$chrom), limits = c(0.5,
                     max(c(length(unique(mas$chrom)),
                           length(unique(me$chrom)),
                           length(unique(mi$chrom)),
                           length(unique(umi$chrom))))+0.5))+
  xlab("") + ylab("Scaffold size")+ ggtitle("Micro chromosomes") +   
  theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
       axis.title.x = element_blank(),
       plot.background = element_blank(),
       panel.grid = element_blank(),
       panel.border = element_blank(),
       axis.line.y = element_line(),
       axis.ticks.x = element_blank(),
       title = element_text(size=8)) +
  guides(fill="none")

#ultra micro
pumi<-ggplot(umi, aes(x = x, y = range)) +
  geom_polygon(data=umi, aes(fill = medianAdmix,  group = interaction(chrom, index)), 
               col="black", lwd=0.1) +
  theme_bw() +
  scale_fill_gradient(high="red", low="white",
                      name="Admix\nproportion",
                      breaks=c(0.25,0.5, 0.75),
                      limits=c(0,1)) +
  geom_point(data=umi, aes(x = ILSpoint, y = position, shape=divergent),
             size =0.5, col="blue")+
  scale_shape_manual(values=shape_guide) +
  guides(shape="none")  + 
  scale_x_continuous(breaks = 1.25:(length(unique(umi$chrom)) + 0.25),
                     labels = unique(umi$chrom), limits = c(0.5,
                    max(c(length(unique(mas$chrom)),
                          length(unique(me$chrom)),
                          length(unique(mi$chrom)),
                          length(unique(umi$chrom))))+0.5))+
  xlab("") + ylab("") + ggtitle("Ultra micro chromosomes") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
  axis.title.x = element_blank(),
  plot.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  axis.line.y = element_line(),
  axis.ticks.x = element_blank(),
  legend.position = c(0.8,0.3),
  legend.text = element_text(size=8),
  title = element_text(size=8),
  legend.key.size = unit(0.5, 'cm'))


plot<-plot_grid(pmas, pme, pmi, pumi, ncol=2, nrow=2)


ggsave("IdeogramSize.pdf",
       plot,
       device = "pdf",
       scale = 1,
       width = 6,
       height = 8,
       units = "in")

#wombo combo
allLabs<-all +
        theme_bw()+theme_void()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8),
              legend.title = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank()) +
        scale_x_continuous(breaks = 1.25:(length(unique(meanAll$chrom)) + 0.25),
                           labels = unique(meanAll$chrom))

WC<-plot_grid(allLabs, plot, ncol=1, nrow=2, rel_heights = c(1, 1))


ggsave("IdeogramCombine.pdf",
       WC,
       device = "pdf",
       scale = 1,
       width = 6,
       height = 8,
       units = "in")


write.table(meanAdmix, file="./plotting.csv", row.names=F, quote = F, col.names = T)


###################################################################################
#Chromosomal patterns
print("starting chromosmal patterns")

all_blue<-all + 
          new_scale_fill() + geom_polygon(data=meanAll, aes(fill = divergent,  group = interaction(chrom, index), alpha=divergent))+
          scale_fill_manual(values = c(yes="white", no="blue")) +
          scale_alpha_manual(values = c(yes=0, no=1))+ 
          guides(fill="none")+
          guides(alpha="none")


micros<-rbind(mi, umi)
fmicros<-micros[micros$divergent %in% "yes",]
allMI<-ggplot()+
          geom_boxplot(data=fmicros, aes(admix, factor(index)),width=0.5)+ 
          coord_flip() +
          theme_bw()+
          ggtitle(paste0("Micro chromosomes")) +
          ylab("Chromosome window") +
          xlab("Admix proportion")+ 
          theme(title = element_text(size=6))


fhme<-me[me$divergent %in% "yes",]
allME<-ggplot()+
        geom_boxplot(data=fhme, aes(admix, factor(index)),width=0.5)+ 
        coord_flip() +
        theme_bw()+
        ggtitle(paste0("Intermediate chromosomes")) +
        ylab("Chromosome window") +
        xlab("Admix proportion")+ 
        theme(title = element_text(size=6))

#remove sex chromosome
hmac<-mas[!mas$chrom %in% "chrZ",]
fhmac<-hmac[mac$divergent %in% "yes",]
allMA<-ggplot()+
          geom_boxplot(data=fhmac, aes(admix,factor(index)),width=0.5)+ 
          coord_flip() +
          theme_bw()+
          ggtitle(paste0("Macro chromosomes")) +
          ylab("Chromosome window") +
          xlab("Admix proportion")+ 
          theme(title = element_text(size=6))



allChrPatterns<-plot_grid(allMA,allME, allMI,all_blue, ncol=2, nrow=2, rel_heights = c(1,1,1,1 ))

all_title<-paste0("Chromosomal_patterns_all.pdf")

ggsave(all_title,
       allChrPatterns,
       device = "pdf",
       scale = 1,
       width = 7,
       height = 6,
       units = "in")


#####################################################################################
#seperate hybrids
#####################################################################################
print("starting seperate hybrid plots loop") 


for(j in 1:length(unique(meanAll$ID))){

    hyb<-meanAll[meanAll$ID %in% unique(meanAll$ID)[j],]

    hybrid<- unique(meanAll$ID)[j]
    
    all<-ggplot(hyb, aes(x = x, y = range)) +
        geom_polygon(data=hyb, aes(fill = admix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
        theme_void() +
        scale_fill_gradient(high="red", low="white", limits=c(0,1))+
        guides(fill="none")

    hmas<-mas[mas$ID %in% unique(meanAll$ID)[j],]

    ###macro and sex chromosomes###
    pmas<-ggplot(hmas, aes(x = x, y = range)) +
        geom_polygon(data=hmas, aes(fill = admix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
        theme_bw() +
        scale_fill_gradient(high="red", low="white",
                      name="Admix proportion",
                      breaks=c(0.25,0.5, 0.75),
                      limits=c(0,1)) +
        geom_point(data=hmas, aes(x = ILSpoint, y = position, shape=divergent), size =0.5, col="blue")+
        scale_shape_manual(values=shape_guide) +
        guides(shape="none") + 
        scale_x_continuous(breaks = 1.25:(length(unique(mas$chrom)) + 0.25),
                                             labels = unique(mas$chrom), limits = c(0.5,
                                            max(c(length(unique(mas$chrom)),
                                                  length(unique(me$chrom)),
                                                  length(unique(mi$chrom)),
                                                  length(unique(umi$chrom))))+0.5))+
        xlab("") + 
        ylab("Scaffold size") + 
        ggtitle("Macro and Z chromosomes") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
                                        axis.title.x = element_blank(),
                                        plot.background = element_blank(),
                                        panel.grid = element_blank(),
                                        panel.border = element_blank(),
                                        axis.line.y = element_line(),
                                        axis.ticks.x = element_blank(),
                                        title = element_text(size=8))+
                                        guides(fill="none")

    ###intermediate chromosomes###

    hme<-me[me$ID %in% unique(meanAll$ID)[j],]
    
    pme<-ggplot(hme, aes(x = x, y = range)) +
        geom_polygon(data=hme, aes(fill = admix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
        theme_bw() +
        scale_fill_gradient(high="red", low="white",
                      name="Admix proportion",
                      breaks=c(0.25,0.5, 0.75),
                      limits=c(0,1)) +
        geom_point(data=hme, aes(x = ILSpoint, y = position, shape=divergent), size =0.5, col="blue") +
        scale_shape_manual(values=shape_guide) +
        guides(shape="none") +
        scale_x_continuous(breaks = 1.25:(length(unique(me$chrom)) + 0.25),
                           labels = unique(me$chrom), limits = c(0.5,
                           max(c(length(unique(mas$chrom)),
                           length(unique(me$chrom)),
                           length(unique(mi$chrom)),
                           length(unique(umi$chrom))))+0.5))+
        xlab("") + 
        ylab("") + 
        ggtitle("Intermediate chromosomes") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
                axis.title.x = element_blank(),
                plot.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(),
                axis.ticks.x = element_blank(),
                title = element_text(size=8)) +
        guides(fill="none")

    ###micro chromosomes###
    
    hmi<-mi[mi$ID %in% unique(meanAll$ID)[j],]

    pmi<-ggplot(hmi, aes(x = x, y = range)) +
        geom_polygon(data=hmi, aes(fill = admix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
        theme_bw() +
        scale_fill_gradient(high="red", low="white",
                              name="Admix proportion",
                              breaks=c(0.25,0.5, 0.75),
                              limits=c(0,1)) +
        geom_point(data=hmi, aes(x = ILSpoint, y = position, shape=divergent), size =0.5, col="blue") +
        scale_shape_manual(values=shape_guide) +
        guides(shape="none") +
        scale_x_continuous(breaks = 1.25:(length(unique(mi$chrom)) + 0.25),
                         labels = unique(mi$chrom), limits = c(0.5,
                         max(c(length(unique(mas$chrom)),
                           length(unique(me$chrom)),
                           length(unique(mi$chrom)),
                           length(unique(umi$chrom))))+0.5)) +
        xlab("") + 
        ylab("Scaffold size") + 
        ggtitle("Micro chromosomes") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
                axis.title.x = element_blank(),
                plot.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(),
                axis.ticks.x = element_blank(),
                title = element_text(size=8)) +
        guides(fill="none")

    ###ultra micro###
    humi<-umi[umi$ID %in% unique(meanAll$ID)[j],]
    
    pumi<-ggplot(humi, aes(x = x, y = range)) +
            geom_polygon(data=humi, aes(fill = admix,  group = interaction(chrom, index)), col="black", lwd=0.1) +
            theme_bw() +
            scale_fill_gradient(high="red", low="white",
                                    name="Admix\nproportion",
                                    breaks=c(0.25,0.5, 0.75),
                                    limits=c(0,1)) +
            geom_point(data=humi, aes(x = ILSpoint, y = position, shape=divergent), size =0.5, col="blue")+
            scale_shape_manual(values=shape_guide) +
            guides(shape="none")  +
            scale_x_continuous(breaks = 1.25:(length(unique(umi$chrom)) + 0.25),
                                labels = unique(umi$chrom), limits = c(0.5,
                                max(c(length(unique(mas$chrom)),
                                length(unique(me$chrom)),
                                length(unique(mi$chrom)),
                                length(unique(umi$chrom))))+0.5))+
            xlab("") + 
            ylab("") + 
            ggtitle("Ultra micro chromosomes") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8, colour ="black"),
                    axis.title.x = element_blank(),
                    plot.background = element_blank(),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    axis.line.y = element_line(),
                    axis.ticks.x = element_blank(),
                    legend.position = c(0.8,0.3),
                    legend.text = element_text(size=8),
                    title = element_text(size=8),
                    legend.key.size = unit(0.5, 'cm'))


    ###plot grid###
    plot<-plot_grid(pmas, pme, pmi, pumi, ncol=2, nrow=2)
    ###wombo combo###
    
    allLabs<-all +
                theme_bw()+theme_void()+
                theme(axis.text.x = element_text(angle = 45, vjust = 1.5, hjust=1, size=8),
                    legend.title = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank()) +
                    scale_x_continuous(breaks = 1.25:(length(unique(meanAll$chrom)) + 0.25),
                                        labels = unique(meanAll$chrom))

    
    ###plot full grid###
    WC<-plot_grid(allLabs, plot, ncol=1, nrow=2, rel_heights = c(1, 1))

    title<-paste0(unique(meanAll$ID)[j],"_IdeogramCombine.pdf")

    ggsave(title,
            WC,
            device = "pdf",
            scale = 1,
            width = 6,
            height = 8,
            units = "in")

    micros<-rbind(hmi, humi)
    fmicros<-micros[micros$divergent %in% "yes",]
    mip<-ggplot()+
        geom_boxplot(data=fmicros, aes(admix, factor(index)),width=0.5)+ 
        coord_flip() +
        theme_bw()+
        ggtitle(paste0(hybrid, "; micro chromosomes")) +
        ylab("Chromosome window") +
        xlab("Admix proportion")+ 
        theme(title = element_text(size=6))


    fhme<-hme[hme$divergent %in% "yes",]
    mep<-ggplot()+
        geom_boxplot(data=fhme, aes(admix, factor(index)),width=0.5)+ 
        coord_flip() +
        theme_bw() +
        ggtitle(paste0(hybrid, "; intermediate chromosomes")) +
        ylab("Chromosome window") +
        xlab("Admix proportion")+ 
        theme(title = element_text(size=6))

    #remove sex chromosome
    hmac<-hmas[!hmas$chrom %in% "chrZ",]
    fhmac<-hmac[hmac$divergent %in% "yes",]
    map<-ggplot()+
        geom_boxplot(data=fhmac, aes(admix,factor(index)),width=0.5)+ 
        coord_flip() +
        theme_bw()+
        ggtitle(paste0(hybrid, "; macro chromosomes")) +
        ylab("Chromosome window") +
        xlab("Admix proportion") + 
        theme(title = element_text(size=6))

    hyb_blue<-all + 
              new_scale_fill() + geom_polygon(data=hyb, aes(fill = divergent,  group = interaction(chrom, index), alpha=divergent))+
              scale_fill_manual(values = c(yes="white", no="blue")) +
              scale_alpha_manual(values = c(yes=0, no=1))+ 
              guides(fill=FALSE)+
              guides(alpha=FALSE)



    chrpat<-plot_grid(map,mep, mip,hyb_blue, ncol=2, nrow=2, rel_heights = c(1,1,1,1 ))

    hyb_title<-paste0(unique(meanAll$ID)[j],"_chromosomal_patterns.pdf")

    ggsave(hyb_title,
            chrpat,
            device = "pdf",
            scale = 1,
            width = 7,
            height = 6,
            units = "in")
}
