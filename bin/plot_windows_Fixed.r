#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(paste0("Loadning packages and dependensies"))
library("tidyverse")
library("ggnewscale")
library("cowplot")
print(paste0("Finnished packages and dependensies"))

#species and populations, hybrids marked with "hybrid"
pop<-read.table(args[1], header=F)
colnames(pop)<-c("ID","pop")

#set value to percentage identity needed to bee considered a high identity window
identity=as.numeric(args[2])

#list of paths to the qopt files
list_tmp=args[3]
list<-read.table(list_tmp)

chromo=args[4]

order<-pop[order(pop$pop),]$ID

#get population names
pop1=unique(pop[!pop$pop == "hybrid",]$pop)[1]
pop2=unique(pop[!pop$pop == "hybrid",]$pop)[2]

#plotting astetics
pop_color_scale<-c("#b87800","#d9ae02", "grey")
names(pop_color_scale)<-c(pop1,pop2, "hybrid")

#first we need to name qmat columns so that the names match across windows. 
print(paste0("Starting loop to filter out nan values and standardise naming"))

melt_df<-data.frame()
for( i in 1:length(list)){
  print(paste0("print iteration: ", i))
  q_path<- paste0(list[i])
  pos_range_tmp<-unlist(str_split(list[i], ":"))[2]
  pos_range<-str_remove(pos_range_tmp, "_best.qopt")
  POS<-as.numeric(unlist(str_split(pos_range, "-")))
  pos<-(POS[1]+POS[2])/2
  qmat<-read.table(as.character(q_path), header=F)
  df<-cbind(pop,qmat)
  chr_range<-paste0(chromo,":",pos_range)
  
  if(length(df[complete.cases(df), ]$ID) > 0){
    if(sum(df[df$pop == pop1,]$V1) > sum(df[df$pop == pop1,]$V2)){
      melt1<-cbind(df[1:2], df$V1,rep(pop1, nrow(qmat)),pos,rep(chr_range, nrow(qmat)))
      colnames(melt1)[3]<-"admix"
      colnames(melt1)[4]<-"ancestry"
      colnames(melt1)[5]<-"position"
      colnames(melt1)[6]<-"region"
    
      melt2<-cbind(df[1:2], df$V2,rep(pop2, nrow(qmat)),pos,rep(chr_range, nrow(qmat)))
      colnames(melt2)[3]<-"admix"
      colnames(melt2)[4]<-"ancestry"
      colnames(melt2)[5]<-"position"
      colnames(melt2)[6]<-"region"
      
      melt_tmp<-rbind(melt1, melt2)
      
      #add point range to   plot with
      melt_tmp_POS1<-cbind(melt_tmp, POS[1])
      colnames(melt_tmp_POS1)[7]<-"range"

      melt_tmp_POS2<-cbind(melt_tmp, POS[2])
      colnames(melt_tmp_POS2)[7]<-"range"

      melt_tmp2<-rbind(melt_tmp_POS1,melt_tmp_POS2)
      melt_df<-rbind(melt_df,melt_tmp2)
        
    }else{
      if(sum(df[df$pop == pop1,]$V1) < sum(df[df$pop == pop1,]$V2)){

      melt1<-cbind(df[1:2], df$V2,rep(pop1, nrow(qmat)),pos,rep(chr_range, nrow(qmat)))
      colnames(melt1)[3]<-"admix"
      colnames(melt1)[4]<-"ancestry"
      colnames(melt1)[5]<-"position"
      colnames(melt1)[6]<-"region"
      
      melt2<-cbind(df[1:2], df$V1,rep(pop2, nrow(qmat)),pos,rep(chr_range, nrow(qmat)))
      colnames(melt2)[3]<-"admix"
      colnames(melt2)[4]<-"ancestry"
      colnames(melt2)[5]<-"position"
      colnames(melt2)[6]<-"region"
      
      melt_tmp<-rbind(melt1, melt2)
      #add point range to   plot with
      melt_tmp_POS1<-cbind(melt_tmp, POS[1])
      colnames(melt_tmp_POS1)[7]<-"range"
        
      melt_tmp_POS2<-cbind(melt_tmp, POS[2])
      colnames(melt_tmp_POS2)[7]<-"range"

      melt_tmp2<-rbind(melt_tmp_POS1,melt_tmp_POS2)

      melt_df<-rbind(melt_df,melt_tmp2)
      }else{
        print(paste0(chr_range, " unexpected error"))
        print(sum(df[df$pop == pop1,]$V1))
        print(sum(df[df$pop == pop1,]$V2))
        
        print(sum(df[df$pop == pop2,]$V1))
        print(sum(df[df$pop == pop2,]$V2))
        break
        }
    }
  }else{
        print(paste0("Window: ",chr_range, " contain nan values and were not included"))
    }
    print(paste0("Window: ",chr_range, "  included"))
}

print(paste0("Finished filtering and standardise naming loop"))

#set the ploting limits to the highest position plus size of one window
xlim<-(max(melt_df$range) + 100)

#micro, intermediate, macro chromosome 
if(max(melt_df$range) > 40000000){
    PATH=paste0("./macro/",chromo,"/")
}else{
    if(max(melt_df$range) > 20000000){
        PATH=paste0("./intermediate/",chromo,"/")
    }else{
        PATH=paste0("./micro/",chromo,"/")    
    }
}


#Admixture proportion where pop1 is always high and pop2 is always low
melt_df_pop1<-melt_df[melt_df$ancestry == pop1,]

#all
a1_title<-paste0("Admixture proportions in fixed windows (", format(length(unique(melt_df_pop1$region)), scientific=F),") on ", chromo)
a1<-ggplot()+
    geom_line(data=melt_df_pop1, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round") + 
    xlim(0,xlim) + 
    theme_bw() +
    ylab("admixture proportion")+
    scale_color_manual(values=pop_color_scale)+
    ggtitle(a1_title)+
    theme(plot.title = element_text(size = 8)) +
    guides(color=guide_legend(title=""))

ggsave(paste0(PATH,"1.AllWindows.png"),
       plot = a1,
       device = "png",
       scale = 1,
       width = 6,
       height = 4,
       units = "in",
       limitsize = F)


#Filter out only windows where the pure species have identity or (1-identity) percent in admix
pure<- melt_df_pop1[!melt_df_pop1$pop %in% c("hybrid" ),]

#plot withouth hybrids
a2_title<-paste0("Admixture proportions in fixed windows (", format(length(unique(melt_df_pop1$region)), scientific=F),") on ", chromo)
a2<-ggplot()+
    geom_line(data=pure, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+ 
    xlim(0,xlim) + 
    theme_bw() +
    ylab("admixture proportion")+
    scale_color_manual(values=pop_color_scale)+
    ggtitle(a2_title, subtitle = "No hybrids")+
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 8) ) +
    guides(color=guide_legend(title=""))

ggsave(paste0(PATH,"2.NoHybridsAllWindows.png"),
       plot = a2,
       device = "png",
       scale = 1,
       width = 6,
       height = 4,
       units = "in",
       limitsize = F)


#select only high identity windows in the pure species
pure_hi<- pure[pure$admix > identity,]
pure_low<- pure[pure$admix < (1-identity),]
pure_diver<-rbind(pure_hi, pure_low)

#Filter out windows that aren't fixed in all pure species
print(paste0("Starting loop to filter out windows that aren't fixed in all pure species"))
upos<-unique(pure_diver$position)
divergent<-data.frame()
NA_list<-data.frame()
for(k in 1:length(upos)){
  subset<-pure_diver[pure_diver$position == upos[k],]
  #check that all species are present
  if(length(subset$ID) == 2*(length(unique(pure$ID))) ){
    subset2sp1<-subset[subset$pop %in% pop1 & subset$ancestry == pop1,]
    subset3sp1<-subset[subset$pop %in% pop2 & subset$ancestry == pop1,]
    #if all species are present check that all populations are within the same range of admixture proportion
    if(length(unique(subset2sp1$admix > identity)) == 1 & length(unique(subset3sp1$admix < (1-identity))) == 1){
      divergent<-rbind(divergent, subset)
      NA_list<-rbind(NA_list, cbind(subset$region, "yes"))
    }
  }else{
    print(paste0("Filtered out: ",unique(subset$region), ", identity of window not fixed"))
    NA_list<-rbind(NA_list, cbind(subset$region, "no"))
  }
}
print(paste0("Finished loop that filtered windows in all pure species"))
colnames(NA_list)<-c("region","divergent")

###########################################
#add regions that are not fixed for any individuals to mask


logic<-unique(melt_df_pop1$region) %in% unique(NA_list$region)

missing<-unique(melt_df_pop1$region)[!logic]

print(missing)

missing<-as.data.frame(missing)
colnames(missing)<-"region"

if(length(missing$region) > 0){
missing$divergent<-"no"

print(colnames(missing))

NA_list<-rbind(NA_list,missing)

print(NA_list)
}

###########################################

#positions with high divergence
udpos<-unique(divergent$position)

if(length(udpos)> 1){
  
  #plot fixed windows no hybrids
  a3_title<-paste0("Admixture proportions in filtered fixed windows on ", chromo)
  a3<-ggplot()+
    geom_line(data=divergent, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+ 
    xlim(0,xlim) + 
    theme_bw() +
    ylab("admixture proportion")+
    scale_color_manual(values=pop_color_scale)+
    ggtitle(a3_title, subtitle = "No hybrids")+
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 8) ) +
    guides(color=guide_legend(title=""))
  
  ggsave(paste0(PATH,"3.NoHybridsFixedWindows.png"),
         plot = a3,
         device = "png",
         scale = 1,
         width = 6,
         height = 4,
         units = "in",
         limitsize = F)
  
  
  #extract divergent positions for hybrids
  hybs<- melt_df_pop1[melt_df_pop1$pop %in% c("hybrid") & melt_df_pop1$position %in% udpos,]
  hupos<-unique(hybs$position)
  
  
  a4_title<-paste0("Admixture proportions in filtered fixed windows on ", chromo)
  a4<-ggplot()+
    geom_line(data=divergent, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
    geom_line(data=hybs, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+ 
    xlim(0,xlim) + 
    theme_bw() +
    ylab("admixture proportion")+
    scale_color_manual(values=pop_color_scale)+
    ggtitle(a4_title)+
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 8) ) +
    guides(color=guide_legend(title=""))
  
  ggsave(paste0(PATH,"4.HybridsFixedWindows.png"),
         plot = a4,
         device = "png",
         scale = 1,
         width = 6,
         height = 4,
         units = "in",
         limitsize = F)
  
  #calculate the median and mean admixture proportion across all hybrids in a window
  meanAdmix<-data.frame()
  for(m in 1:length(hupos)){
    subset_hupos<-hybs[hybs$position %in% hupos[m],]
    subset_hupos$mean<-sum(subset_hupos$admix)/length(subset_hupos$admix)
    subset_hupos$median<-median(subset_hupos$admix)
    subset_hupos$chromo<-chromo  

    mean_col<-cbind(subset_hupos$ID,
                    subset_hupos$region,
                    subset_hupos$position, 
                    subset_hupos$pop,
                    subset_hupos$admix, 
                    subset_hupos$ancestry,
                    subset_hupos$mean,
                    subset_hupos$median,
                    subset_hupos$chromo)
    
    print(mean_col)

    dist_col<-mean_col

    colnames(dist_col)<-c("ID","region","position","pop","admix","ancestry","meanAdmix","medianAdmix", "chrom")
    meanAdmix<-rbind(meanAdmix, dist_col)
  }
  
  meanAdmix$meanAdmix<-as.numeric(meanAdmix$meanAdmix)
  meanAdmix$medianAdmix<-as.numeric(meanAdmix$medianAdmix)
  meanAdmix$position<-as.numeric(meanAdmix$position)
  
  #plot the mean and median admixture proportions across all hybrids
  a5_title<-paste0("Average Admixture proportions in fixed windows on ", chromo)
  a5<-ggplot()+
      geom_line(data=divergent, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=hybs, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=meanAdmix, aes(position, meanAdmix))+ 
      geom_point(data=meanAdmix, aes(position, meanAdmix), size=1)+ 
      geom_line(data=meanAdmix, aes(position, medianAdmix), col="blue")+
      geom_point(data=meanAdmix, aes(position, medianAdmix), size=1, col="blue")+   
      xlim(0,xlim) + 
      theme_bw() +
      ylab("admixture proportion")+
      scale_color_manual(values=pop_color_scale)+
      ggtitle(a5_title)+
      theme(plot.title = element_text(size = 8),
            plot.subtitle = element_text(size = 8) ) +
      guides(color=guide_legend(title=""))
  
  ggsave(paste0(PATH,"5.AverageAdmixtureProportionsFixedWindows.png"),
         plot = a5,
         device = "png",
         scale = 1,
         width = 6,
         height = 4,
         units = "in",
         limitsize = F)

    ####with ab line############

  a5_5_title<-paste0("Average Admixture proportions in fixed windows on ", chromo)
  a5_5<-ggplot()+
      geom_line(data=divergent, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=hybs, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=meanAdmix, aes(position, meanAdmix))+
      geom_point(data=meanAdmix, aes(position, meanAdmix), size=1)+
      geom_hline(yintercept=mean(meanAdmix$meanAdmix), col = "red", linetype = "dashed" )+
      xlim(0,xlim) +
      theme_bw() +
      ylab("admixture proportion")+
      scale_color_manual(values=pop_color_scale)+
      ggtitle(a5_5_title)+
      theme(plot.title = element_text(size = 8),
            plot.subtitle = element_text(size = 8) ) +
      guides(color=guide_legend(title=""))
  
  ggsave(paste0(PATH,"5.5.AblineAverageAdmixtureProportionsFixedWindows.png"),
         plot = a5_5,
         device = "png",
         scale = 1,
         width = 6,
         height = 4,
         units = "in",
         limitsize = F)


    ###########################


  print("Plotting combined plots with cowplot packages")  
  #single individual subset
  for(n in 1:length(unique(hybs$ID))){
    
    single<-hybs[hybs$ID == hybs$ID[n],]
    
    p_title<-paste0("Hybrid: ",unique(single$ID))
    p<-ggplot()+
        geom_line(data=divergent, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
        scale_color_manual(values=pop_color_scale)+
        guides(color=guide_legend(title=""))+
        new_scale_color()+ 
        geom_line(data=single, aes(position, admix, color=pop))+
        geom_point(data=single, aes(position, admix, color=pop), size=1)+
        xlim(0,xlim) + 
        theme_bw() +
        ylab("admixture proportion")+
        scale_color_manual(values=pop_color_scale)+
        ggtitle(p_title)+
        theme(plot.title = element_text(size = 8),
              plot.subtitle = element_text(size = 8) ) +
        guides(color=guide_legend(title=""))
    
    nam<-paste0("s_", n)
    assign( nam, p)
    
  }
  
  singles<-plot_grid(plotlist = mget(paste0("s_", 1:length(unique(hybs$ID)))), nrow = length(unique(hybs$ID))) 
  
  ggsave(paste0(PATH,"6.SingleHybrid.pdf"),
         plot = singles,
         device = "pdf",
         scale = 1,
         width = 6,
         height = (4*length(unique(hybs$ID))),
         units = "in",
         limitsize = F)

  #plot admixture windows
  for(o in 1:length(unique(melt_df$region))){
    #get positions from dataframe containing all windows and indv
    melt_hupos<-melt_df[melt_df$region == unique(melt_df$region)[o],]
    
    win<-ggplot(data= melt_hupos, aes(x=admix, y=ID, fill=ancestry))+
          geom_bar(stat="identity")+
          theme_classic()+ 
          ylab("")+ 
          xlab("")+
          ggtitle(unique(melt_hupos$region)) +
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x=element_blank(),
                axis.line.y=element_line(),
                axis.text.y = element_text(size=8, angle=0),
                title = element_text( size=8))+
          guides(fill="none")+
          scale_y_discrete(limits=order)+
          scale_fill_manual(values = pop_color_scale)
    
    nam<-paste0("pl_", o)
    assign( nam, win)
    

  }
  

  test<-plot_grid(plotlist = mget(paste0("pl_", 1:length(unique(melt_df$region)))),  nrow = 1, col=length(unique(melt_df$region)))

  ggsave(paste0(PATH,"7.AdmixturePlots.pdf"),
        plot = test,
        device = "pdf",
        scale = 1,
        width = (3*length(unique(melt_df$region))),
        height = 3,
        units = "in",
        limitsize = F)



  #######All window mean admixture########
  UFmeanAdmix<-data.frame()
  UFhybs<-melt_df_pop1[melt_df_pop1$pop == "hybrid",]
  for(x in 1:length(unique(melt_df_pop1$region))){
    subset_UF<-UFhybs[UFhybs$region == unique(melt_df_pop1$region)[x],]
    subset_UF$meanUF<-sum(subset_UF$admix)/length(subset_UF$admix)
    subset_UF$median<-median(subset_UF$admix)
    subset_UF$chromo<-chromo


    mean_colUF<-cbind(subset_UF$ID,
                    subset_UF$region,
                    subset_UF$position,
                    subset_UF$pop,
                    subset_UF$admix,
                    subset_UF$ancestry,
                    subset_UF$meanUF,
                    subset_UF$median,
                    subset_UF$chromo)

    colnames(mean_colUF)<-c("ID","region","position","pop","admix","ancestry","meanAdmix", "medianAdmix","chrom")
    UFmeanAdmix<-rbind(UFmeanAdmix, mean_colUF)
  }

  UFmeanAdmix$meanAdmix<-as.numeric(UFmeanAdmix$meanAdmix)
  UFmeanAdmix$position<-as.numeric(UFmeanAdmix$position)
    
    ##########################################################################
    FmeanAdmix<-merge(UFmeanAdmix, NA_list, by="region")

    write.csv(FmeanAdmix,paste0(chromo,"_AverageAdmixture.csv"), row.names=F)
    ##########################################################################

  a8_title<-paste0("Average Admixture proportions in unfiltered windows on ", chromo)
  a8<-ggplot()+
      geom_line(data=pure, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=UFhybs, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=UFmeanAdmix, aes(position, meanAdmix))+
      geom_point(data=UFmeanAdmix, aes(position, meanAdmix), size=1)+
      xlim(0,xlim) +
      theme_bw() +
      ylab("admixture proportion")+
      scale_color_manual(values=pop_color_scale)+
      ggtitle(a8_title)+
      theme(plot.title = element_text(size = 8),
            plot.subtitle = element_text(size = 8) ) +
      guides(color=guide_legend(title=""))

  ggsave(paste0(PATH,"8.UnfilteredAverageAdmixtureProportions.png"),
         plot = a8,
         device = "png",
         scale = 1,
         width = 6,
         height = 4,
         units = "in",
         limitsize = F)

  
  ########################################


  print("Finished plotting combined plots with cowplot packages")
  }else{
    print(paste0("No windows present after filtering for fixed windows, plot all windows with cowplot"))

      #plot admixture windows
  for(o in 1:length(unique(melt_df$region))){
    #get positions from dataframe containing all windows and indv
    melt_hupos<-melt_df[melt_df$region == unique(melt_df$region)[o],]

    win<-ggplot(data= melt_hupos, aes(x=admix, y=ID, fill=ancestry))+
          geom_bar(stat="identity")+
          theme_classic()+
          ylab("")+
          xlab("")+
          ggtitle(unique(melt_hupos$region)) +
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x=element_blank(),
                axis.line.y=element_line(),
                axis.text.y = element_text(size=8, angle=0),
                title = element_text( size=8))+
          guides(fill="none")+
          scale_y_discrete(limits=order)+
          scale_fill_manual(values = pop_color_scale)

    nam<-paste0("pl_", o)
    assign( nam, win)


  }


  test<-plot_grid(plotlist = mget(paste0("pl_", 1:length(unique(melt_df$region)))),  nrow = 1, col=length(unique(melt_df$region)))

  ggsave(paste0(PATH,"7.AdmixturePlots.pdf"),
        plot = test,
        device = "pdf",
        scale = 1,
        width = (3*length(unique(melt_df$region))),
        height = 3,
        units = "in",
        limitsize = F)
  print("Finished plotting combined plots with cowplot packages")

 #######All window mean admixture########
  UFmeanAdmix<-data.frame()
  UFhybs<-melt_df_pop1[melt_df_pop1$pop == "hybrid",]
  for(x in 1:length(unique(melt_df_pop1$region))){
    subset_UF<-UFhybs[UFhybs$region == unique(melt_df_pop1$region)[x],]
    subset_UF$meanUF<-sum(subset_UF$admix)/length(subset_UF$admix)
    subset_UF$median<-median(subset_UF$admix)
    subset_UF$chromo<-chromo


    mean_colUF<-cbind(subset_UF$ID,
                    subset_UF$region,
                    subset_UF$position,
                    subset_UF$pop,
                    subset_UF$admix,
                    subset_UF$ancestry,
                    subset_UF$meanUF,
                    subset_UF$median,
                    subset_UF$chromo)

    colnames(mean_colUF)<-c("ID","region","position","pop","admix","ancestry","meanAdmix", "medianAdmix","chrom")
    UFmeanAdmix<-rbind(UFmeanAdmix, mean_colUF)
  }

  UFmeanAdmix$meanAdmix<-as.numeric(UFmeanAdmix$meanAdmix)
  UFmeanAdmix$position<-as.numeric(UFmeanAdmix$position)

  a8_title<-paste0("Average Admixture proportions in unfiltered windows on ", chromo)
  a8<-ggplot()+
      geom_line(data=pure, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=UFhybs, aes(range, admix, color=pop, group=interaction(pop,position, ID)), lwd=3, alpha=0.5, lineend = "round")+
      geom_line(data=UFmeanAdmix, aes(position, meanAdmix))+
      geom_point(data=UFmeanAdmix, aes(position, meanAdmix), size=1)+
      xlim(0,xlim) +
      theme_bw() +
      ylab("admixture proportion")+
      scale_color_manual(values=pop_color_scale)+
      ggtitle(a8_title)+
      theme(plot.title = element_text(size = 8),
            plot.subtitle = element_text(size = 8) ) +
      guides(color=guide_legend(title=""))

  ggsave(paste0(PATH,"8.UnfilteredAverageAdmixtureProportions.png"),
         plot = a8,
         device = "png",
         scale = 1,
         width = 6,
         height = 4,
         units = "in",
         limitsize = F)


  ######################################## 

FmeanAdmix<-merge(UFmeanAdmix, NA_list, by="region")
write.csv(FmeanAdmix,paste0(chromo,"_AverageAdmixture.csv"), row.names=F)
}
print(paste0("Plotting script finished")) 

