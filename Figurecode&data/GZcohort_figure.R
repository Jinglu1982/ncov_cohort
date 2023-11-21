setwd("C:/Users/54437/Desktop/GZcohort")
library(tidyverse)
library(ggplot2)
library(pheatmap)

#### Input cohort data ####
#Ct value
ct_2023throat <- read_xlsx("2023cohort.xlsx",sheet = 1)
ct_2023nose <- read_xlsx("2023cohort.xlsx",sheet = 2)

ct_2022throat <- read_xlsx("2022cohort.xlsx",sheet = 1)
ct_2022nose <- read_xlsx("2022cohort.xlsx",sheet = 2)

#Infection cases
case <- read_xlsx("cases＆heatmap.xlsx",sheet = 1)

#Neutralization titer
titer <- read_xlsx("cases＆heatmap.xlsx",sheet = 2)

#### Data pairing integration ####
#Increase classification
ct_2022nose$type <- c("Nasal swab")
ct_2022nose$year <- c("2022")

ct_2023nose$type <- c("Nasal swab")
ct_2023nose$year <- c("2023")

ct_2022throat$type <- c("Oropharyngeal swab")
ct_2022throat$year <- c("2022")

ct_2023throat$type <- c("Oropharyngeal swab")
ct_2023throat$year <- c("2023")

#Merge four databases
all_ct <- rbind(ct_2022throat[,c(1:5)],ct_2023throat[,c(1,2,4:6)],ct_2022nose[,1:5],ct_2023nose[,c(1,2,4:6)])
all_ct_pos <- filter(all_ct,all_ct$Ct<38)

#### Figure 1B ####
case$`Time periods of onset` <- factor(case$`Time periods of onset`,levels = c("2022/12/7-12/11" ,    "2022/12/12-12/18" ,   "2022/12/19-12/25" ,  
                                                                         "2022/12/26-2023/1/1" ,"2023/4/17-4/23"  ,    "2023/5/1-5/7" ,      
                                                                         "2023/5/8-5/14"   ,    "2023/5/15-5/21"  ,    "2023/5/22-5/28"   ,  
                                                                         "2023/6/5-6/11"   ,    "2023/6/19-6/25"  ,    "2023/7/3-7/9"),ordered = T)

fig_case <- ggplot(case,aes(x=`Time periods of onset`,y=Cases,fill=group,group=1))+
  geom_col(width = 0.6)+
  coord_cartesian(ylim = c(0,25))+
  geom_text(aes(label=Cases,vjust=-0.8,hjust=0.5),size=1.8)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("#8DC4CB","#C1D1EE"))+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 1,angle = 45, hjust = 1),
        legend.position = "none")

fig_case

ggsave("fig_case.png",height = 1.7,width = 7.5)

#### Figure 3A ####
#Intraclass significance(throat vs. nose)
p_val1 <- all_ct_pos %>% 
  group_by(year) %>% 
  wilcox_test(formula = Ct ~ type) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>% 
  add_xy_position(x='year')

#Significance between groups(throat vs. throat/nose vs. nose)
p_val2 <- all_ct_pos %>%
  group_by(type) %>%
  wilcox_test(formula = Ct ~ year) %>%
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>% 
  add_xy_position(x='year')

#Modify the position of p_val2
p_val2$xmin[1] <- 0.8
p_val2$xmin[2] <- 1.2
p_val2$xmax[1] <- 1.8
p_val2$xmax[2] <- 2.2

boxplot_2 <- ggplot(all_ct_pos,aes(year,Ct))+
  geom_boxplot(aes(fill=type,color=type),width=0.5,cex=0.7,
               notch = T,notchwidth = 0.5,key_glyph = draw_key_label,
               position = position_dodge(width =.8))+
  geom_point(aes(fill=type,color=type,shape=type),
             position=position_jitterdodge(jitter.width=0.25,dodge.width=0.7),
             size=1.3,alpha=0.4)+
  scale_color_manual(values = c("#525285","#7D345F"))+
  scale_fill_manual(values = c("#6D6CA3","#963F71"))+
  scale_x_discrete(labels = c("BA.5.2.48-wave", "XBB.1*-wave"))+
  stat_pvalue_manual(p_val1,y.position = -15,tip.length = 0)+
  stat_pvalue_manual(p_val2[1,],y.position = -16,tip.length = 0)+
  stat_pvalue_manual(p_val2[2,],y.position = -16.9,tip.length = 0)+
  scale_y_continuous(breaks=seq(14,38,6),
                     trans = "reverse",labels = c("14","20","26","32","(-)"),
                     sec.axis = sec_axis(~.*-0.2976+14.631,
                                         name = 'Log10 viral loads(copies/ml)'))+
  theme_classic(base_size = 12)+
  labs(x="", y="Ct Value")+
  theme(axis.text.x = element_text(size = 12,face = "bold.italic",hjust=0.5,color="black"),
        axis.line.x = element_blank(),#element_line(color="red",size = 0.8),
        axis.text.y = element_text(size = 9,face = "plain",hjust=0.5),
        axis.line.y = element_blank(),#element_line(color="black",size = 0.8),
        axis.title.y = element_text(size = 10,face = "bold",hjust=0.5),
        panel.spacing.y = unit(0.8, "lines"),
        strip.text = element_text(size = 10,face = "bold",family="Arial"), 
        strip.background = element_rect(color="black",size = 1),
        legend.position = "top",legend.title = element_blank(),legend.text=element_text(size=14),
        legend.background=element_rect(fill='transparent'),
        panel.border = element_rect(color="black",fill = 'transparent',size = 1),
        plot.title = element_text(hjust = 0.5))

boxplot_2

#### Figure 3B ####
#The viral load of BA5
all_2022_1 <- filter(all_ct,all_ct$year==2022,all_ct$DAY >=0)

#Identifying names requires variation
all_2022_1 <- all_2022_1 %>%
  mutate(ID = ifelse(row_number() >= 290, paste0(ID, "a"), ID))

fig_ct_2022 <- all_2022_1 %>% 
  ggplot(aes(x=DAY, y=Ct)) + 
  geom_line(aes(color=type,group=ID), alpha=0.15) + 
  stat_smooth(aes(color=type),alpha=.2,
              method="gam",formula = y~s(x,bs="cs"),
              size=1.1,se=T,fullrange=TRUE)+
  scale_y_reverse(limits=c(40,15),
                  breaks=c(40,35,30,25,20,15),
                  labels=c("40","35","30","25","20","15"))+
  #sec.axis = sec_axis(~.*-0.2976+14.631,name = 'Log10 viral loads(copies/ml)'))+
  scale_x_continuous(breaks=seq(from=0,to=14,by=1)) + 
  scale_color_manual(values = c("#6D6CA3","#963F71"),labels = c("Nasal swab","Oropharyngeal swab"))+
  theme_minimal() + 
  labs(x="Days since Onset of Symptoms", y="Ct Value")+
  theme(legend.position ="top",
        legend.direction = 'vertical',
        legend.title = element_blank(),legend.text=element_text(size=10),
        #Adjust axes
        axis.text.x = element_text(size = 10,face = "bold.italic",hjust=0.1),
        axis.line = element_line(color = "black"),
        axis.text.y = element_text(size = 9,face = "plain",hjust=0.5),
        axis.title.x = element_text(size = 12,face = "bold",hjust=0.5),
        #axis.line.y = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold",hjust=0.5),
        axis.ticks = element_line(size = 0.5, color = "black", linetype = "solid"),
        #Adjustment panel
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 38, colour = 'grey', linetype = 2, size = 0.9)+
  annotate("text",x=14,y=38,hjust=1.5, vjust= 1.5,label="Detection Threshold",size=4,col="black")

fig_ct_2022

#View parameters(fig_ct_2022)
shuzhi_2022 <-ggplot_build(fig_ct_2022)$data[[2]]
shuzhi_2022$y <- -shuzhi_2022$y

shuzhi_2022_nose <- as.data.table(filter(shuzhi_2022,shuzhi_2022$group==1))
shuzhi_2022_throat <- as.data.table(filter(shuzhi_2022,shuzhi_2022$group==2))

#Plot confidence intervals
ribbon_2022 <- shuzhi_2022
ribbon_2022$ymin <- -ribbon_2022$ymin
ribbon_2022$ymax <- -ribbon_2022$ymax
ribbon_2022$y <- -ribbon_2022$y

names(ribbon_2022)[2] <- c("DAY")
names(ribbon_2022)[3] <- c("Ct")
names(ribbon_2022)[1] <- c("type")

ribbon_2022$ymin[is.na(ribbon_2022$ymin)] <- 40

ribbon_2022$type[ribbon_2022$type == "#6D6CA3"] <- "nasal swab"
ribbon_2022$type[ribbon_2022$type == "#963F71"] <- "oropharyngeal swab"

#To regenerate the confidence interval, you need to remove the confidence interval from the original graph（Change“se=F”）
fig_ct_2022_figure2 <- fig_ct_2022+geom_ribbon(
  data = ribbon_2022,aes(x=DAY,ymin=ymin,ymax=ymax,group=type,fill=type),alpha=0.2)+
  scale_fill_manual(values = c("#6D6CA3","#963F71"),
                    labels = c("Nasal swab","Oropharyngeal swab"))+
  annotate("text",x=14,y=15,hjust=1.5, vjust= 1.5,label="BA.5.2.48-wave",size=4.5,col="black",fontface = "bold")

fig_ct_2022_figure2

#### Figure 3C ####
#The viral load of XBB
all_2023_1 <- filter(all_ct,all_ct$year==2023,all_ct$DAY >=0)

#Identifying names requires variation
all_2023_1 <- all_2023_1 %>%
  mutate(ID = ifelse(row_number() >= 78, paste0(ID, "a"), ID))

fig_ct_2023 <- all_2023_1 %>% 
  ggplot(aes(x=DAY, y=Ct)) + 
  geom_line(aes(color=type,group=ID), alpha=0.15) + 
  stat_smooth(aes(color=type),alpha=.2,
              method="gam",formula = y~s(x,bs="cs"),
              size=1.1,se=T,fullrange=TRUE)+
  scale_y_reverse(limits=c(40,15),
                  breaks=c(40,35,30,25,20,15),
                  labels=c("40","35","30","25","20","15"))+
  #sec.axis = sec_axis(~.*-0.2976+14.631,name = 'Log10 viral loads(copies/ml)'))+
  scale_x_continuous(breaks=seq(from=0,to=14,by=1)) + 
  scale_color_manual(values = c("#6D6CA3","#963F71"),labels = c("Nasal swab","Oropharyngeal swab"))+
  theme_minimal() + 
  labs(x="Days since Onset of Symptoms", y="Ct Value")+
  theme(legend.position = c(1,1),legend.justification = c(1,1),
        legend.direction = 'vertical',
        legend.title = element_blank(),legend.text=element_text(size=10),
        #Adjust axes
        axis.text.x = element_text(size = 10,face = "bold.italic",hjust=0.1),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 12,face = "bold",hjust=0.5),
        axis.text.y = element_text(size = 9,face = "plain",hjust=0.5),
        #axis.line.y = element_blank(),
        axis.title.y = element_text(size = 13,face = "bold",hjust=0.5),
        axis.ticks = element_line(size = 0.5, color = "black", linetype = "solid"),
        #Adjust panel
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 38, colour = 'grey', linetype = 2, size = 0.9)+
  annotate("text",x=12,y=38,hjust=1.5, vjust= 1.5,label="Detection Threshold",size=4,col="black")

fig_ct_2023

#View parameters(fig_ct_2023)
shuzhi_2023 <-ggplot_build(fig_ct_2023)$data[[2]]
shuzhi_2023$y <- -shuzhi_2023$y

shuzhi_2023_nose <- as.data.table(filter(shuzhi_2023,shuzhi_2023$group==1))
shuzhi_2023_throat <- as.data.table(filter(shuzhi_2023,shuzhi_2023$group==2))

#Plot confidence intervals
ribbon_2023 <- shuzhi_2023
ribbon_2023$ymin <- -ribbon_2023$ymin
ribbon_2023$ymax <- -ribbon_2023$ymax
ribbon_2023$y <- -ribbon_2023$y

names(ribbon_2023)[2] <- c("DAY")
names(ribbon_2023)[3] <- c("Ct")
names(ribbon_2023)[1] <- c("type")

ribbon_2023$ymin[is.na(ribbon_2023$ymin)] <- 40

ribbon_2023$type[ribbon_2023$type == "#6D6CA3"] <- "nasal swab"
ribbon_2023$type[ribbon_2023$type == "#963F71"] <- "oropharyngeal swab"


#To regenerate the confidence interval, you need to remove the confidence interval from the original graph（Change“se=F”）
fig_ct_2023_figure2 <- fig_ct_2023+geom_ribbon(
  data = ribbon_2023,aes(x=DAY,ymin=ymin,ymax=ymax,group=type,fill=type),alpha=0.2)+
  scale_fill_manual(values = c("#6D6CA3","#963F71"),labels = c("Nasal swab","Oropharyngeal swab"))+
  annotate("text",x=12,y=15,hjust=1.5, vjust= 1.5,label="XBB.1*-wave",size=4.5,col="black",fontface = "bold")

fig_ct_2023_figure2

#### Figure 3 ####
figure_all <- ggarrange(boxplot_2,fig_ct_2022_figure2,fig_ct_2023_figure2
                        ,ncol=3,common.legend = TRUE, legend="top",
                        labels=c("A", "B","C"),font.label = list(size = 17, face = "bold"))

figure_all

#### Figure 4C ####
fig_heatmap <- pheatmap(titer_refined[,c(1:9)],scale = "none",
               color = c(colorRampPalette(colors = c("#4EAC26","white"))(length(bk1)/2),colorRampPalette(colors = c("white","#E375B9"))(length(bk1)/2)),
               fontsize_number = 12,number_color = "black",number_format = "%.0f",
               breaks = bk1,fontsize = 14,
               cellwidth = 40,cellheight = 40,cluster_rows = F,cluster_cols = F,border_color = "black",
               display_numbers = T,
               angle_col = 45,
               main = "People infected with XBB (n=9)")

#### Figure 4D ####
titer <- as.data.frame(titer)
row.names(titer) <- titer[,1]
titer_refined <- titer[7:8,-1]

#Setting the legend range##
bk1 <- c(seq(-6,-0.1,by=0.01),seq(0,6,by=0.01))

#Transforming the data into a matrix##
titer_refined[,c(10:33)] <- as.matrix(titer_refined[,c(10:33)])

fig_heatmap_1 <- pheatmap(titer_refined[,c(10:21)],scale = "none",
                          color = c(colorRampPalette(colors = c("#4EAC26","white"))(length(bk1)/2),colorRampPalette(colors = c("white","#E375B9"))(length(bk1)/2)),
                          breaks = bk1,fontsize = 14,
                          cellwidth = 40,cellheight = 40,cluster_rows = F,cluster_cols = F,border_color = "black",
                          display_numbers = T,fontsize_number = 12,number_color = "black",number_format = "%.0f",
                          angle_col = 45,main = "People uninfected with XBB (n=24)")

fig_heatmap_2 <- pheatmap(titer_refined[,c(22:33)],scale = "none",
                          color = c(colorRampPalette(colors = c("#4EAC26","white"))(length(bk1)/2),colorRampPalette(colors = c("white","#E375B9"))(length(bk1)/2)),
                          
                          breaks = bk1,fontsize_number = 12,number_color = "black",number_format = "%.0f",
                          cellwidth = 40,cellheight = 40,cluster_rows = F,cluster_cols = F,border_color = "black",
                          display_numbers = T,fontsize = 14,
                          angle_col = 45)

ggsave(fig_heatmap,
       filename="heatmap.png",bg = "transparent",
       width = 10,
       height = 2.5,
       dpi = 300)
ggsave(fig_heatmap_1,
       filename="heatmap1.png",bg = "transparent",
       width = 10,
       height = 2.5,
       dpi = 300)
ggsave(fig_heatmap_2,
       filename="heatmap2.png",bg = "transparent",
       width = 10,
       height = 5,
       dpi = 300)
