
library(ape)
library(readxl)
library(data.table)
library(treeio)
library(ggtree)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(ggh4x)
library(stringr)
library(cowplot)

tr_full <- read.tree("time_tree.treefile")
tr_thai <- read.tree("thai_time_tree.treefile")
metadata <- read_xlsx("ZIKV_sequence_metadata.xlsx")
dataset_with_nodes <- readRDS("dataset_with_nodes_phylowave.rds")

# fitness estimates
r_fitness_result_start <- readRDS("time_varying_fitness.rds")
relative_fitness_data_parent <- readRDS("relative_fitness_refparent.rds")
relative_fitness_data_ref2 <- readRDS("relative_fitness_refgroup2.rds")
figuredata_fitness_est <- readRDS("figuredata_fitness_est.rds")
figuredata_fitness_data <- readRDS("figuredata_fitness_data.rds")

source(file = '2_1_Index_computation_20251129.R') 

color1 <-  rev(c("#6D2F20", "#BC7524","#BF3626","#1F5A86","#4E6D58","#2D223C"))
group_num <- 6

# Rename the identified lineages using the labels in `new_group`.
# The analyses were run using the default group names stored in `groups`,
# whereas figures use `new_group` as the lineage labels.
groups_new <- data.table(groups=factor(1:6),new_group=factor(c(4L,6L,5L,3L,2L,1L)))

#### tree preprocess ####
metadata <- as.data.table(metadata)
colnames(metadata)
date_region_info <- metadata[,c(1,2,3,4)]

unique(date_region_info$Country)
removeid <- tr_full$tip.label[which(!tr_full$tip.label%in%tr_thai$tip.label)]
removeid <- removeid[which(!removeid%in%c("NC_035889.1_2015.5",subset(date_region_info,Country%in%c("Brazil"))$id_new) )]

n_b = 6
set.seed(2025)
braz_keep <- c("NC_035889.1_2015.5",sample(subset(date_region_info,(date_num<2017)&(Country%in%c("Brazil"))&!id_new%in%c("NC_035889.1_2015.5"))$id_new,n_b) )
removeid <- c(removeid,subset(date_region_info,(Country%in%c("Brazil"))&!id_new%in%braz_keep)$id_new)
remove_tips <- tr_full$tip.label[which(tr_full$tip.label%in%removeid)]
tr_forfigure <- treeio::drop.tip(tr_full,remove_tips)

#### nodedata preprocess ####
nodedata_tmp <- as.data.frame(dataset_with_nodes)

tip_nodedata <- subset(nodedata_tmp,name_seq%in%tr_forfigure$tip.label)
tip_nodedata <- as.data.table(tip_nodedata)

nontip_nodedata <- data.table(name_seq=length(tr_forfigure$tip.label)+(1:(length(tr_forfigure$tip.label)-1)),
                              is.node="yes",
                              clade_info=NA)
tip_nodedata_tmp <- data.table(name_seq=tr_forfigure$tip.label[which(!tr_forfigure$tip.label%in%tip_nodedata$name_seq)],
                               is.node="no",
                               clade_info=NA)

genetic_distance_mat = dist.nodes.with.names(tr_forfigure)
nroot = length(tr_forfigure$tip.label) + 1 ## Root number
distance_to_root =genetic_distance_mat[nroot,]
root_height = tip_nodedata[which(tip_nodedata[,name_seq == dataset_with_nodes$name_seq[1] ]),]$time - distance_to_root[which(names(distance_to_root) == dataset_with_nodes$name_seq[1])]
nodes_height = root_height + distance_to_root[length(tr_forfigure$tip.label)+(1:(length(tr_forfigure$tip.label)-1))]

nontip_nodedata[,time:=nodes_height]
nontip_nodedata[,time1:=round(time,4)]

date_region_info_tmp <- date_region_info[,2:3]
colnames(date_region_info_tmp)[1] <- "name_seq"
colnames(date_region_info_tmp)[2] <- "time"
tip_nodedata_tmp <- merge(tip_nodedata_tmp,date_region_info_tmp,by="name_seq")

for_match <- subset(nodedata_tmp,is.node=="yes")
for_match <- as.data.table(for_match)
for_match[,time1:=round(time,4)]
colnames(nontip_nodedata)
colnames(for_match)

nontip_nodedata <- merge(nontip_nodedata,for_match[,7:9],by="time1",all.x=TRUE)
nontip_nodedata <- nontip_nodedata[,-1]
nontip_nodedata[,ID:=name_seq]

tip_id <- data.table(name_seq = tr_forfigure$tip.label,ID = 1:length(tr_forfigure$tip.label))
tip_nodedata_tmp$index=NA
tip_nodedata_tmp$groups=NA
tip_nodedata_tmp$new_group=NA
tip_nodedata <- merge(rbind(tip_nodedata[,-5],tip_nodedata_tmp),tip_id,by="name_seq")
nontip_nodedata[,index:=NA]

dataset_with_nodes_forfigure <-  rbind(tip_nodedata,nontip_nodedata)
dataset_with_nodes_forfigure <- as.data.table(dataset_with_nodes_forfigure)
dataset_with_nodes_forfigure <- dataset_with_nodes_forfigure[order(dataset_with_nodes_forfigure[,ID]),]
colnames(dataset_with_nodes_forfigure)
colnames(nodedata_tmp)
dataset_with_nodes_forfigure <- dataset_with_nodes_forfigure[,c(1:4,8,5:7)]

#### figure 1 ####
f1 <- as.data.frame( dataset_with_nodes_forfigure )
f1 <- as.data.table(f1)
f1$node <- f1$ID

# to merge in the node on the branch to group 4, for figure linetype for brazil sequences
set(f1,which(f1[,ID==133]),"groups",factor(4,levels = 1:6))

p_tree1 <- ggtree(tr_forfigure,mrsd=lubridate::date_decimal(max(f1$time))) %<+% f1 +
  geom_hilight(node=143, fill="#2C6582", alpha=0.6,type = "roundrect")+
  geom_point2(aes(subset=ID%in%c(109)),shape=10,size=3,stroke = 1)+
  geom_point2(aes(subset=ID%in%c(143)),shape=8,size=3,stroke = 1)+
  aes( linetype = is.na(groups) ) +
  geom_tippoint(aes(subset=!is.na(groups)),size = 2,color = "#A6767A") +
  theme_tree2()+
  scale_x_continuous(breaks = seq(2000,2024,by=5),limits = c(1999.5,2024),expand = c(0,0))+
  theme(legend.position = c(0.1,0.8))+
  guides(linetype = "none")

p_tree <- p_tree1 + 
  annotate("text",label="Brazil",x=2017.8,y=26.8,color = "#164464",fontface="bold",size=3) +
  annotate("point",shape=10,x=2002,y=90,size=3) +
  annotate("text",label="Feb 2000 [Oct 1998-Apr 2001]",x=2002+5.5,y=90+1,size=2.5) +
  annotate("text",label="PP=1.00",x=2002+2.4,y=90-3+1,size=2.5) +
  annotate("point",shape=8,x=2002,y=80,size=3) +
  annotate("text",label="Dec 2012 [Jul 2012-Apr 2013]",x=2002+5.4,y=80+1,size=2.5)+
  annotate("text",label="PP=0.95",x=2002+2.4,y=80-3+1,size=2.5)

stopifnot(dim(subset(date_region_info,id_new%in%tr_thai$tip.label))[1] == length(tr_thai$tip.label))
p_samplesize <- ggplot(NULL)+
  geom_bar(data=subset(date_region_info,id_new%in%tr_thai$tip.label),aes(x=floor(date_num)),
           alpha=0.95,stat="count",color="grey40",fill="#D5B5BB")+
  scale_x_continuous(breaks = seq(2000,2024,by=5),limits = c(1999.5,2024),expand = c(0,0))+
  scale_y_continuous(limits = c(0,30),breaks=seq(0,28,by=5),expand = c(0,0))+
  labs(x="Time (year)",y="Number of sequences",title = "Sequence collection by time")+
  theme_bw()+
  theme(axis.line = element_line(),panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        plot.title = element_text(size=11),
        legend.key.size = unit(0.5, 'cm'),
        strip.background = element_rect(fill = "white"))

psave <- ggarrange(p_samplesize,p_tree,ncol=2,widths = c(1,1),labels = c("a","b"))
ggsave("figure1.png",psave,width = 8,height = 4.5,bg="white")

#### figure 2 ####
t_start_data <- NULL
t_end_data <- NULL

for (k in 1:group_num){
  t_start_this <- min(subset(dataset_with_nodes,groups == k)$time)
  
  t_end_this <- max(subset(dataset_with_nodes,groups==k)$time)
  t_start_data1 <- data.frame(groups=k,t_start_m=t_start_this)
  t_end_data1 <- data.frame(groups=k,t_end_m=t_end_this)
  t_start_data <- rbind(t_start_data,t_start_data1)
  t_end_data <- rbind(t_end_data,t_end_data1)
}

t_start_data <- merge(t_start_data,groups_new,by="groups")
t_end_data <- merge(t_end_data,groups_new,by="groups")

##### the tree figure #####
p_tree1 <- ggtree(tr_forfigure,mrsd=lubridate::date_decimal(max(f1$time))) %<+% f1 +
  geom_hilight(node=143, fill="#2C6582", alpha=0.6,type = "roundrect")+
  aes( color = new_group,linetype = is.na(groups) ) +
  geom_point2(aes(subset = ((!is.na(groups)) & (!ID%in%c(132,133))), color = new_group), size = 2) +
  geom_point2(aes(subset=ID%in%c(132)),shape=0,size=3,stroke = 1,color=color1[3])+
  geom_point2(aes(subset=ID%in%c(133)),shape=2,size=3,stroke = 1,color=color1[4])+
  scale_color_manual(name="Group",
                     breaks = 1:6,
                     values = color1,na.value = 'black')+
  theme_tree2()+
  scale_x_continuous(breaks = seq(2000,2024,by=5),limits = c(1999.5,2024),expand = c(0,0))+
  theme(legend.position = c(0.1,0.7))+
  guides(linetype = "none")

p_tree <- p_tree1 + 
  annotate("text",label="Brazil",x=2017.8,y=26.8,color = "#164464",fontface="bold",size=3)+
  annotate("point",shape=0,size=3,color=color1[3],x=2007,y=85) +
  annotate("text",label="PP=1.00",x=2007+2.4,y=85,size=3) +
  annotate("point",shape=2,size=3,color=color1[4],x=2007,y=80) +
  annotate("text",label="PP=0.77",x=2007+2.4,y=80,size=3)

p_index <- ggplot(dataset_with_nodes)+
  geom_point(aes(x=time,y=index,color=new_group,group=new_group),size=0.8)+
  scale_color_manual(name="",breaks = 1:6,
                     values = color1)+
  scale_y_continuous(breaks = seq(0,1,by=0.2),limits = c(0.35,1.01),expand = c(0,0))+
  scale_x_continuous(breaks = seq(2000,2024,by=5),limits = c(1999.5,2024),expand = c(0,0))+
  labs(x = "Time (year)",y="Index")+
  theme_bw()+
  theme(legend.position = "",axis.line = element_line(),panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))

##### fitness figures #####
p_fitness <- ggplot()+
  geom_line(data=figuredata_fitness_est,aes(x=t,y=est_m,color=factor(new_group),group=factor(new_group)))+
  geom_ribbon(data=figuredata_fitness_est,aes(x=t,ymin=est_l,ymax=est_u,fill=factor(new_group),group=factor(new_group)),alpha=0.5)+
  geom_point(data=figuredata_fitness_data,aes(x=t,y=data_m,color=factor(new_group),group=factor(new_group)),size=0.8)+
  geom_errorbar(data=figuredata_fitness_data,aes(x=t,ymin=data_l,ymax=data_u,color=factor(new_group),group=factor(new_group)),width=0.2)+
  facet_wrap(.~factor(new_group),ncol=3,scales = "free",
             labeller = as_labeller(c("1"="Group 1","2"="Group 2",
                                      "3"="Group 3","4"="Group 4","5"="Group 5","6"="Group 6")))+
  scale_color_manual(name="",
                     breaks=1:group_num,
                     values=color1)+
  scale_fill_manual(name="",
                    breaks=1:group_num,
                    values=color1)+
  scale_x_continuous(breaks = seq(2000,2024,by=5),limits = c(1999.5,2024),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.02,1.02),breaks = seq(0,1,by=0.2),expand = c(0,0))+
  labs(x="Time (year)",y="Proportion",title="Fitness of lineage")+
  theme_bw()+
  theme(legend.position = "",axis.line = element_line(),panel.border = element_blank(),
        strip.text.x = element_text(hjust = 0,size = 11),
        plot.title = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())

p_beta_parent <- ggplot()+
  geom_point(data=subset(relative_fitness_data_parent,new_group!=1),aes(x=factor(new_group),y=exp(est_m),color=factor(new_group),group=factor(new_group)))+
  geom_errorbar(data=subset(relative_fitness_data_parent,new_group!=1),aes(x=factor(new_group),ymin=exp(est_l),ymax=exp(est_u),color=factor(new_group),group=factor(new_group)),width=0.2)+
  scale_color_manual(name="",
                     breaks=1:group_num,
                     values=color1)+
  scale_fill_manual(name="",
                    breaks=1:group_num,
                    values=color1)+
  scale_x_discrete(limits = factor(1:6))+
  scale_y_continuous(trans = "log",breaks = c(1,2,4,6),limits = c(0.9,8))+
  geom_abline(slope=0,intercept = 0,linetype="dashed")+
  coord_flip()+
  labs(x="Groups",y="Relative fitness growth rate\ncompared to the parental lineage")+
  theme_bw()+
  theme(legend.position = "",axis.line = element_line(),panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=9),
        strip.background = element_rect(fill = "white"))
p_beta_parent <- p_beta_parent + 
  annotate("text",label="Not applicable",x=1,y=1.5,size=2.5)

p_beta_ref2 <- ggplot()+
  geom_point(data=subset(relative_fitness_data_ref2,new_group!=1),aes(x=factor(new_group),y=exp(est_m),color=factor(new_group),group=factor(new_group)))+
  geom_errorbar(data=subset(relative_fitness_data_ref2,!new_group%in%c(1,2)),aes(x=factor(new_group),ymin=exp(est_l),ymax=exp(est_u),color=factor(new_group),group=factor(new_group)),width=0.2)+
  scale_color_manual(name="",
                     breaks=1:group_num,
                     values=color1)+
  scale_fill_manual(name="",
                    breaks=1:group_num,
                    values=color1)+
  scale_x_discrete(limits = factor(1:6))+
  scale_y_continuous(trans = "log",breaks = c(1,2,4,6,8,10),limits = c(0.9,10))+
  geom_abline(slope=0,intercept = 0,linetype="dashed")+
  coord_flip()+
  labs(x="Groups",y="Relative fitness growth rate\ncompared to Group 2")+
  theme_bw()+
  theme(legend.position = "",axis.line = element_line(),panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=9),
        strip.background = element_rect(fill = "white"))
p_beta_ref2 <- p_beta_ref2 + 
  annotate("text",label="Not applicable",x=1,y=1.58,size=2.5)+
  annotate("text",label="Reference",x=2,y=1.42,size=2.5)

pal = c(colorspace::sequential_hcl(8, palette = "Blues"), 'white', rev(colorspace::sequential_hcl(8, palette = "OrRd")))

p_relative_fitness_time <- ggplot(r_fitness_result_start)+
  geom_tile(aes(x=time,y=factor(new_group),fill=exp(r_fitness_m)))+
  scale_fill_gradientn(colours = pal, na.value = 'grey92',
                       limits = c(1/1.28,1.25), trans = "log",
                       breaks = (c(0.8,1,1.20)),name="Relative\nfitness")+
  geom_point(data=t_start_data,aes(x=t_start_m,y=factor(new_group)),size=1.5,shape=2,stroke = 1)+
  geom_point(data=t_end_data,aes(x=t_end_m,y=factor(new_group)),size=1.5,shape=1,stroke = 1)+
  scale_x_continuous(breaks = seq(2000,2024,by=5),limits = c(1999.5,2024),expand = c(0,0))+
  labs(x="Time (year)",y="Group")+
  theme_bw()+
  theme(legend.background = element_blank(), legend.key = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), plot.background = element_blank(), 
        legend.title = element_text(size=9),
        axis.line = element_line(), panel.grid = element_blank())

p_relative_fitness_time <- p_relative_fitness_time + annotate('rect', xmin = 2015.01, xmax = 2016.99, ymin = 0.5, ymax = 6.5,
                                                              fill = NA, color = "#2C6582", linewidth = 1)

ptree_index <- ggarrange(p_tree,p_index,ncol=1,heights = c(2,1),labels = c("a","b"),align = "v")
p_relative_f <- ggarrange(p_beta_ref2,p_beta_parent,p_relative_fitness_time,ncol=3,labels = c("d","e","f"),widths = c(1,1,1.3),align = "h")

p_fitness_1 <- ggarrange(p_fitness,labels = "c")
p_fitness_combine <- ggarrange(p_fitness_1,p_relative_f,ncol=1,heights = c(2,1))
p_all <- ggarrange(ptree_index,p_fitness_combine,ncol=2,widths = c(1.1,2))

ggsave("figure2.png",p_all,width = 12.5,height = 8,bg="white")

#### figure 3 ####

aa_table <- readRDS("aa_mutation_table.rds")
utr_table <- readRDS("utr_mutation_table.rds")
aa_table_columnname <- readRDS("aa_table_columnname.rds")
utr_table_columnname <- readRDS("utr_table_columnname.rds")

geno_name_aa <- as.data.frame(aa_table_columnname)
geno_name1_aa <- as.data.frame(aa_table_columnname)
geno_name_aa$new_group <- "Gene1"
geno_name1_aa$new_group <- "Gene2"

geno_name_utr <- as.data.frame(utr_table_columnname)
geno_name1_utr <- as.data.frame(utr_table_columnname)
geno_name_utr$new_group <- "Gene1"
geno_name1_utr$new_group <- "Gene2"

geno_name_aa <- as.data.table(geno_name_aa)
geno_name1_aa <- as.data.table(geno_name1_aa)
geno_name_utr <- as.data.table(geno_name_utr)
geno_name1_utr <- as.data.table(geno_name1_utr)

y_order <- c(levels(unique(utr_table$mutation_name1)),levels(unique(aa_table$mutation_name1)))
utr_table <- as.data.table(utr_table)
set(utr_table,NULL,"mutation_name1",utr_table[,factor(mutation_name1,levels = y_order)])
aa_table <- as.data.table(aa_table)
set(aa_table,NULL,"mutation_name1",aa_table[,factor(mutation_name1,levels = y_order)])

set(geno_name1_utr,NULL,"mutation_name1",geno_name1_utr[,factor(mutation_name1,levels = y_order)])
set(geno_name_utr,NULL,"mutation_name1",geno_name_utr[,factor(mutation_name1,levels = y_order)])
set(geno_name1_aa,NULL,"mutation_name1",geno_name1_aa[,factor(mutation_name1,levels = y_order)])
set(geno_name_aa,NULL,"mutation_name1",geno_name_aa[,factor(mutation_name1,levels = y_order)])

posdata_figure <- rbind(utr_table,
                        aa_table)
set(posdata_figure,NULL,"mutation_name1",posdata_figure[,factor(mutation_name1,levels = y_order)])

geno_name_utr[,protein:="3′ UTR"]
geno_name_figure <- rbind(geno_name_aa[,c(4,7,2,5,6)],geno_name_utr[,c(1,5,6,3,4)])

colnames(geno_name1_utr)[2] = "pos1"
geno_name1_figure <- rbind(geno_name1_aa[,c(4,7,3,5,6)],geno_name1_utr)

##### table #####
p_mutation_table <- ggplot()+
  geom_tile(data=posdata_figure,
            aes(y=mutation_name1,x=new_group,fill=new_color_group),alpha=0.9,color="white")+
  geom_text(data=geno_name_figure,
            aes(y=mutation_name1,x=new_group,label=protein),color="black",size=5)+
  geom_text(data=geno_name1_figure,
            aes(y=mutation_name1,x=new_group,label=pos1),color="black",size=5)+
  geom_text(data=subset(posdata_figure,name1==name_after),
            aes(y=mutation_name1,x=new_group,label=str_to_upper(name1)),color="white",fontface="bold",size=5)+
  geom_text(data=subset(posdata_figure,name1!=name_after),
            aes(y=mutation_name1,x=new_group,label=str_to_upper(name1)),color="black",fontface="bold",size=5)+
  facet_grid(facet2~facet1,scales="free")+
  force_panelsizes(cols = c(8,1,1),rows=c(17,4))+
  facetted_pos_scales(x=list(
    facet1 == 1 ~ scale_x_discrete(position = "top",
                                   limits = c("Gene1","Gene2",as.character(1:group_num)),
                                   breaks = c("Gene1","Gene2",as.character(1:group_num)),
                                   label = c(" \n "," \n ",as.character(1:group_num))),
    facet1 == 3 ~ scale_x_discrete(position = "top",
                                   limits = c("brazil"),
                                   breaks = c("brazil"),
                                   label = c("Brazil")),
    facet1 == 2 ~ scale_x_discrete(position = "top",
                                   limits = c("fp"),
                                   breaks = c("fp"),
                                   label = c("FP"))
  ))+
  scale_fill_manual(name="",
                    breaks=c("no",as.character(1:group_num),"brazil"),
                    values=c("grey80",color1,"#5B5878"))+
  theme_bw()+
  theme(legend.position = "",axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size=11),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_blank(),panel.spacing.y=unit(0.05, "lines"))

p_mutation_table1 <- ggdraw(p_mutation_table) +
  draw_line(x = c(0.215, 0.775), y = c(0.975, 0.975), linewidth = 0.5)+
  draw_text("Group",x = mean(c(0.21, 0.775)), y = 0.988,size=14)+
  draw_text("Genome\nregion",x = 0.068, y = 0.975,size=14) +
  draw_text("Position",x = 0.163, y = 0.975,size=14)

##### simplified tree #####
keep_id <- NULL
for (k in 1:group_num){
  d1 <- subset(dataset_with_nodes,groups==k)
  id_this <- d1[which.max(d1[,time]),]$name_seq
  keep_id <- c(keep_id,id_this)
}

keep_id <- c(keep_id,"NC_035889.1_2015.5","KX447511.1_2014.042")

remove_tips_simplified <- tr_full$tip.label[which(!tr_full$tip.label%in%keep_id)]

tr_sub_simplified <- treeio::drop.tip(tr_full,remove_tips_simplified)

sim_tip <- subset(dataset_with_nodes,name_seq %in% keep_id)

genetic_distance_mat_sim = dist.nodes.with.names(tr_sub_simplified)
sim_nontip <- data.table(name_seq=length(tr_sub_simplified$tip.label)+(1:(length(tr_sub_simplified$tip.label)-1)),
                         is.node="yes",
                         clade_info=NA)
nroot = length(tr_sub_simplified$tip.label) + 1 ## Root number
distance_to_root_sim =genetic_distance_mat_sim[nroot,]
root_height = sim_tip[which(sim_tip[,name_seq == names(distance_to_root_sim[2])]),]$time - distance_to_root_sim[2]
nodes_height = root_height + distance_to_root_sim[(group_num+2)+(1:(group_num+2-1))]
sim_nontip[,time:=nodes_height]
sim_nontip[,time1:=round(time,4)]

for_match_sim <- subset(dataset_with_nodes,is.node=="yes")
for_match_sim[,time1:=round(time,4)]

sim_nontip <- merge(sim_nontip,for_match_sim[,7:9],by="time1",all.x=TRUE)
sim_nontip <- sim_nontip[,-1]
sim_nontip[,ID:=name_seq]

tip_id <- data.table(name_seq = tr_sub_simplified$tip.label,ID = 1:(group_num+2))
sim_tip <- merge(sim_tip[,-5],tip_id,by="name_seq")
sim_nontip[,index:=NA]

nodedata_sim <- rbind(sim_tip,sim_nontip)
nodedata_sim <- as.data.table(nodedata_sim)
nodedata_sim <- nodedata_sim[order(nodedata_sim[,ID]),]

f1 <- as.data.frame( nodedata_sim )
f1 <- as.data.table(nodedata_sim)
f1$node <- f1$ID

set(f1,which(f1[,ID==12]),"new_group",factor(4))

psimtree <- ggtree(tr_sub_simplified,layout = "roundrect",mrsd=lubridate::date_decimal(max(f1$time))) %<+% f1 +
  aes( color = new_group,linetype = is.na(new_group) ) +
  geom_point2(aes(color = new_group), size = 2) +
  scale_color_manual(breaks = 1:6,
                     values = color1,na.value = "grey40")+
  theme(legend.margin=margin(c(t=0,r=0,b=0,l=0)))+
  guides(color="none",linetype = "none")

psimtree1 <- ggarrange(psimtree,NULL,ncol=2,widths = c(1,0.2))
psimtree2 <- ggdraw(psimtree1) +
  draw_text("Schematic phylogeny of lineages",x = 0.4, y = 0.975,size=14)+
  draw_text("Group 1",size=14,color=color1[1],x=0.36,y=0.1) +
  draw_text("Group 2",size=14,color=color1[2],x=0.73,y=0.1+0.116) +
  draw_text("Group 4",size=14,color=color1[4],x=0.87,y=0.1+0.116*2) +
  draw_text("Brazil",size=14,color="grey40",x=0.62,y=0.1+0.116*3) +
  draw_text("FP",size=14,color="grey40",x=0.54,y=0.1+0.116*4) +
  draw_text("Group 3",size=14,color=color1[3],x=0.7,y=0.1+0.116*5) +
  draw_text("Group 5",size=14,color=color1[5],x=0.865,y=0.1+0.116*6) +
  draw_text("Group 6",size=14,color=color1[6],x=0.87,y=0.1+0.116*7) 

psimtree3 <- ggarrange(NULL,psimtree2,NULL,nrow=3,heights = c(0.07,1.6,0.93))
psave <- ggarrange(p_mutation_table1,psimtree3,ncol=2,widths = c(1,0.5))

ggsave("figure3.png",psave,width = 12.75,height = 9,bg="white")

