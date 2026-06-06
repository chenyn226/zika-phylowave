
library(phytools)
library(data.table)

#### load function from phylowave code base ####
# reference: Lefrancq, Noémie, et al. "Learning the fitness dynamics of pathogens from phylogenies." Nature 637.8046 (2025): 683-690.
source(file = '2_1_Index_computation_20251129.R') 
source(file = '2_2_Lineage_detection_20260127.R') 

#### parameter setting ####
min_descendants_per_tested_node = 10
min_group_size = 10
group_count_threshold= 15
timescale = 5
wind = 2

genome_length = 11171
mutation_rate = 6.8e-4 

time_window_initial = 2030
time_window_increment = 100
p_value_smooth = 0.05
weight_by_time = 0.5
k_smooth = 3

plot_screening = F
weighting_transformation = c('inv_sqrt')
parallelize_code = F
number_cores = 1

max_stepwise_deviance_explained_threshold = 0.005
max_groups_found = 30
stepwise_AIC_threshold = 0
keep_track = T

#### lineage detection algorithm ####
tr1 <- read.tree("thai_time_tree.treefile")

##### index calculation #####

genetic_distance_mat = dist.nodes.with.names(tr1)

nodedata1_tips <- data.table(name_seq=tr1$tip.label,
                             is_node="no",
                             clade_info=NA)

k=1
tmp1 <- strsplit(nodedata1_tips[k,]$name_seq,"_")[[1]]
if (length(tmp1)==2){
  set(nodedata1_tips,k,"time",as.numeric(tmp1[2]))
}
if (length(tmp1)==3){
  set(nodedata1_tips,k,"time",as.numeric(tmp1[3]))
}
for (k in 2:dim(nodedata1_tips)[1]){
  tmp1 <- strsplit(nodedata1_tips[k,]$name_seq,"_")[[1]]
  if (length(tmp1)==2){
    set(nodedata1_tips,k,"time",as.numeric(tmp1[2]))
  }
  if (length(tmp1)==3){
    set(nodedata1_tips,k,"time",as.numeric(tmp1[3]))
  }
}


times_seqs <- nodedata1_tips$time
n_seq <- length(tr1$tip.label)

nodedata1_nontips <- data.table(name_seq=length(tr1$tip.label)+(1:(length(tr1$tip.label)-1)),
                                is_node="yes",
                                clade_info=NA)
nroot = length(tr1$tip.label) + 1 ## Root number
distance_to_root = genetic_distance_mat[nroot,]
root_height = nodedata1_tips[which(nodedata1_tips[,name_seq == names(distance_to_root[1])]),]$time - distance_to_root[1]
nodes_height = root_height + distance_to_root[n_seq+(1:(n_seq-1))]

nodedata1_nontips[,time:=nodes_height]

nodedata1 <- rbind(nodedata1_tips,nodedata1_nontips)
nodedata1[,ID:=1:dim(nodedata1)[1]]
index1 <- compute.index(time_distance_mat = genetic_distance_mat,
                        timed_tree = tr1,
                        time_window = wind,
                        metadata = nodedata1,
                        mutation_rate = mutation_rate,
                        timescale = timescale,
                        genome_length = genome_length)

nodedata1$index <- index1


stopifnot(nodedata1$name_seq==names(index1))

##### find clade based on index dynamics #####
edge1 <- tr1$edge.length
location_index <- match((n_seq+1):(2*n_seq-1),tr1$edge[,2])
edge1 <- edge1[location_index]

colnames(nodedata1)
colnames(nodedata1)[2] = "is.node"

set(nodedata1,which(nodedata1[,index==0]),"index",0.0000000001)

nodedata1 <- as.data.table(nodedata1)
nodedata1 <- nodedata1[order(nodedata1[,ID]),]

potential_splits = 
  find.groups.by.index.dynamics(timed_tree = tr1,
                                metadata = as.data.frame(nodedata1),
                                node_support = edge1,
                                threshold_node_support = 1/(genome_length*mutation_rate),
                                time_window_initial = time_window_initial, 
                                time_window_increment = time_window_increment,
                                min_descendants_per_tested_node =min_descendants_per_tested_node,
                                min_group_size = min_group_size,
                                p_value_smooth = p_value_smooth,
                                stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                stepwise_AIC_threshold = stepwise_AIC_threshold, 
                                weight_by_time = weight_by_time, 
                                weighting_transformation = weighting_transformation,
                                k_smooth = k_smooth,
                                parallelize_code = parallelize_code, 
                                number_cores = number_cores, 
                                plot_screening = plot_screening, 
                                max_groups_found = max_groups_found, 
                                keep_track = keep_track)


df_explained_dev = data.frame('N_groups' = 0:length(potential_splits$best_dev_explained), 
                              'Non_explained_deviance' = (1-c(potential_splits$first_dev,
                                                              potential_splits$best_dev_explained)),
                              'Non_explained_deviance_log' = log(1-c(potential_splits$first_dev,
                                                                     potential_splits$best_dev_explained)))
df_explained_dev$Non_explained_deviance_log = df_explained_dev$Non_explained_deviance_log- min(df_explained_dev$Non_explained_deviance_log)

# optimize the number of groups

split <- merge.groups(timed_tree = tr1,metadata = as.data.frame(nodedata1),
                      initial_splits = potential_splits$potential_splits,
                      group_count_threshold = group_count_threshold, # 30 in the example 
                      group_freq_threshold = 0.01)

##### get the results #####
nodedata1[,groups:=factor(split$groups)]
name_groups = levels(nodedata1$groups)


i = 1
time_groups_world = min(nodedata1[which(nodedata1[,groups==name_groups[i]]),]$time)

for (i in 2:length(name_groups)){
  time_groups_world = c(time_groups_world,
                        min(nodedata1[which(nodedata1[,groups==name_groups[i]]),]$time))
}

# reorder group number by time of emerge
group_time_data <- data.table(group1=1:length(name_groups),grouptime=time_groups_world)
group_time_data <- group_time_data[order(group_time_data[,grouptime],decreasing = TRUE),]
group_time_data[,group2:=1:dim(group_time_data)[1]]
group_time_data <- group_time_data[order(group_time_data[,group1]),]

colnames(nodedata1)[7] = "group1"
set(group_time_data,NULL,"group1",group_time_data[,factor(group1)])
nodedata1 <- merge(nodedata1,group_time_data[,c(1,3)],by="group1")
nodedata1 <- nodedata1[,-1]
colnames(nodedata1)[7] <- "groups"
set(nodedata1,NULL,"groups",nodedata1[,as.factor(groups)])
levels(nodedata1$groups)

split_nodegroup <- data.frame(group1=split$tip_and_nodes_groups)
split_nodegroup$orderindex=1:dim(split_nodegroup)[1]
split_nodegroup <- merge(split_nodegroup,group_time_data[,c(1,3)],by="group1",all.x = TRUE)
split_nodegroup <- as.data.table(split_nodegroup)
split_nodegroup <- split_nodegroup[order(split_nodegroup[,orderindex])]
split$tip_and_nodes_groups <- split_nodegroup$group2
names(split$tip_and_nodes_groups) <- 1:length(split$tip_and_nodes_groups)

split$groups = as.factor(split$groups)
split_groupdata <- data.frame(group1=split$groups)
split_groupdata$orderindex = 1:dim(split_groupdata)[1]
split_groupdata <- merge(split_groupdata,group_time_data,by="group1")
split_groupdata <- as.data.table(split_groupdata)
split_groupdata <- split_groupdata[order(split_groupdata[,orderindex]),]
summary(split_groupdata)
split$groups <- split_groupdata$group2

# saveRDS(split,"split_phylowave.rds")
# saveRDS(nodedata1,"dataset_with_nodes_phylowave.rds")
