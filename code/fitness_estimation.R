
library(phytools)
library(data.table)
library(cmdstanr)
library(stringr)

tr1 <- read.tree("thai_time_tree.treefile")
nodedata <- readRDS("dataset_with_nodes_phylowave.rds")

min_year1 = 2000
window1 = 1

# function for stan model implementation
read.chains.from.table = function(table){
  if(typeof(table) != 'double') {
    print('Changing type of table to matrix')
    table = as.matrix(table)
  }
  Chains = list()
  col_variables = sapply(colnames(table), function(x)str_split(x, pattern = "[.]")[[1]][1])
  variable_names = unique(col_variables)
  nchains = nrow(table)
  for(i in 1:length(variable_names)){
    a = match(col_variables, variable_names[i])
    if(length(which(is.na(a) == F)) == 1) {
      Chains[[i]] = table[,match(variable_names[i], col_variables)]
    }
    else{
      a = match(col_variables, variable_names[i])
      tmp = colnames(table)[which(is.na(a) == F)]
      ndims = 0
      dims = NULL
      empty = F
      while(empty == F & ndims < 10){
        l = sapply(tmp, function(x)str_split(x, pattern = "[.]")[[1]][1+ndims+1])
        if(all(is.na(l))) empty = T
        if(all(is.na(l)) == F) {
          ndims = ndims + 1
          dims = c(dims, max(as.numeric(l)))
        }
      }
      if(ndims >5) print('Error: this function only supports arrays of <=5 dimensions')
      Chains[[i]] = array(NA, dim = c(nchains, dims))
      a = which(is.na(match(col_variables, variable_names[i]))==F)
      
      if(ndims == 1){
        Chains[[i]] = table[,a]
        colnames(Chains[[i]]) = NULL
      }
      
      if(ndims == 2){
        j = 1
        for(d2 in 1:dims[2]){
          for(d1 in 1:dims[1]){
            Chains[[i]][,d1,d2] = table[,a[j]]
            j = j+1
          }
        }
      }
      
      if(ndims == 3){
        j = 1
        for(d3 in 1:dims[3]){
          for(d2 in 1:dims[2]){
            for(d1 in 1:dims[1]){
              Chains[[i]][,d1,d2,d3] = table[,a[j]]
              j = j+1
            }
          }
        }
      }
      
      if(ndims == 4){
        j = 1
        for(d4 in 1:dims[4]){
          for(d3 in 1:dims[3]){
            for(d2 in 1:dims[2]){
              for(d1 in 1:dims[1]){
                Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                j = j+1
              }
            }
          }
        }
      }
      
      if(ndims == 5){
        j = 1
        for(d5 in 1:dims[5]){
          for(d4 in 1:dims[4]){
            for(d3 in 1:dims[3]){
              for(d2 in 1:dims[2]){
                for(d1 in 1:dims[1]){
                  Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                  j = j+1
                }
              }
            }
          }
        }
      }
    }
    names(Chains)[i] = variable_names[i]
  }
  return(Chains)
}

estimate_rel_fitness_groups_with_branches_tstart = function(dataset_with_nodes, tree, min_year = 1950, window = NULL, N = NULL, 
                                                            model_compiled, iter_warmup = 250, iter_sampling = 500, refresh = 50, seed = 1){
  ## Retrieve branches from tree
  branches = tree$edge
  
  ## Branch times
  branches_times = tree$edge
  branches_times[,1] = branches_times[,2] = NA
  times_tips_nodes = dataset_with_nodes$time
  branches_times[,1] = times_tips_nodes[match(branches[,1], dataset_with_nodes$ID)]
  branches_times[,2] = times_tips_nodes[match(branches[,2], dataset_with_nodes$ID)]
  
  ## Branch groups
  branches_group = tree$edge
  branches_group[,1] = branches_group[,2] = NA
  branches_group[,1] = as.numeric(dataset_with_nodes$groups[branches[,1]])
  branches_group[,2] = as.numeric(dataset_with_nodes$groups[branches[,2]])
  
  ## Group names
  groups = table(factor(dataset_with_nodes$groups))
  groups = names(groups)
  
  ## Set time window, to count number of sequences and nodes within each group
  max_year = max(dataset_with_nodes$time)
  if(is.null(window) == F){
    time_windows = seq(min_year, max(dataset_with_nodes$time), window)
  }else if(is.null(N) == F){
    time_windows=seq(min_year, max_year, length.out=N)
  }
  mid_time = time_windows[-length(time_windows)]+(time_windows[2]-time_windows[1])/2
  
  ## Set data frame to store counts
  count_time = matrix(NA, nrow = length(groups), 
                      ncol = length(mid_time))
  
  ## Compute counts through time
  for(i in 1:length(mid_time)){
    t_min = time_windows[i]
    t_max = time_windows[i+1]
    a = which(branches_times[,1] <= t_min & branches_times[,2] >= t_max) ## branches alive, with no sampled individuals 
    b = which(branches_times[,1] <= t_min & (branches_times[,2] >= t_min & branches_times[,2] < t_max)) ## branches born before interval and died within interval
    c = which((branches_times[,1] >= t_min & branches_times[,1] < t_max) & branches_times[,2] >= t_max) ## branches born in interval and died after interval
    d = which((branches_times[,1] >= t_min & branches_times[,1] < t_max) & (branches_times[,2] >= t_min & branches_times[,2] < t_max) ) ## branches born and died within interval
    
    e = c(a,b,c,d)
    e = e[which(branches_group[e,1] - branches_group[e,2] == 0)] ## Filter, to keep only branches that do not have a switch of group
    
    indiv_time = branches_group[e,1] ## Only consider 1 group per branch
    
    tmp = table(factor(indiv_time, levels = as.numeric(groups))) ## Counts per group
    
    count_time[,i] = tmp # Store results
  }
  
  ## Add information on parents groups
  parents=rep(0, length(as.numeric(groups)))
  for (i in as.numeric(groups)){
    beginning_node=dataset_with_nodes$ID[dataset_with_nodes$groups==i & !is.na(dataset_with_nodes$groups)][1]
    path=nodepath(phy=tree,from=beginning_node, to=tree$Nnode+2)
    group_path=dataset_with_nodes$groups[path]
    unique_path=unique(group_path)
    unique_path=unique_path[!is.na(unique_path)]
    if (length(unique_path)>1) parents[i]=unique_path[2]
  }
  # lineage_pres_abs[parents==0,] = 1 # The ancestral lineages are always here
  # parents[which(lineage_pres_abs[,1] == 1)] = 0 ## Simplify stan run and consider that all the lineages that do exist from the strating point are considered as ancestral in the stan code
  
  ## Compute starting times each group
  t_start = rep(NA, length(groups))
  t_start_upper_bound = rep(NA, length(groups))
  t_start_index = rep(NA, length(groups))
  
  for(j in 1:length(t_start)){
    tmp = which(dataset_with_nodes$groups == groups[j])
    tmp2 = dataset_with_nodes$ID[tmp[which.min(dataset_with_nodes$time[tmp])]]
    # if(parents[j] > 0){
    #   m = branches_times[which(tree$edge[,2] == tmp2),1]
    # }else if(parents[j] == 0){
    #   m = branches_times[which(tree$edge[,2] == tmp2),2]
    # }
    m = subset(dataset_with_nodes,ID==tmp2)$time
    if(length(m) == 0) m = min(mid_time) ## Specific case of the root
    t_start[j] = m
    ## Compute index of the starting time (technically not used in the code anymore)
    index = tail(which(mid_time <= t_start[j]), 1)
    if(length(index) == 0) index = 1
    t_start_index[j] = index
  }
  
  ## Compute lineage presence/absence
  lineage_pres_abs =  matrix(0, nrow = length(groups),
                             ncol = length(mid_time))
  for (i in 1:length(groups)){
    lineage_pres_abs[i,t_start_index[i]:dim(lineage_pres_abs)[2]]=1
  }
  parents[which(lineage_pres_abs[,1] == 1)] = 0
  t_start[which(lineage_pres_abs[,1] == 1)] = mid_time[1]-2 ## Make sure the first starting time is before the strat of the time series
  t_start[which.min(t_start)] = mid_time[1]-2
  
  ## Number of ancestral groups that are here from the beginning
  G = length(which(t_start <= mid_time[1] & parents==0))
  
  ## Number of ancestral groups that appear
  GA = length(which(t_start > mid_time[1] & parents==0))
  
  tmp_p <- data.table(group1=1:length(groups),parents1=parents)
  set(tmp_p,which(tmp_p[,parents1!=0]),"t1",tmp_p[which(tmp_p[,parents1!=0]),t_start[parents1]])
  set(tmp_p,which(tmp_p[,parents1==0]),"t1",min_year)
  gap_max = pmin(5,t_start - tmp_p$t1) # 5 is my fixed gap_max
  gap_max = pmax(gap_max,0.0001)
  
  ## Build data object, comprising with all the necessary numbers, vectors and matrices
  data <- list(N=length(mid_time), G=G, GA = GA, K = length(groups), Y = t(count_time),
               parents=parents, t_wind =  window,
               t = mid_time, t_start_approx = t_start, t_start_index = t_start_index, 
               t_new = seq(min_year, max_year, length.out = 150), N_new = 150,
               lin_presence = t(lineage_pres_abs),
               gap_max = gap_max,group_index = which(parents!=0))
  
  initial_values = function(){
    return(list('beta' = rnorm(n = length(groups)-1, mean = 0, sd = 0.05), ## Small betas
                'alpha_true' = rmultinom(n=1, size= 100, prob = rep(1/data$G, data$G))[,1]/100, ## Equal-ish starting frequencies of ancestral groups
                'gamma_true' = abs(rnorm(n=length(groups)-data$G, mean=0, sd=0.001)))) ## Small starting frequencies
  }
  
  fit <- model_compiled$sample(data = data, refresh = 50, #seed=24,
                               chains = 3, parallel_chains = 3,
                               iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                               max_treedepth = 12, adapt_delta = 0.97,
                               init = list(initial_values(), initial_values(), initial_values()))
  
  ## Diagnostic
  fit$cmdstan_diagnose()#check this
  
  ## Extract chains
  t <- list()
  for (f in fit$output_files()) t[[f]] <- as.matrix(data.table::fread(cmd= paste0("grep -v '^#' ", f)))
  t = do.call(rbind, t) ## Combine chains
  t = as.matrix(t)
  Chains = read.chains.from.table(t)
  remove(t)
  
  return(list('fit' = fit,
              'chains' = Chains,
              'data' = data))
}



model_compiled_tstart <- cmdstan_model(stan_file = 'Model_lineage_fitness_tstart.stan')

nodedata <- as.data.table(nodedata)
nodedata <- nodedata[order(nodedata[,ID]),]
nodedata <- as.data.frame(nodedata)
stopifnot(order(nodedata$ID)==nodedata$ID)


res_fitness <- estimate_rel_fitness_groups_with_branches_tstart(
  dataset_with_nodes = nodedata,
  tree = tr1,
  min_year = min_year1,
  window = window1,
  model_compiled = model_compiled_tstart,
  iter_warmup = 1000, iter_sampling = 2000, refresh = 50, seed = 2
)

