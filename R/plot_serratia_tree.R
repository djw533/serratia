


plot_serratia_tree <- function() {


  tree_cols <- c("0" = "black",
                 "1" = "#bab0ac",
                 "2" = "#59a14f",
                 "3" = "#e15759",
                 "4" = "#76b7b2",
                 "5" = "#4e79a7",
                 "6" = "#f28e2b",
                 "7" = "#9c755f",
                 "8" = "#edc948",
                 "9" = "#b07aa1",
                 "10" = "#D7B5A6",
                 "11" = "#ff9da7",
                 "12" = "#a0a0a0",
                 "13" = "#B3E2CD",
                 "14" = "black",
                 "0" = "black")

  names <- c("0" = "Internal node",
             "1" =    "S. quinovorans",
             "2" =    "S. marcescens",
             "3" =    "S. fonticola",
             "4" =    "S. liquefacies",
             "5" =    "S. entomophila",
             "6" =    "S. ficaria",
             "7" =    "S. grimesii",
             "8" =    "S. plymuthica",
             "9" =    "S. proteamaculans",
             "10" =    "S. odorifera",
             "11" =    "S. rubidaea",
             "12" =   "S. odorifera-like",
             "13" =   "S. marcescens-like")


  rooted_serratia_tree <- phytools::midpoint.root(serratia_tree)
  t1 <- ggtree::ggtree(rooted_serratia_tree, ladderize = T, right = T)


  ## fastani clusters:
  ##### dplyr example to get the groups from this dataframe to put into groupOTU:

  dplyr_genus_groups <- fastani_data %>%
    select(Name, X95) %>% # only take the name of the strains and the clusters when cut off at 95% ANI
    mutate(Name = gsub("#","_",Name)) %>% # remove hashes from WSI lane id's
    group_by(X95) %>%
    mutate(strains = list(Name)) %>% # get list of strains in each 95% ANI cluster
    select(X95,strains) %>%
    unique() %>% # remove the strain name/ lane id column and remove duplicate rows
    pull(strains, name = X95) # take the list of the lists of strain names/lane ids and name each item of the list by the fastANI cluster



  genus_groups_tree <- ggtree::groupOTU(rooted_serratia_tree, dplyr_genus_groups) #, overlap='abandon',connect = T)
  #get tibble data:
  tree_tibble <- as_tibble(genus_groups_tree)


  ## plot tree second time with colours:

  t2 <-  ggtree::ggtree(genus_groups_tree, ladderize = T, right = T, size=0.5, aes(color=group)) +
    scale_color_manual(values = c(tree_cols),
                       labels = c(names)) +
    ggtree::geom_treescale(x=0, y=0, offset = 10, width = 0.1, linesize = 0.5) +
    theme(legend.position = "none")



  ### add in fastbaps data:
  ## use level three
  fastbaps_l3 <- fastbaps %>%
    select(Level.3) %>%
    rename(cluster = Level.3) %>%
    tibble::rownames_to_column(var = "strain") %>%
    mutate(cluster = as.factor(cluster))

  ##get the groups of the fastbaps:

  fastbaps_groups <- fastbaps_l3 %>%
    group_by(cluster) %>%
    summarise(strains = list(strain)) %>%
    pull(strains, name = cluster)

  #get the clade names:
  clade_labels <- micro.gen.extra::get_clade_nodes(genus_groups_tree, fastbaps_l3)


  #draw the tree:
  new_tree <- micro.gen.extra::add_clades_to_tree(t2,clade_labels,gradient = F)


  return(new_tree)
}
