



plot_marcescens_tree <- function() {


  serratia.marc.tree <- micro.gen.extra::remove_tips_from_tree(tree = serratia_tree,
                                                               tips = serratia.metadata %>%
                                                                 filter(fastani_derived_species != "Serratia marcescens") %>%
                                                                 pull(Genome_id_in_tree)
  )


  #add lineages on to this
  serratia.marc.lineages <- serratia.metadata %>%
    filter(fastani_derived_species == "Serratia marcescens") %>%
    select(Genome_id_in_tree,Lineage) %>%
    rename(cluster = Lineage,
           strain = Genome_id_in_tree) %>%
    mutate(cluster = as.factor(cluster))

  ##get the groups of the fastbaps:
  serratia.marc.lineages.groups <- serratia.marc.lineages %>%
    group_by(cluster) %>%
    summarise(strains = list(strain)) %>%
    pull(strains, name = cluster)

  #get the clade names:
  serratia.marc.tree.clade_labels <- micro.gen.extra::get_clade_nodes(serratia.marc.tree, serratia.marc.lineages)

  serratia.marc.tree.ggtree <- ggtree::ggtree(serratia.marc.tree, ladderize = T, right = T)
  #draw the tree:
  serratia.marc.tree.ggtree.lineages <- micro.gen.extra::add_clades_to_tree(serratia.marc.tree.ggtree,serratia.marc.tree.clade_labels)

  return(serratia.marc.tree.ggtree.lineages)

}
