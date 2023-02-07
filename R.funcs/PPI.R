
# get interaction of a gene list by STRING PPI database
get_interactions_by_PPI <- function(genes, species="human") {
  all.df <- data.frame()
  #
  if (species=="human") {
    species <- "Homo%20sapiens"
  } else if (species=="mouse") {
    species <- "Mus%20musculus"
  }
  for (i in genes) {
    tmpgene <- i
    url <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=",
                 tmpgene,"&species=",species, sep="")
    # 404 error
    options(warn=-1)
    options(show.error.messages= FALSE)
    options(stringsAsFactors = F)
    webDf <- try(read.table(url, header=T))
    # head(webDf)
    if (is.na(webDf)) next
    all.df <- rbind(all.df, webDf)
  }
  all.df
}

# get the candidate regulators in previous PPI dataframe
# stringId_A	STRING identifier (protein A)
# stringId_B	STRING identifier (protein B)
# preferredName_A	common protein name (protein A)
# preferredName_B	common protein name (protein B)
# ncbiTaxonId	NCBI taxon identifier
# score	combined score
# nscore	gene neighborhood score
# fscore	gene fusion score
# pscore	phylogenetic profile score
# ascore	coexpression score
# escore	experimental score
# dscore	database score
# tscore	textmining score
get_regulators_by_PPI <- function(all.df=interaction.df, min.num=5, filter=F) {
  # filter by score
  if (filter==T) {
    all.df <- subset(all.df, escore>0 | dscore>0 | tscore>0)
  }
  #
  sort.regulators <- sort(table(interaction.df$preferredName_B), decreasing = T)
  sort.regulators <- names(sort.regulators)[sort.regulators >= min.num]
  #
  regulators.df <- data.frame()
  for (i in sort.regulators) {
    targets <- paste(subset(all.df, preferredName_B==i)$preferredName_A, collapse = ",")
    regulators.df <- rbind(regulators.df, data.frame(regulator=i, targets=targets))
  }
  regulators.df
}
