# ECHO Native Circadian Ontological Rhythmicity Exporer (ENCORE) 2D version
# By Hannah De los Santos
# Originated on: 1/14/19

# set directory and clear workspace ----

#https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script/35842176#35842176
# set working directory - only works in RStudio (with rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# clear workspace (sorry!)
rm(list=ls())

# data versions ----

vers_encore <- "3.0.4"
vers_string <- "11.0"

# preload ----

# neurospora library
library(AnnotationHub)

# libraries for ontology explorer
library(shiny)
library(topGO)
library(STRINGdb)
library(r2d3)
library(data.table)
library(jsonlite)
library(ggplot2)
library(igraph)

# libraries for ontology enrichment
library(topGO)
library(stringr)
library(mygene)
library(data.table)
library(STRINGdb)

# ontology libraries;
library(AnnotationDbi)
library(org.Ag.eg.db)
library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.EcK12.eg.db)
library(org.Sc.sgd.db)

# increase input size
options(shiny.maxRequestSize=150*1024^2)

options(shiny.port = 9886)

# functions and variables for generating ontology file ----

# DROSOPHILA

library(DBI)

annFun.drosophila <- function (whichOnto, feasibleGenes = NULL, mapping, ID = "entrez") 
{
  tableName <- c("genes", "accessions", "alias", "ensembl", 
                 "gene_info", "gene_info", "unigene")
  keyName <- c("gene_id", "accessions", "alias_symbol", "ensembl_id", 
               "symbol", "gene_name", "unigene_id")
  names(tableName) <- names(keyName) <- c("entrez", "genbank", 
                                          "alias", "ensembl", "symbol", "genename", "unigene")
  mapping <- paste(sub(".db$", "", mapping), ".db", sep = "")
  require(mapping, character.only = TRUE) || stop(paste("package", 
                                                        mapping, "is required", sep = " "))
  mapping <- sub(".db$", "", mapping)
  geneID <- keyName[tolower(ID)]
  .sql <- paste("SELECT DISTINCT ", geneID, ", go_id FROM ", 
                tableName[tolower(ID)], " INNER JOIN ", paste("go", tolower(whichOnto), 
                                                              sep = "_"), " USING(_id)", sep = "")
  retVal <- dbGetQuery(get(paste(mapping, "dbconn", sep = "_"))(), 
                       .sql)
  retVal <- retVal[!retVal$go_id %in% c("GO:0110109","GO:0120176","GO:0120177", "GO:0120170", "GO:0062023"),]
  if (!is.null(feasibleGenes)) 
    retVal <- retVal[retVal[[geneID]] %in% feasibleGenes, 
                     ]
  return(split(retVal[[geneID]], retVal[["go_id"]]))
}  


# prune dag

prune_dag <- function(ex_graph, pvals, ont_pval_cut){
  sg <- ex_graph
  
  sg_levels <- buildLevels(sg)
  
  sg_relationships <- names(sg@edgeData@data)
  sg_relationships <- strsplit(sg_relationships,"|", fixed = T)
  
  # get children
  children <- as.character(lapply(sg_relationships, `[[`, 1))
  parents <- as.character(lapply(sg_relationships, `[[`, 2))
  
  # children parent dataframe
  cp.df <- data.frame(matrix(0,2,length(sg_relationships)))
  cp.df[1,] <- children
  cp.df[2,] <- parents
  rownames(cp.df) <- c("child","parent")
  # sort the dataframe by the children, most to least
  cp.df <- cp.df[,order(as.numeric(mget(as.character(cp.df[1,]), envir=sg_levels[["nodes2level"]])), decreasing = T)]
  
  # parallel array to the cp.df, keeps track of whether children are significant
  child_sig <- rep(F, ncol(cp.df))
  names(child_sig) <- cp.df[1,]
  keep <- rep(F, sg_levels$noOfNodes)
  names(keep) <- sg@nodes
  
  # get significance
  # pvals <- fc_results[["Damped"]][["BP"]][["go_test"]]@score # all pvalues
  pvals <- pvals[sg@nodes] # relevant pvalues
  keep <- pvals < ont_pval_cut # set only the significant pvalues
  child_sig <- pvals[as.character(cp.df[1,])] < .05
  
  # go through
  for (i in 1:ncol(cp.df)){
    if (child_sig[i]){
      keep[cp.df["parent",i]] <- T
      child_sig[cp.df["parent",i]] <- T
    }
  }
  
  prune_sg <- subGraph(names(keep)[keep],sg)
  
  return(prune_sg)
}

# ont_tax: string
# pval_type: string
# pval_cut: number
# note: has side affects, inherits from getting ontology file section
get_fc_results <- function(ont_tax, pval_type, pval_cut, ont_pval_type, ont_pval_cut, gene_focus, ont_group){
  
  incProgress(1/26, detail = paste("Setting up ontology packages. Started on:",Sys.time()))
  
  gene_ids <- map_sub.df$entrezgene
  
  # get mapping for neurospora, if necessary
  if (ont_tax == "5141"){
    # neurospora crassa
    # neurospora library:
    ah <- AnnotationHub()
    nc_name <- which(AnnotationHub::query(ah,"OrgDb")$species == "Neurospora crassa")

    # now we load in our species! nc will serve as our org.XX.eg.db package
    nc <- AnnotationHub::query(ah,"OrgDb")[[names(AnnotationHub::query(ah,"OrgDb"))[nc_name]]]
    
    kt <- "GID" # ncu numbers
    ont.df <- select(nc,
                     keys = keys(nc,keytype = kt),
                     columns = c("GID","GO_ID","GO_TERM_NAME"),
                     keytype = kt)
    
    # convert to a gene2go list
    gene2GO <- rep(list(c()),length(unique(ont.df$GID)))
    names(gene2GO) <- unique(ont.df$GID)
    for (i in 1:nrow(ont.df)){
      gene2GO[[ont.df$GID[i]]] <- c(gene2GO[[ont.df$GID[i]]], ont.df$GO_ID[i]) 
    }
  }
  
  # keep none significant AC category for each ont category
  no_sig <- list("BP"=c(),"CC"=c(),"MF"=c())
  
  fc.cats <- c("Damped", "Forced", "Harmonic", "Overexpressed", "Repressed", "All.Circ.wo.OE.RE", "All.Circ")
  no_sig[!c("BP","CC","MF") %in% ont_group] <- fc.cats
  # now do analysis for each AC coefficient category
  # including all circadian without overexpressed/repressed, and all circadian
  fc_results <- list()
  for (f in fc.cats){
    # now create factor for gene of interest
    if (is_all_range){
      if (f == "All.Circ.wo.OE.RE"){
        gene_list <- factor(as.integer(gene_ids %in% 
                                         gene_ids[map_sub.df$Osc.Type %in% c("Damped", "Forced", "Harmonic") &
                                                    map_sub.df[pval_type] < pval_cut &
                                                    map_sub.df$query %in% gene_focus])) # to pass to topgo
      } else if (f == "All.Circ"){
        gene_list <- factor(as.integer(gene_ids %in% 
                                         gene_ids[map_sub.df[pval_type] < pval_cut &
                                                    map_sub.df$query %in% gene_focus])) # to pass to topgo
      } else {
        gene_list <- factor(as.integer(gene_ids %in% 
                                         gene_ids[map_sub.df$Osc.Type == f & 
                                                    map_sub.df[pval_type] < pval_cut &
                                                    map_sub.df$query %in% gene_focus])) # to pass to topgo
        
      }
    } else {
      if (f == "All.Circ.wo.OE.RE"){
        gene_list <- factor(as.integer(gene_ids %in% 
                                         gene_ids[map_sub.df$Osc.Type %in% c("Damped", "Forced", "Harmonic") &
                                                    map_sub.df[pval_type] < pval_cut &
                                                    map_sub.df$Period >= low_range &
                                                    map_sub.df$Period <= high_range &
                                                    map_sub.df$query %in% gene_focus])) # to pass to topgo
        
      } else if (f == "All.Circ"){
        gene_list <- factor(as.integer(gene_ids %in% 
                                         gene_ids[map_sub.df[pval_type] < pval_cut &
                                                    map_sub.df$Period >= low_range &
                                                    map_sub.df$Period <= high_range &
                                                    map_sub.df$query %in% gene_focus])) # to pass to topgo
      } else {
        gene_list <- factor(as.integer(gene_ids %in% 
                                         gene_ids[map_sub.df$Osc.Type == f & 
                                                    map_sub.df[pval_type] < pval_cut &
                                                    map_sub.df$Period >= low_range &
                                                    map_sub.df$Period <= high_range &
                                                    map_sub.df$query %in% gene_focus])) # to pass to topgo
        
      }
    }
    
    names(gene_list) <- gene_ids
    levels(gene_list) <- c(0,1)
    gene_list <<- gene_list
    # choices: Biological Process, Cellular Component, Molecular Function
    # now need to do it for each choice
    ontologies <- ont_group
    ontology_list <- list()
    for (ont in ontologies){
      print(paste(f,ont))
      
      incProgress(1/26, detail = paste("Calculating enrichments for",f,ont,"ontology. Started on:",Sys.time()))
      if (ont_tax == "5141"){
        # make a topgo data object to find annotations
        go_data <- new("topGOdata",
                       description = paste(f,ont,"GO Data"),
                       ontology = ont, 
                       allGenes = gene_list,
                       nodeSize = 10,
                       annotationFun = annFUN.gene2GO,
                       gene2GO = gene2GO)
      } else if (ont_tax == "4932"){
        # make a topgo data object to find annotations
        go_data <- new("topGOdata",
                       description = paste(f,ont,"GO Data"),
                       ontology = ont, 
                       allGenes = gene_list,
                       nodeSize = 10,
                       annotationFun = annFUN.org,
                       mapping="org.Sc.sgd.db",
                       ID = "entrez")
      } else if (ont_tax == "7227"){
        # make a topgo data object to find annotations
        go_data <- new("topGOdata",
                       description = paste(f,ont,"GO Data"),
                       ontology = ont, 
                       allGenes = gene_list,
                       nodeSize = 10,
                       annotationFun = annFun.drosophila,
                       mapping=paste0("org.",org_info$Short.Organism.Name[org_info$Taxonomy.Number==ont_tax], ".eg.db"),
                       ID = "entrez")
      } else {
        # make a topgo data object to find annotations
        go_data <- new("topGOdata",
                       description = paste(f,ont,"GO Data"),
                       ontology = ont, 
                       allGenes = gene_list,
                       nodeSize = 10,
                       annotationFun = annFUN.org,
                       mapping=paste0("org.",org_info$Short.Organism.Name[org_info$Taxonomy.Number==ont_tax], ".eg.db"),
                       ID = "entrez")
      }
      
      if (length(go_data@graph@nodes) > 0){
        # run Fisher test for enrichment
        go_test <- runTest(go_data, algorithm = "classic", statistic = "fisher")
        
        go_levels <- buildLevels(go_data@graph) # build the dag graph so we can get levels
        
        # adjust for multiple hypothesis testing
        # before adjusting, PANTHER doesn't take into account terms that don't have at least 2 terms in sig
        go_results <- GenTable(go_data, classicFisher = go_test,
                               # orderBy = "classicFisher", 
                               # ranksOf = "classicFisher", 
                               topNodes = sum(score(go_test) <= 1))
        rownames(go_results) <- go_results$GO.ID
        go_results <- go_results[names(go_test@score),]
        # only keep the ones with at least 2 in significant
        has_2 <- go_results$Significant >= 2
        
        
        if (ont_pval_type == "BH.Adj.P.Value"){
          go_test@score[has_2] <- p.adjust(go_test@score[has_2], method = "BH")
        } else if (ont_pval_type == "BY.Adj.P.Value"){
          go_test@score[has_2] <- p.adjust(go_test@score[has_2], method = "BY")
        } # do not adjust if only pvalue is specified
        
        # prune graph based on pvalue
        prune_sg <- prune_dag(go_data@graph, go_test@score, ont_pval_cut)
        
        if (numNodes(prune_sg) > 0){
          # build new levels
          prune_sg_levels <- buildLevels(prune_sg)
          prune_levels <- as.numeric(as.list(prune_sg_levels$nodes2level))
          names(prune_levels) <- names(as.list(prune_sg_levels$nodes2level))
          # sort prune levels by name
          prune_levels <- prune_levels[order(names(prune_levels))]
          
          # generate a table of results for fisher test
          go_results <- GenTable(go_data, classicFisher = go_test,
                                 # orderBy = "classicFisher", 
                                 # ranksOf = "classicFisher", 
                                 topNodes = sum(score(go_test) <= 1))
          # remove less thans
          go_results$classicFisher <- gsub("< ", "", go_results$classicFisher, fixed=TRUE)
          go_results$classicFisher <- as.numeric(go_results$classicFisher)
          
          # prune the graph
          go_data@graph <- prune_sg
          
          # remove results unrelated to the pruned nodes
          go_results <- go_results[go_results$GO.ID %in% prune_sg@nodes,]
          
          # sort by names
          go_results <- go_results[order(go_results$GO.ID),]
          
          # add level and sort by it
          go_results[,"Level"] <- 0
          go_results$Level <- prune_levels
          
          # sort table by level
          go_results <- go_results[order(go_results$Level),]
          
          # get fold enrichment
          go_results[["Fold.Enrichment"]] <- go_results$Significant/go_results$Expected
          
          # figure out whether everyone has a child
          # preallocate keep
          hasChild <- rep(F,nrow(go_results))
          names(hasChild) <- go_results$GO.ID
          
          # get graph
          sg <- go_data@graph
          
          # get parents and children -- parallel arrays
          sg_relationships <- names(sg@edgeData@data)
          sg_relationships <- strsplit(sg_relationships,"|", fixed = T)
          
          children <- as.character(lapply(sg_relationships, `[[`, 1))
          parents <- as.character(lapply(sg_relationships, `[[`, 2))
          
          hasChild <- names(hasChild) %in% parents
          names(hasChild) <- go_results$GO.ID
          
          go_results$hasChild <- hasChild
          
          # aggregate results
          go_list <- list(title = ont, go_data = go_data, 
                          go_test = go_test, go_results = go_results)
          # put in an overall list for each ontology
          ontology_list[[ont]] <- go_list
        } else {
          no_sig[[ont]] <- c(no_sig[[ont]],f)
        }
      } else {
        no_sig[[ont]] <- c(no_sig[[ont]],f)
      }
    }
    
    fc_results[[f]] <- ontology_list
  }
  
  return(list(fc_results,no_sig))
}

org_info <- read.csv("data//org_map_types.csv", stringsAsFactors = F)
org_id_ex <- read.csv("data//org_id_example.csv", stringsAsFactors = F)

# formulate available organisms
org_avail <- as.character(org_info$Taxonomy.Number)
names(org_avail) <- org_info$Scientific.Organism
id_types <- c("Gene Symbol" = "SYMBOL",
                "Entrez ID"="ENTREZID",
                "Ensembl"="ENSEMBL",
                "Ensembl_Protein"="ENSEMBLPROT",
                "Ensembl_Transcript"="ENSEMBLTRANS",
                "Uniprot" = "UNIPROT",
                "Unigene" = "UNIGENE",
                "NCU Number" = "GID",
                "IPI" = "IPI",
                "SGD" = "SGD",
                "Yeast Symbol" = "COMMON",
                "Flybase"="FLYBASE"
                )

# setting up names for instructions
id_to_common <- names(sort(id_types))
names(id_to_common) <- sort(id_types)
avail_id_types <- data.frame(matrix(0,
                                    length(org_avail),
                                    2+length(id_to_common)))
colnames(avail_id_types) <- c("Organism","Scientific.Name",id_to_common)
avail_id_types$Organism <- org_info$Organism
avail_id_types$Scientific.Name <- org_info$Scientific.Organism

avail_id_types[,-c(1,2)] <- org_info[,names(id_to_common)]

# functions and variables for visualization ----

# preallocate links
links <- NULL

# maximum interactions ??
max_interact <- 100

# color bar for replicates
color_bar <- c("Rep. 1"="red","Rep. 2"="blue","Rep. 3"="green",
               "Rep. 4"="yellow","Rep. 5"="purple","Rep. 6"="pink",
               "Rep. 7"="orange","Rep. 8"="magenta","Fit"="black",
               "Original"="grey")

all_ont_parents <- c("BP"="GO:0008150",
                     "CC"="GO:0005575",
                     "MF"="GO:0003674")

# function for json conversion for go.df
data_to_json <- function(data) {
  jsonlite::toJSON(data, dataframe = "rows", auto_unbox = FALSE, rownames = TRUE, digits = NA)
}

# function to get new table for child values
# par === curr
get_child_table <- function(par, serverValues, sig_in){
  fc_map <- list("Damped" = "DM", "Harmonic" = "HA", "Forced" = "FR",
                 "Overexpressed" = "OE", "Repressed" = "RP",
                 "All.Circ.wo.OE.RE" = "ACW", "All.Circ" = "AC")
  # preallocate where we get frequency data
  go.df <- data.frame(matrix(0,0,3+(length(serverValues$fc_groups) * 3)))
  colnames(go.df) <- c("GO_ID","GO_Term", "Tot_Genes", paste0("Tot_",serverValues$fc_groups), 
                       paste0("Pval_",serverValues$fc_groups), paste0("FE_",serverValues$fc_groups))
  all_pvals <- c() # storing pvalues to see if there's anything significant
  for (f in serverValues$fc_groups){
    if (f != ""){
      if (par != ""){
        
        # get graph
        sg <- fc_results[[f]][[serverValues$ont]][["go_data"]]@graph
        
        # get parents and children -- parallel arrays
        sg_relationships <- names(sg@edgeData@data)
        sg_relationships <- strsplit(sg_relationships,"|", fixed = T)
        
        children <- as.character(lapply(sg_relationships, `[[`, 1))
        parents <- as.character(lapply(sg_relationships, `[[`, 2))
        
        # combine into a dataframe
        cp.df <- data.frame(matrix(0,length(sg_relationships),2))
        cp.df[,1] <- children
        cp.df[,2] <- parents
        colnames(cp.df) <- c("child","parent")
        
        # remove rows not related to our parent
        cp.df <- cp.df[cp.df$parent == par,]
        
        # now cycle through and add to the go data frame as normal
        # note that totals are percentages
        
        #for (ont in ont.types){
        # subset the dataframe with the level we care about
        f_results <- fc_results[[f]][[serverValues$ont]]$go_results # for ease
        # we are only looking at the children
        f_results <- f_results[f_results$GO.ID %in% cp.df$child,]
        if (sig_in == "Significant"){
          level.df <- f_results[f_results$classicFisher < user_input_ont["ont_sig_level"],]
          if (nrow(level.df)==0){ # its only children are insignificant
            level.df <- f_results[f_results$classicFisher >= user_input_ont["ont_sig_level"],]
          }
        } else {
          level.df <- f_results[f_results$classicFisher >= user_input_ont["ont_sig_level"],]
          if (nrow(level.df)==0){ # its only children are significant
            level.df <- f_results[f_results$classicFisher < user_input_ont["ont_sig_level"],]
          }
        }
      } else {
        f_results <- fc_results[[f]][[serverValues$ont]]$go_results # for ease
        if (sig_in == "Significant"){
          level.df <- f_results[f_results$Level == 2 & f_results$classicFisher < user_input_ont["ont_sig_level"],]
        } else {
          level.df <- f_results[f_results$Level == 2 & f_results$classicFisher >= user_input_ont["ont_sig_level"],]
        }
      }
      
      if (nrow(level.df) > 0){
        for (i in 1:nrow(level.df)){
          # if GO category is already accounted for
          if (sum(go.df$GO_ID == level.df$GO.ID[i]) > 0){
            # add to total genes and to that forcing coefficient category
            go.df[go.df$GO_ID == level.df$GO.ID[i],"Tot_Genes"] <-
              go.df[go.df$GO_ID == level.df$GO.ID[i],"Tot_Genes"] +
              (level.df[i,"Significant"]/level.df[i,"Annotated"])
            go.df[go.df$GO_ID == level.df$GO.ID[i],paste0("Tot_",f)] <-
              go.df[go.df$GO_ID == level.df$GO.ID[i],paste0("Tot_",f)] +
              (level.df[i,"Significant"]/level.df[i,"Annotated"])
            go.df[go.df$GO_ID == level.df$GO.ID[i],paste0("Pval_",f)] <-
              go.df[go.df$GO_ID == level.df$GO.ID[i],paste0("Pval_",f)] +
              level.df$classicFisher[i]
            go.df[go.df$GO_ID == level.df$GO.ID[i],paste0("FE_",f)] <-
              go.df[go.df$GO_ID == level.df$GO.ID[i],paste0("FE_",f)] +
              level.df$Fold.Enrichment[i]
          } else { # add to the go.df dataframe
            addrow <- c(level.df$GO.ID[i], level.df$Term[i], 
                        as.numeric(level.df$Significant[i]/level.df$Annotated[i]),as.numeric(rep(0,3*length(serverValues$fc_groups))))
            addrow.df <- data.frame(matrix(0,1,length(addrow)))
            addrow.df[1,] <- addrow
            addrow.df[,-c(1,2)] <- as.numeric(addrow.df[,-c(1,2)]) # add numerics
            colnames(addrow.df) <- colnames(go.df)
            go.df <- rbind(go.df, addrow.df)
            go.df[nrow(go.df),paste0("Tot_",f)] <- level.df$Significant[i]/level.df$Annotated[i]
            go.df[nrow(go.df),paste0("Pval_",f)] <- level.df$classicFisher[i]
            go.df[nrow(go.df),paste0("FE_",f)] <- level.df$Fold.Enrichment[i]
            
            all_pvals <- c(all_pvals,level.df$classicFisher[i])
          }
        }
      }
    }
    #}
    
  }
  
  if (nrow(go.df) > 0){
    # already have checked that there is at least one pvalue
    # adjust if we've only found not significant children
    if (sig_in == "Significant"){
      if (all(all_pvals >= user_input_ont["ont_sig_level"])){
        sig_in <- "Not Significant"
      }
    } else {
      if (all(all_pvals < user_input_ont["ont_sig_level"])){
        sig_in <- "Significant"
      }
    }
    
    # sort by total genes, so colors look pretty
    go.df <- go.df[order(go.df$Tot_Genes, decreasing = T),]
    
    # for the stacked bar, we're going to add another column that is the width
    go.df$Perc_Wid <- go.df$Tot_Genes/sum(go.df$Tot_Genes)
    
    fc_abbrev <- rep(as.character(fc_map[serverValues$fc_groups]),
                     each=nrow(go.df)) #OE, RP
    
    serverValues[["go.df"]] <- go.df
    serverValues[["fc_abbrev"]] <- fc_abbrev
    serverValues$sig <- sig_in
    # serverValues[["corr_color_pal"]] <- color_pal[temp_group]
  } else {
    serverValues[["go.df"]] <- go.df
  }
  
  # if there are no children, the pie chart won't change
  return(serverValues)
}

# color palette
# app
color_pal <-  c("Repressed"="#5D92C0","Damped"="#4B4FAE","Harmonic"="#7D4FB2", "Forced"="#B354B7", "Overexpressed"="#B65388",
                "All.Circ.wo.OE.RE" = "#00D055", "All.Circ" = "#00D055")
color_pal_dark <- c("Repressed"="#427FB2","Damped"="#2E349F","Harmonic"="#6532A2", "Forced"="#A338A8", "Overexpressed"="#A73774",
                    "All.Circ.wo.OE.RE" = "#00A044", "All.Circ" = "#00A044")
# campfire
# color_pal <-  c("Repressed"="#7EBBD0","Damped"="#6F81C1","Harmonic"="#856BBD", "Forced"="#B46FC2", "Overexpressed"="#C16F84")
# color_pal_dark <-  c("Repressed"="#5D92C0","Damped"="#4B4FAE","Harmonic"="#7D4FB2", "Forced"="#B354B7", "Overexpressed"="#B65388")
ont_map <- c("BP" = "Biological Process",
             "CC" =  "Cellular Component",
             "MF" = "Molecular Function")

# information about data and websites
all_org <- c("Mus musculus","Neurospora crassa")
all_org_info <- data.frame(matrix(0,length(all_org),3))
colnames(all_org_info) <- c("Organism","R Package","Gene Website Lookup")
rownames(all_org_info) <- all_org
all_org_info[all_org[1],] <- c(all_org[1],"org.Mm.eg.db","https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
all_org_info[all_org[2],] <- c(all_org[2],"INSERT HERE","https://www.ncbi.nlm.nih.gov/gene/?term=")


prev_parents <- c("")
prev_sig <- c("")
curr <- ""
go_name_curr <- ""
path_trace <- c("")

# ui, ontology explorer ----

ui <- navbarPage("ENCORE: ECHO Native Circadian Ontological Rhythmicity Explorer",
  tabPanel("Explore",
    sidebarLayout(
      sidebarPanel(
        tags$head(
          tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: middle; } 
                #inline .form-group { display: table-row;}")
        ),
        
        tags$p(HTML("<b>Note: Instructions can be found in the 'Instructions' tab.</b> For selections, * indicates inputs that will automatically change visualizations upon change. <b>It is also strongly recommended that your computer is plugged in to run this application.</b>")),
        
        tags$p(paste("ENCORE Version:", vers_encore)),
        
        hr(),
        # leave for ease of testing:
        
        # selectInput(inputId = "dat_select",
        #             # label = "Load Dataset Ontology File",
        #             # accept = ".RData"
        #             # ),
        #             label = "Which dataset would you like to work with?",
        #             choices = c("Emily_Proteome_Quantnorm_20_to_28.RData",
        #                         "Emily_Transcriptome_Quantnorm_20_to_28.RData",
        #                         "Neurospora_RNA_20_to_26.RData"#,"v2pt1_mouse_liver_enrichment_ontology.RData"
        #             )),
        
        div(style="display: inline-block;",
          fileInput(inputId = "dat_select",
                    label = "Upload Dataset ENCORE File (.RData):",
                    accept = ".RData"
                    )),
        div(style="display: inline-block; vertical-align: 3.8em; width: 5px;",
          actionButton(inputId = "dat_load",
                       "Load Data*")),
        tags$br(),
        # choose organism = encodes taxonomy id
        div(style="display: inline-block;",
          selectInput(inputId = "org_select",
                      label = "Choose Dataset Organism:",
                      choices = org_avail,
                      width = "auto")),
        div(style="display: inline-block; vertical-align: .95em; width: 5px;",
          actionButton(inputId = "org_load",
                       "Load Organism*")),
        hr(),
        
        selectInput(inputId = "ont",
                    label = "Which type of ontology would you like to explore?",
                    choices = c("Biological Process" = "BP",
                                "Cellular Component" = "CC",
                                "Molecular Function"= "MF")),
        selectInput(inputId = "fc_map",
                    label = "Which group's Ontology Map would you like to explore?",
                    choices = c("Damped", "Harmonic", "Forced", 
                                "Overexpressed", "Repressed",  
                                "All Circadian (without Overexpressed/Repressed)" = "All.Circ.wo.OE.RE",
                                "All Circadian" = "All.Circ")),
        div(style="display:inline-block; width:400px;",
            selectInput(inputId = "ont_in_map",
                        label = "Which ontology term would you like to see the path to in the Ontology Map?",
                        choices = c(""))
        ),
        div(style="display: inline-block; vertical-align: top; width: 5px;",
            actionButton("ont_term_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_ont_term"),
        
        div(style="display:inline-block; width:400px;",
            checkboxGroupInput(inputId = "fc_groups",
                               label = "Which AC categories would you like to see in the Ontology Explorer/Group Comparison?",
                               choices = c("Repressed", "Damped", "Harmonic", "Forced", 
                                           "Overexpressed", 
                                           "All Circadian (without Overexpressed/Repressed)" = "All.Circ.wo.OE.RE",
                                           "All Circadian" = "All.Circ"))
        ),
        div(style="display: inline-block; vertical-align: top; width: 5px;",
            actionButton("fc_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_fc"),
        
        actionButton(inputId = "update_map",
                     "Update Map*"),
        actionButton(inputId = "restart",
                     "Restart Explorer!*"),
        tags$p(),
        numericInput("num_chord", "Max number of protein-protein interactions:*",
                     min = 0, max = 500, step = 1, value = 50
        ),
        checkboxInput(inputId = "togg_path",
                      label = "Show path at top right?*",
                      value = T),
        checkboxInput(inputId = "togg_border",
                      label = "Turn on visualization border?*",
                      value = T),
        checkboxInput(inputId = "togg_label",
                      label = "Turn off visualization labels?*",
                      value = F),
        checkboxInput(inputId = "togg_dark",
                      label = "Turn on dark mode visualizations?*",
                      value = F),
        checkboxInput(inputId = "togg_fit",
                      label = "Show ECHO fitted values in heatmaps?",
                      value = F),
        checkboxInput(inputId = "togg_sig_groups",
                      label = "Only show genes in groups significantly enriched for the GO Term in chord diagrams?",
                      value = T),
        div(class="header", checked=NA,
            tags$b("Visualization window size (in pixels)*:")),
        
        div(style="display: inline-block; width: 70px;",
            textInput("viz_wid", "width:")),
        
        div(style="display: inline-block; width: 70px;",
            textInput("viz_height", "height:")),
        
        tags$div(id = "inline", numericInput(inputId = "font_size",
                                             label = "Font Size?* :",
                                             value = 16,
                                             min = 0,
                                             step = 1,
                                             width = "100px"))
        
        
      , width = 4),
      mainPanel(
        tabsetPanel(
          tabPanel("Ontology Map", 
            uiOutput("ont_map_ui"),
            tags$br(),
            fluidRow(style = "border: 1px #e3e3e3; border-style: solid; border-radius: 10px; background: #f5f5f5; padding: 10px;",
                     uiOutput("ont_map_below")
            )
          ),
          tabPanel("Ontology Explorer", 
            uiOutput("ont_nav_ui"),
            tags$br(),
            fluidRow(style = "border: 1px #e3e3e3; border-style: solid; border-radius: 10px; background: #f5f5f5; padding: 10px;",
                     uiOutput("ont_nav_below")
            )
          ),
          tabPanel("Group Comparison", 
            uiOutput("group_comp_ui"),
            tags$br(),
            fluidRow(style = "border: 1px #e3e3e3; border-style: solid; border-radius: 10px; background: #f5f5f5; padding: 10px;",
                     uiOutput("group_comp_below")
            )
          ),
          tabPanel("Gene/Term Explorer", fluidPage(
            fluidRow(htmlOutput("frame")),
            fluidRow( column(6,plotOutput("expr", height = "500px"),
                            verbatimTextOutput("expr_text")),
                     column(6,dataTableOutput("gene_list"))
            ) #,
            # style="background: rgb(0, 0, 0);"
            
          )),
          # tabPanel("Chord Genes"), # rename this??
          
          tabPanel("Data Information", verbatimTextOutput("user_data_info")),
          
          # instructions, ontology explorer ----
          
          tabPanel("Instructions",
                   tags$div(class="header", checked = NA,
                            list(
                              tags$p(),
                              tags$p("Welcome to the Gene Ontology Explorer! Here you'll find information on what visualizations are available, and how to use navigate through this app! Note: This is still in beta testing. For contact information, see the end of this tab, or the data information tab of the 'Create ENCORE File' section."),
                              HTML('<center><h2>Navigation</h2></center>'),
                              HTML('<center>'),tags$b("Quick Start:"),HTML('</center>'),
                              tags$p("Navigation: To begin, load data, then select desired ontologies, ontological terms, and categories. To see the path, click on 'Update Map' in the 'Ontology Map' tab. To jump to a specific ontological term's children in the 'Ontology Explorer', click on the term's bar in the map."), 
                              
                              tags$p("In the explorer, to see a ontological term's children and the protein-protein interactions for that ontological term, click on warm-colored ontology rectangle below each bar. To go back, click on the arrow to the left. To switch between significant terms and nonsignificant terms, click the star to the left. More interaction information appears below"),
                              HTML('<center>'),tags$b("Initial Navigation:"),HTML('</center>'),
                              tags$p(HTML("After uploading results from an ENCORE file derived from ECHO results, one begins by selecting their organism, choosing which ontology to explore, which processes to display, and which amplitude change (AC) coefficient categories to look at. All images created using output from <a href='https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000442' target='_blank'>Hughes et al. (2009)</a>")),
                              
                              HTML('<center><img src="ont_map.PNG" style="width:400px"></center><br>'),
                              HTML('<center>'),tags$b("Ontology Map"),HTML('</center>'),
                              tags$p("This is a sankey, or flow, diagram showing the path to a selected ontological term for a specific AC category. Possible paths from the ontological type to the selected term flow from left to right. Not significant terms are grey and significant terms are colored. To jump to a specific ontological category's children, click its bar node; this will update the 'Ontology Explorer' and 'Group Comparison' tabs accordingly, giving you a place to start! Note that this clicked path will be the shortest path available to that term."),
                              tags$p("Further, we can find more information through different interactions:"),
                              
                              HTML('<center><img src="hover_ont_map_node.png" style="width:200px">   <img src="hover_ont_map_link.png" style="width:200px"></center><br>'),
                              tags$p("Hover on a bar corresponding to an ontological term to see its title. Hover on a link to see the source ontological term to the target ontological term. Click on the GO Term for more informaton about that term and the selected AC Category's enrichment below the visualization, and a link for more information about the GO Term appears in Gene/Term Explorer Tab."),
                              
                              HTML('<center><img src="ont_nav.PNG" style="width:400px"></center><br>'),
                              HTML('<center>'),tags$b("Ontology Navigator"),HTML('</center>'),
                              tags$p("This is a layered bar graph representing the fraction of amplitude change coefficient categories that contribute to each ontological catgory, with the width representing the total fraction annotated. To explore the children of a certain category, click on the corresponding colored bar below each bar. For example, clicking on the 'macromolecule localization' bar shows its children, 'RNA localization', 'cellular macromolecule localization', and 'protein localization'."),
                              tags$p("The path taken and dataset appears in the top right corner, with the current parent category selected in bold. To go back, click on the arrow on the left. To switch between looking at significant and nonsignificant categories, click on the star on the left. If nothing changes, this means that categories of the opposite significance do not exist."),
                              tags$p("Further, we can find more information through different interactions:"),
                              
                              HTML('<center><img src="hover_fc_bar.png" style="width:200px">   <img src="hover_ont_bar.png" style="width:200px"></center><br>'),
                              tags$p("Hover on each bar layer for fraction annotated by each category, p-value of enrichment, and fold enrichment information. Hover on ontology bar for total fraction annotated. Click on the GO Term for more informaton about that term below the GO term, and a link for more information appears in Gene/Term Explorer Tab. Click on an AC group's bar text to get enrichment information for that term for the selected AC category below the visualization."),
                              
                              HTML('<center><img src="group_comp.PNG" style="width:400px"></center><br>'),
                              HTML('<center>'),tags$b("Group Comparison:"),HTML('</center>'),
                              tags$p("Once you've selected a GO Term in the Ontology Explorer, you can explore protein-protein ineractions between different genes belonging to that category, as well as their mean-centered, normalized heat maps of gene expression sorted by phase. In the chord diagram, each chord indicates a protein-protein interaction between two genes. Darker chord colors indicate connections within category. Chords are sorted by highest strength to lowest strength connections, and more connections, if any, can be found by increasing the maximum number of protein-protein interactions on the sidebar. Note: if only one gene corresponds to the group, a chord diagram will not be drawn. Also, some ids may not have any STRING mapping and will be removed. Information about genes without identifiers, genes removed for not having a high enough score or connections, and the minimum combined STRING score for the represented chords appears below the visualization."),
                              tags$p("We can also find more information about these connections through different interactions:"),
                              
                              HTML('<center><img src="hover_fc_group.png" style="width:200px">   <img src="hover_gene_group.png" style="width:200px">   <img src="hover_gene_connect.png" style="width:200px">   <img src="hover_hs.png" style="width:200px"></center><br>'),
                              tags$p("Hover on outermost ring to highlight only one AC coefficient category. Hover on any inner arc to highlight a specific gene and find the number of connections on that gene. Hover on chord to highlight the connection to that gene. Hover on the heat map to see the phase shift, in hours, for that gene. To find out more specific gene information, click on gene arc or chord to find out information about that gene or its connecting gene in the Gene/Term Explorer tab. Information about websites used appears at the end of instructions."),
                              
                              HTML('<center><img src="gene_expr.PNG" style="width:400px"></center><br>'),
                              HTML('<center>'),tags$b("Gene Explorer:"),HTML('</center>'),
                              tags$p("When you click on an ontological term in the Ontology Explorer or a gene in the Group Comparison, information about the term appears here. If you click on an ontological term, a link to information about that term will appear above. If you click on a gene from the Group Comparison Tool, the expression with summary information will appear on the left, with a link connecting to gene information on above. Genes appearing in the chord diagram also appear here, along with their amplitude change coefficient categories and hours shifted (phase) values."),
                              
                              # HTML('<center>'),tags$b("Chord Genes:"),HTML('</center>'),
                              # tags$p("Once an ontological term is selected to generate a diagram in the Group Comparison Tool, genes in the chord diagram can be found here, along with their amplitude change coefficient categories and hours shifted (phase) values."),
                              
                              HTML('<center><h2>Information Used<br></h2></center>'),
                              HTML(paste0("For ontologies, this <a href='https://www.ebi.ac.uk/QuickGO/' target='_blank'>website</a>  is used to look up GO terms. Information on protein-protein information is sourced from STRINGdb v", vers_string," for each organism. For specific gene information, <a href='https://www.uniprot.org/' target='_blank'>UniProtKB</a> is used. More version information can be found in the 'Data Version Information' tab in the Create ENCORE file section.")),
                              
                              HTML('<center><h2>Contact and Version Information</h2></center>'),
                              tags$p(HTML("If you are using the results from this app or want to learn about its methods, please <a href='https://github.com/delosh653/ENCORE' target='_blank'>cite us</a>. Additionally, please cite <a href='https://dl.acm.org/citation.cfm?id=3107420&CFID=826084181&CFTOKEN=52238765' target='_blank'>ECHO</a>.")),
                              ("If you run into any errors, please email delosh@rpi.edu with the following (subject 
                               line: ENCORE Error):"),tags$br(),
                              "- a short desciption of your problem" ,tags$br(),
                              "- ENCORE version number",tags$br(),
                              "- your dataset/file(s)",tags$br(),
                              "- your exact settings for the run (a screenshot will do)",tags$br(),
                              "- your exact error from the console window (a screenshot will do)",tags$br(),tags$br(),
                              HTML("All images created by ENCORE using data from: <a href='https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000442' target='_blank'>Hughes et al. (2009)</a>."),
                              tags$br(),tags$br(),
                              tags$p(paste("ENCORE Version", vers_encore))
                              
                              )
                              )
                   )
        )
      , width = 8)
    )
  ),
  # ui, create ontology file ----
  tabPanel("Create ENCORE File", 
    tabsetPanel(
      tabPanel("Create ENCORE File",
        fluidRow(
          column(width = 4,offset = 4,
            HTML("<center><h2>Create ENCORE File</h2>"),
            HTML("Create your ontology file to run ENCORE here! <b>It is highly recommended you plug in your computer before running.</b> Upload an .RData results file from <a href='https://github.com/delosh653/ECHO' target='_blank'>ECHO</a> to begin.<br>"),
            fileInput(inputId = "dat_ont",
                      label = "Upload ECHO Results File (.RData):",
                      accept = ".RData"
            ),
            textInput(inputId = "save_filen",
                      label = "Enter save file name:"),
            selectInput(inputId = "org_tax", # organism taxonomy number, as character
                        label = "Choose the organism for this dataset: ",
                        choices = org_avail),
            div(style="display: inline-block;",
              selectInput(inputId = "id_type",
                          label = "Choose Gene ID type:",
                          choices = id_types)),
            div(style="display: inline-block; vertical-align:top;  width: 20px;",
                actionButton("id_help", icon("question", lib="font-awesome"))),
            uiOutput("Help_id"),
            
            div(style="display: inline-block;",
              textAreaInput("gene_focus","Enter genes to focus enrichment on (line separated):",width="300px",height = "100px")
            ),
            div(style="display: inline-block; vertical-align:top;  width: 20px;",
                actionButton("focus_help", icon("question", lib="font-awesome"))),
            uiOutput("Help_focus"),
            
            checkboxGroupInput("ont_group", "Choose which ontology types to compute:",
                               choices = c("Biological Process" ="BP",
                                           "Molecular Function" = "MF",
                                           "Cellular Component" = "CC"),
                               selected = c("Biological Process" ="BP",
                                            "Molecular Function" = "MF",
                                            "Cellular Component" = "CC")),
            
            selectInput("pval_cat",
                        "Choose P-Value adjustment to use for ECHO significance:",
                        c("Benjamini-Hochberg"="BH.Adj.P.Value",
                          "Benjamini-Yekutieli"="BY.Adj.P.Value",
                          "None"="P.Value")),
            numericInput(inputId = "sig_level",
                         label = "Enter significance level for ECHO significance:",
                         value = .05,
                         min = 0, max = 1, step = .01),
            
            selectInput("ont_pval_cat",
                        "Choose P-Value adjustment to use for Gene Ontology significance:",
                        c("Benjamini-Hochberg"="BH.Adj.P.Value",
                          "Benjamini-Yekutieli"="BY.Adj.P.Value",
                          "None"="P.Value")),
            numericInput(inputId = "ont_sig_level",
                         label = "Enter significance level for Gene Ontology significance:",
                         value = .05,
                         min = 0, max = 1, step = .01),
            
            div(class="header", checked=NA,
                tags$b("Restrict significant rhythms to be between:")),
            
            div(style="display: inline-block; width: 70px;",
                textInput("low_range", "(lower:)")),
            
            div(style="display: inline-block; width: 70px;",
                textInput("high_range", "(upper:)")),
            
            div(style="display: inline-block; vertical-align:top;  width: 20px;",
                actionButton("limit_help", icon("question", lib="font-awesome"))),
            uiOutput("Help_limit"),
            
            actionButton(inputId = "run_ont", label = "Create File!"),
            
            verbatimTextOutput("download_ont_file"),
            
            actionButton(inputId = "ont_file", label = "Download File"),
            
            verbatimTextOutput("dropped_genes"),
            
            HTML("</center>")
          )
        )
      ), # tabpanel ENCORE file
      # instructions, ENCORE file ----
      tabPanel("ENCORE File Help",
               fluidRow(
                 column(width = 12,
                   HTML("<center><h2>ENCORE Help:</h2></center>"),
                   HTML("In order to create visualizations using ENCORE, we first generate a file that has calculated all the ontology enrichments for speed during visualization. After uploading <a href='https://github.com/delosh653/ECHO' target='_blank'>ECHO</a> output and dataset information, simply click start to generate the ENCORE ontology file. A progress bar in the bottom right corner indicates which group is currently being calculated. After calculations are done, download results (note: depending on the size of the organism's genome, the file may be large).<br><br>"),
                   tags$p("Calculations for ontology enrichments are calculated using the whole genome as the background, as specified by PANTHER. For more information on data accession, see the 'Package and Data Information' tab."),
                   tags$p("Note that only certain id types are available for certain organisms:"),
                   tableOutput("org_id_tbl"),
                   tags$p("Further examples of certain id types for each organism appear below. Note that Symbol and Yeast Symbol can be particularly inconsistent in naming schemes for specific organisms."),
                   tableOutput("org_id_ex_tbl")
                   # HTML("</center>")
               ))),
      tabPanel("Data Version Information",
               fluidRow(
                 column(width = 10, offset = 1,
                   HTML("<h2>Data Versioning:</h2>"),
                   HTML(paste0("Protein link information was downloaded from <a href='https://string-db.org/'>STRING v",vers_string,"</a>, and names were mapped to STRING identifiers using the STRINGdb package listed below. Whole genome backgrounds were sourced from the <a href='http://www.pantherdb.org/' target='_blank'>PANTHER</a> database using v14.0, accessed on 01/25/2019.<br><br>")),
                   ("CRAN Package information:"),tags$br(),
                   ("- shiny v1.2.0"),tags$br(),
                   ("- r2d3 v0.2.3"),tags$br(),
                   ("- data.table v1.11.8"),tags$br(),
                   ("- jsonlite v1.5"),tags$br(),
                   ("- ggplot2 v2.3.1.0"),tags$br(),
                   ("- stringr v1.3.1"),tags$br(),
                   ("- igraph v1.2.2"),tags$br(),
                   tags$p("Bioconductor Package Information:"), tags$br(),
                   ("- topGO v2.32.0"),tags$br(),
                   ("- STRINGdb v1.20.0"),tags$br(),
                   ("- mygene v1.16.2"),tags$br(),
                   ("- AnnotationHub v2.14.2"),tags$br(),
                   ("- AnnotationDbi v1.44.0"),tags$br(),
                   ("- org.Ag.eg.db v3.7.0"),tags$br(),
                   ("- org.Dm.eg.db v3.7.0"),tags$br(),
                   ("- org.Hs.eg.db v3.7.0"),tags$br(),
                   ("- org.Mm.eg.db v3.6.0"),tags$br(),
                   ("- org.EcK12.eg.db v3.7.0"),tags$br(),
                   ("- org.Sc.sgd.db v3.7.0"),tags$br(),
                   
                   HTML('<h2>Contact and Version Information</h2>'),
                   tags$p(HTML("If you are using the results from this app or want to learn about its methods, please <a href='https://github.com/delosh653/ENCORE' target='_blank'>cite us</a>. Additionally, please cite <a href='https://dl.acm.org/citation.cfm?id=3107420&CFID=826084181&CFTOKEN=52238765' target='_blank'>ECHO</a>.")),
                   ("If you run into any errors, please email delosh@rpi.edu with the following (subject 
                               line: ENCORE Error):"),tags$br(),
                   "- a short desciption of your problem" ,tags$br(),
                   "- ENCORE version number",tags$br(),
                   "- your dataset/file(s)",tags$br(),
                   "- your exact settings for the run (a screenshot will do)",tags$br(),
                   "- your exact error from the console window (a screenshot will do)",tags$br(),tags$br(),
                   "All images created by ENCORE using data from:",tags$br(),
                   "INSERT CITATION HERE",
                   tags$br(),tags$br(),
                   tags$p(paste("ENCORE Version", vers_encore))
                   # HTML("</center>")
                 )
               ))
    )
  )
)

# server ----

server <-  function(input, output, session) {
  # preallocation ----
  
  serverValues <- reactiveValues()
  
  # help ----
  
  output$Help_sig=renderUI({ # time inputs help
    if(input$sig_help%%2){
      helpText("Sets the intial threshold of significance for those displayed, with a Benjamini-Hochberg adjusted p-value. P-value cutoff is dictated by intial ontology set-up; default cutoff is .05.")
    }
    else{
      return()
    }
  })
  
  output$Help_fc=renderUI({ # time inputs help
    if(input$fc_help%%2){
      helpText("Currently All Circadian is all categories, and All Circadian (without Overexpressed/ Repressed) is Damped, Forced, and Harmonic. If either of the aforementioned categories is chosen, due to overlap with other AC categories, other AC categories (such as Damped) will be unselected for visualizations. If a category has no significant terms for that ontology category, it will be automatically unselected and will not appear in visualizations. ")
    }
    else{
      return()
    }
  })
  
  output$Help_limit=renderUI({ # upper and lower limits for rhythms help
    if(input$limit_help%%2){
      helpText("Upper and lower limits to restrict significant rhythms, in hours. For example, one could look for hours between 20 and 26. If either space is left blank, all significant periods found by ECHO by the specified cutoff will be considered significant.")
    }
    else{
      return()
    }
  })
  
  output$Help_id=renderUI({ # upper and lower limits for rhythms help
    if(input$id_help%%2){
      helpText("Available ids for each organism are specified in the 'Help' tab.")
    }
    else{
      return()
    }
  })
  
  output$Help_focus=renderUI({ # upper and lower limits for rhythms help
    if(input$focus_help%%2){
      helpText("Enter a subset of of ECHO-significant genes to consider enrichments for, with each gene on a separate line. If nothing is entered, all genes will be considered. These must be entered exactly as in the gene names for the dataset, with not additional whitespace. For example, if genes A, B, and C are considered to be significant by ECHO, and the user enters A and B to focus on, C will be left out of all enrichment groups.")
    }
    else{
      return()
    }
  })
  
  output$Help_ont_term=renderUI({ # upper and lower limits for rhythms help
    if(input$ont_term_help%%2){
      helpText("Choose an ontological term by either scrolling through the list, or by pressing backspace and typing a specific term. Terms with asterisks next to them indicate significance by the ontology cutoff specified when creating the ENCORE file. If NA appears, there are no significant terms for the given AC category and ontology type. Choosing the root node of the ontology (such as 'biological_process') will not result in a map, since there is only 1 node.")
    }
    else{
      return()
    }
  })
  
  # ontology explorer back end ----
  
  observeEvent(input$dat_load,{
    withProgress(message = "Loading Data!", value = 0, {
      incProgress(1/2, detail = paste("Loading selected data. Started on:",Sys.time()))
      load(input$dat_select$datapath, envir = globalenv())
      
      end_num <<- 16
      if (user_input$run_conf){
        end_num <<- 26
      }
      
      map_sub.df <- map_sub.df[!is.na(map_sub.df$query),]
      map_sub.df$Osc.Type <- as.character(map_sub.df$Osc.Type)
      
      # making backwards compatible (<v3)
      
      if (is.null(user_input_ont$ont_group)){
        user_input_ont$ont_group <<- "NOT SAVED BEFORE v3 (likewise for focused genes; both default to all)"
        ont_avail <<- c("Biological Process" = "BP",
                        "Cellular Component" = "CC",
                        "Molecular Function"= "MF")
      } else {
        ont_avail <<- user_input_ont$ont_group
        names(ont_avail) <<- ont_map[ont_avail]
      }
      
      if (length(user_input_ont$gene_focus) < 1){
        user_input_ont$gene_focus <<- map_sub.df$entrezgene
      } 
      
      updateSelectInput(session, "ont",
                        label = "Which type of ontology would you like to explore?",
                        choices = ont_avail)
      
      if (exists("fc_results") & !is.null(links)){
        if (!input$fc_map %in% no_sig[[input$ont]]){
          all_terms <- fc_results[[input$fc_map]][[input$ont]][["go_results"]]$Term
          all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level] <-
            paste(all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level],"*")
          
          updateSelectInput(session, "ont_in_map",
                            label = "Which ontology term would you like to see the path to?",
                            choices = all_terms
          )
        } else {
          updateSelectInput(session, "ont_in_map",
                            label = "Which ontology term would you like to see the path to?",
                            choices = NA
                            )
        }
      }
      
      incProgress(1/2, detail = paste("Finished! Finished on:",Sys.time()))
    })
    
  })
  
  observeEvent(input$org_load,{
    withProgress(message = "Loading Organism!", value = 0, {
      incProgress(1/2, detail = paste("Loading selected organism. Started on:",Sys.time()))
      links <<- fread(paste0("links//",input$org_select,".protein.links.v",vers_string,".txt"))
      
      bg <<- read.csv(paste0("data//",input$org_select,"_background.csv"), stringsAsFactors = F)
      
      updateSelectInput(session, "ont",
                        label = "Which type of ontology would you like to explore?",
                        choices = ont_avail)
      
      if (exists("fc_results") & !is.null(links)){
        if (!input$fc_map %in% no_sig[[input$ont]]){
          all_terms <- fc_results[[input$fc_map]][[input$ont]][["go_results"]]$Term
          all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level] <-
            paste(all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level],"*")
          
          updateSelectInput(session, "ont_in_map",
                            label = "Which ontology term would you like to see the path to?",
                            choices = all_terms
          )
        } else {
          updateSelectInput(session, "ont_in_map",
                            label = "Which ontology term would you like to see the path to?",
                            choices = NA
          )
        }
      }
      
      incProgress(1/2, detail = paste("Finished! Finished on:",Sys.time()))
    })
    
  })
  
  observeEvent(input$ont, {
    if (exists("fc_results") & !is.null(links)){
      if (!input$fc_map %in% no_sig[[input$ont]]){
        all_terms <- fc_results[[input$fc_map]][[input$ont]][["go_results"]]$Term
        all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level] <-
          paste(all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level],"*")
        
        updateSelectInput(session, "ont_in_map",
                          label = "Which ontology term would you like to see the path to?",
                          choices = all_terms
        )
      } else {
        updateSelectInput(session, "ont_in_map",
                          label = "Which ontology term would you like to see the path to?",
                          choices = NA
        )
      }
    }
  })
  
  observeEvent(input$fc_map, {
    if (exists("fc_results") & !is.null(links)){
      if (!input$fc_map %in% no_sig[[input$ont]]){
        all_terms <- fc_results[[input$fc_map]][[input$ont]][["go_results"]]$Term
        all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level] <-
          paste(all_terms[fc_results[[input$fc_map]][[input$ont]][["go_results"]]$classicFisher < user_input_ont$ont_sig_level],"*")
        
        updateSelectInput(session, "ont_in_map",
                          label = "Which ontology term would you like to see the path to?",
                          choices = all_terms
        )
      } else {
        updateSelectInput(session, "ont_in_map",
                          label = "Which ontology term would you like to see the path to?",
                          choices = NA
        )
      }
    }
  })
  
  observeEvent(input$update_map, {
    if (exists("fc_results") & !is.null(links) & input$ont_in_map!= "NA" & !input$ont_in_map %in% c("biological_process", "molecular_function", "cellular_component")){
      
      #get example graph
      ex_graph <- fc_results[[input$fc_map]][[input$ont]][["go_data"]]@graph
      care_group <- fc_results[[input$fc_map]][[input$ont]][["go_results"]]
      child <- unlist(strsplit(input$ont_in_map,split = " *",fixed = T))[1] # child of interest
      child <- care_group[which(child==care_group$Term)[1],"GO.ID"]
      
      # get subgraph
      sg <- inducedGraph(ex_graph,child)
      
      # get first path
      ont_parent <- all_ont_parents[input$ont]
      # shortest_path <- rev(shortest_paths(igraph.from.graphNEL(sg), child, to=ont_parent, weights = NA)$vpath[[1]])
      
      # get sig color
      sig_color <- color_pal[input$fc_map]
      no_sig_color <- "#afb1b5"
      
      # get connections -- parents are second col, children are first
      conn <- as.data.frame(t(sapply(strsplit(names(sg@edgeData@data),"|",fixed=T), function(x){return(x)})),
                            stringsAsFactors = F)
      colnames(conn) <- c("child","parent")
      
      # get their significance
      rownames(care_group) <- care_group$GO.ID
      all_sig <- data.frame("name" = unique(c(conn$parent,conn$child)), 
                            "go_name"=unique(c(conn$parent,conn$child)),
                            "sig" = care_group[unique(c(conn$parent,conn$child)), "classicFisher"],
                            "is_sig" = rep(F, length(unique(c(conn$parent,conn$child)))))
      all_sig$is_sig[all_sig$sig < user_input_ont$ont_sig_level] <- T
      all_sig$name<- care_group[unique(c(conn$parent,conn$child)),"Term"]
      
      # rename them to be their casual names
      conn$parent <- care_group[conn$parent,"Term"]
      conn$child <- care_group[conn$child,"Term"]
      
      
      nodes <- data.frame("name"=unique(c(conn$parent,conn$child)), stringsAsFactors = F)
      # using minus 1 to convert to start counting from 0 system
      links <-  data.frame("source"=(sapply(conn$parent, function(x){which(x==nodes$name)}))-1,
                           "target"=(sapply(conn$child, function(x){which(x==nodes$name)}))-1,
                           stringsAsFactors = F)
      links <- cbind(links, cbind("value"=rep(1, nrow(links))))
      
      # combining these in a list so that we can visualize them
      dat_sankey <- list("nodes"=nodes,"links"=links)
      
      serverValues[["dat_sankey"]] <- dat_sankey
      serverValues[["all_sig"]] <- all_sig
      serverValues[["sig_color"]] <- sig_color
      serverValues[["no_sig_color"]] <- no_sig_color
    }
  })
  
  observeEvent(input$new_path, withProgress(message = "Computing!",detail = paste("Started on:",Sys.time()),value = 0, {
    if(!input$new_path %in% all_ont_parents){
      incProgress(1/3, detail = paste("Getting information for Ontology Explorer. Started on:",Sys.time()))
      
      #get example graph
      ex_graph <- fc_results[[input$fc_map]][[input$ont]][["go_data"]]@graph
      care_group <- fc_results[[input$fc_map]][[input$ont]][["go_results"]]
      child <- input$new_path # child of interest -- it is a go id
      # child <- care_group[which(child==care_group$Term)[1],"GO.ID"]
      
      # get subgraph
      sg <- inducedGraph(ex_graph,child)
      
      # get shortest path
      ont_parent <- all_ont_parents[input$ont]
      shortest_path <- names(rev(shortest_paths(igraph.from.graphNEL(sg), input$new_path, to=ont_parent, weights = NA)$vpath[[1]]))
      
      orig_group <- c("Repressed","Damped", "Harmonic","Forced", "Overexpressed","All.Circ.wo.OE.RE","All.Circ")
      
      for (inputId in names(input)) {
        serverValues[[inputId]] <- input[[inputId]]
      }
      
      if (length(serverValues$fc_groups)==0){
        serverValues$fc_groups = c(input$fc_map)
      } else {
        if (!input$fc_map %in% serverValues$fc_groups){
          serverValues$fc_groups <- c(serverValues$fc_groups, input$fc_map)
        }
        
        if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups ){ # overwrites all other groups
          if (input$fc_map != "All.Circ.wo.OE.RE"){
            serverValues$fc_groups <- serverValues$fc_groups[serverValues$fc_groups != "All.Circ.wo.OE.RE"]
          } else {
            serverValues$fc_groups <- c("All.Circ.wo.OE.RE")
          }
        } else if ("All.Circ" %in% serverValues$fc_groups) { # overwrites all other groups except wo
          if (input$fc_map != "All.Circ"){
            serverValues$fc_groups <- serverValues$fc_groups[serverValues$fc_groups != "All.Circ"]
          } else {
            serverValues$fc_groups <- c("All.Circ")
          }
        }
        
        serverValues$fc_groups <- serverValues$fc_groups[order(match(serverValues$fc_groups,orig_group))]
        
      }
      
      # remove if in the "no significant"
      serverValues$fc_groups <- base::setdiff(serverValues$fc_groups,no_sig[[serverValues$ont]])
      
      sig_vect <- c()
      for (ch in shortest_path[-1]){
        if (care_group$classicFisher[care_group$GO.ID==ch] < user_input_ont$ont_sig_level){
          sig_vect[length(sig_vect)+1] <- "Significant"
        } else {
          sig_vect[length(sig_vect)+1] <- "Not Significant"
        }
      }
      # child is the last one in path
      serverValues[["sig"]] <- sig_vect[length(sig_vect)]
      
      # make stuff for the floor
      # groups of interest, provided by shiny app
      # currently fixed, will be updated
      
      # get the table and such
      #serverValues <- get_child_table("",serverValues, serverValues$sig)
      
      serverValues[["color_pal"]] <- color_pal[serverValues$fc_groups]
      serverValues[["color_pal_dark"]] <- color_pal_dark[serverValues$fc_groups]
      serverValues[["go_name"]] <- serverValues$ont
      
      # update the record global variables -- these are now set based on the shortest path
      prev_parents <<- c("", "",shortest_path[-c(1,length(shortest_path))])
      prev_sig <<- c("Significant",sig_vect)
      curr <<- child
      
      rownames(care_group) <- care_group$GO.ID
      c_g <- care_group[shortest_path[-1],]
      common <- c_g$Term
      common[c_g$classicFisher < user_input_ont$ont_sig_level] <- paste(common[c_g$classicFisher < user_input_ont$ont_sig_level],"*")
      
      if (is.null(input$dat_select$name)){
        good_name <- user_input_ont$orig_filen
      } else {
        good_name <- input$dat_select$name
      }
      
      path_trace <<- c(good_name, unname(c(ont_map[serverValues$ont])), common)
      
      if (serverValues$sig == "Significant"){
        path_trace[2] <<- paste(path_trace[2],"*")
      }
      
      orig_sig <- serverValues$sig
      # serverValues$update_val <- F
      serverValues <- get_child_table(curr,serverValues, serverValues$sig)
      prev_sig[length(prev_sig)] <- serverValues$sig
      
      go_term <- curr
      go_name <- character(0)
      for (f in serverValues$fc_groups){
        if (length(go_name) == 0){
          go_name <- fc_results[[f]][[serverValues$ont]][["go_results"]][fc_results[[f]][[serverValues$ont]][["go_results"]]$GO.ID == go_term,"Term"]
        }
      }
      
      serverValues[["go_name"]] <- go_name
      
      incProgress(1/3, detail = paste("Getting information for Chord Diagram. Started on:",Sys.time()))
      
      if (curr != ""){
        # then we also have to get everything for the chord diagram and such
        go_term <- curr
        go_name <- character(0)
        
        for (f in serverValues$fc_groups){
          if (length(go_name) == 0){
            go_name <- fc_results[[f]][[serverValues$ont]][["go_results"]][fc_results[[f]][[serverValues$ont]][["go_results"]]$GO.ID == go_term,"Term"]
          }
        }
        
        int_genes <- c()
        fc_cat_genes <- c()
        for (f in serverValues$fc_groups){
          go_data <- fc_results[[f]][[serverValues$ont]][["go_data"]]
          go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
          
          # grab the genes connected to the term
          if (go_term %in% go_results$GO.ID) {
            go_genes <- genesInTerm(go_data, go_term)
            int_genes <- c(int_genes,go_genes[[go_term]])
          }
        }
        
        # int_genes <- int_genes[!duplicated(int_genes)]
        
        # string db ---
        
        # grab gene stringdb ids, original names
        
        gene_map.df <- data.frame("int_genes"=int_genes, stringsAsFactors = F)#, "fc_cat_genes" = fc_cat_genes)
        rownames(map_sub.df) <- map_sub.df$entrezgene
        gene_map.df$fc_cat_genes <- map_sub.df[gene_map.df$int_genes, "Osc.Type"]
        gene_map.df$STRING_id <- map_sub.df[gene_map.df$int_genes,"STRING_id"]
        gene_map.df$int_genes_orig <- map_sub.df[gene_map.df$int_genes,"query"]
        gene_map.df$period <- map_sub.df[gene_map.df$int_genes,"Period"]
        gene_map.df$pval <- map_sub.df[gene_map.df$int_genes,user_input_ont$pval_cat]
        # no factors!
        i <- sapply(gene_map.df, is.factor)
        gene_map.df[i] <- lapply(gene_map.df[i], as.character)
        
        # remove any mappings not found
        gene_map.df <- gene_map.df[!is.na(gene_map.df$int_genes),]
        gene_map.df <- gene_map.df[!is.na(gene_map.df$fc_cat_genes),]
        
        # remove not selected groups
        # need to change names for if we're looking at all circadian
        if ("All.Circ" %in% serverValues$fc_groups){
          gene_map.df$fc_cat_genes[gene_map.df$fc_cat_genes %in% c("Overexpressed","Repressed","Harmonic","Damped","Forced")] <- "All.Circ"
        } else if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){
          gene_map.df$fc_cat_genes[gene_map.df$fc_cat_genes %in% c("Harmonic","Damped","Forced")] <- "All.Circ.wo.OE.RE"
        }
        gene_map.df <- gene_map.df[gene_map.df$fc_cat_genes %in% serverValues$fc_groups,]
        # remove not selected periods if not all range
        if (!is_all_range){
          gene_map.df <- gene_map.df[gene_map.df$period <= high_range & gene_map.df$period >= low_range,]
        }
        gene_map.df <- gene_map.df[gene_map.df$pval < user_input_ont$sig_level & gene_map.df$int_genes %in% user_input_ont$gene_focus,]
        
        # if we toggle that we only want to see specified signifiance of group, remove others
        if (input$togg_sig_groups){
          keep <- c()
          for (f in serverValues$fc_groups){
            go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
            pval <- go_results[go_results$GO.ID == go_term,"classicFisher"]
            if (length(pval) > 0){
              if (orig_sig == "Significant" & pval < user_input_ont$ont_sig_level){
                keep <- c(keep,f)
              } else if (orig_sig != "Significant" & pval >= user_input_ont$ont_sig_level){
                keep <- c(keep,f)
              }
            }
          }
          
          gene_map.df <- gene_map.df[gene_map.df$fc_cat_genes %in% keep,]
        }
        
        # KEEP GENES THAT DO NOT HAVE STRING_ID MAPPING
        no_string <- gene_map.df[is.na(gene_map.df$STRING_id), "int_genes_orig"]
        
        gene_map.df <- gene_map.df[!is.na(gene_map.df$STRING_id),]
        gene_map.df <- gene_map.df[order(gene_map.df$fc_cat_genes),]
        gene_map.df <- gene_map.df[!duplicated(gene_map.df$STRING_id),]
        
        # get connections - hmmmm
        # these are connections that only appear in the gene ontology
        temp <- links[links$protein1 %in% gene_map.df$STRING_id & links$protein2 %in% gene_map.df$STRING_id,]
        # sort by highest to lowest combined score
        temp <- temp[order(-combined_score),]
        # get rid of reverse duplicates, since doing birectional
        if (nrow(temp) > 0){
          temp <- temp[!duplicated(t(apply(temp[,1:2], 1, sort))), ]
        }
        
        # at the moment, the threshold is 0
        thresh <- as.numeric(temp[nrow(temp),"combined_score"])
        num_low <- 0
        if (nrow(temp) > as.numeric(serverValues$num_chord)){
          # count the number of links that were kicked out for not being above the threshold
          num_low <- nrow(temp) - as.numeric(serverValues$num_chord)
          thresh <- as.numeric(temp[as.numeric(serverValues$num_chord),"combined_score"])
          
          temp <- temp[1:as.numeric(serverValues$num_chord),]
        }
        interacts <- as.data.frame(temp)
        colnames(interacts)[1:2] <- c("from","to") # the square matrix for chords
        #interacts <- interacts[!is.na(interacts$from) & !is.na(interacts$to),]
        
        # map interacts back
        # remove duplicates
        rownames(gene_map.df) <- gene_map.df$STRING_id
        interacts$from <- as.character(gene_map.df[interacts$from,"int_genes_orig"])
        interacts$to <- as.character(gene_map.df[interacts$to,"int_genes_orig"])
        
        # now add the type of gene, based on the from
        rownames(gene_map.df) <- gene_map.df$int_genes_orig
        interacts$fc_type <- gene_map.df[interacts$from, "fc_cat_genes"]
        # sort interacts by type
        interacts <- interacts[order(interacts$fc_type),]
        
        # now make a matrix
        connect.df <- data.frame(matrix(0,length(gene_map.df$int_genes_orig),length(gene_map.df$int_genes_orig)))
        colnames(connect.df) <- rownames(connect.df) <- gene_map.df$int_genes_orig
        
        if (nrow(interacts) > 0){
          for (i in 1:nrow(interacts)){
            connect.df[interacts$from[i] , interacts$to[i]] <- 1
          }
          # remove rows and columns with no connections
          # keep their names
          connect.df <- connect.df[rowSums(connect.df) > 0 | colSums(connect.df) > 0,
                                   rowSums(connect.df) > 0 | colSums(connect.df) > 0]
        }
        # make a symmetric matrix
        connect.df <- (t(connect.df)>0 & connect.df==0)+connect.df
        
        # now get the ones that didn't make the cut
        kicked_out <- gene_map.df$int_genes_orig[!gene_map.df$int_genes_orig %in% rownames(connect.df)]
        
        # heat maps ---
        
        int_genes_care <- rownames(connect.df)
        
        int_tr <- total_results[total_results$`Gene Name` %in% int_genes_care,]
        # need to change names for if we're looking at all circadian
        if ("All.Circ" %in% serverValues$fc_groups){
          int_tr$`Oscillation Type`[int_tr$`Oscillation Type` %in% c("Overexpressed","Repressed","Harmonic","Damped","Forced")] <- "All.Circ"
        } else if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){
          int_tr$`Oscillation Type`[int_tr$`Oscillation Type` %in% c("Harmonic","Damped","Forced")] <- "All.Circ.wo.OE.RE"
        }
        
        hm_total <- matrix(0,nrow(int_tr),length((end_num+1):(end_num+length(timen)*1)))
        hm_order <- rep(0,nrow(hm_total))
        hm_names <- rep("",nrow(hm_total))
        hm_hours_shifted <- rep(0,nrow(hm_total))
        hm_fc <- rep("",nrow(hm_total))
        count_next <- 0
        
        for (f in serverValues$fc_groups){
          int_sub_tr <- int_tr[int_tr$`Oscillation Type` == f & !is.na(int_tr$`Oscillation Type`),]
          if (nrow(int_sub_tr) > 0){
            # there should be no na rows in this data, by default
            # adjust phase
            int_sub_tr$`Phase Shift`[int_sub_tr$Initial.Amplitude < 0] <- int_sub_tr$`Phase Shift`[int_sub_tr$Initial.Amplitude < 0]+pi
            int_sub_tr$Initial.Amplitude[int_sub_tr$Initial.Amplitude < 0] <- -1*int_sub_tr$Initial.Amplitude[int_sub_tr$Initial.Amplitude < 0]
            
            # fixing the phase shift
            int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` > 2*pi] <- int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` > 2*pi]-2*pi
            int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` <0] <- int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` < 0]+2*pi
            
            # HEAT MAP STUFF
            if (!input$togg_fit){
              #get matrix of just the relative expression over time
              hm_mat <- as.matrix(int_sub_tr[,(end_num+1):(end_num+length(timen)*num_reps)])
              
              #if there are replicates, average the relative expression for each replicate
              mtx_reps <- list() # to store actual matrix
              mtx_count <- list() # to store how many are NA
              for (i in 1:num_reps){
                mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
                mtx_count[[i]] <- is.na(mtx_reps[[i]])
                mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
              }
              repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
              hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
              for (i in 1:num_reps){
                hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
                repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
              }
              repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
              hm_mat <- hm_mat/repmtx
              
              # center rows around mean
              # vector of row means
              all_row_mean <- rowMeans(hm_mat, na.rm = TRUE)
              hm_mat <- hm_mat - all_row_mean
              
              
            } else {
              # get fitted values for heatmap instead
              hm_mat <- as.matrix(int_sub_tr[,(end_num+length(timen)*num_reps+1):ncol(int_sub_tr)])
              
            }
            rownames(hm_mat) <- int_sub_tr[,1]
            #normalize each row to be between -1 and 1
            for (i in 1:nrow(int_sub_tr)){
              gene_max <- max(abs((hm_mat[i,])),na.rm = TRUE)
              hm_mat[i,] <- hm_mat[i,]/gene_max 
            }
            #sort by phase shift, if more than 1 gene
            if (nrow(hm_mat) > 1){
              ord <- order(int_sub_tr$`Phase Shift`)
              hm_mat <- hm_mat[ord,]
            } else {
              ord <- 1
            }
            hm_total[(count_next+1):(count_next+nrow(hm_mat)),] <- hm_mat
            hm_order[(count_next+1):(count_next+nrow(hm_mat))] <- which(int_tr$`Oscillation Type` == f)[ord]
            hm_names[(count_next+1):(count_next+nrow(hm_mat))] <- int_sub_tr$`Gene Name`[ord]
            hm_hours_shifted[(count_next+1):(count_next+nrow(hm_mat))] <- int_sub_tr$`Hours Shifted`[ord]
            hm_fc[(count_next+1):(count_next+nrow(hm_mat))] <- f
            count_next <- count_next+nrow(hm_mat)
          }
        }
        
        # making the heat map a data frame
        heat.df <- as.data.frame(hm_total)
        colnames(heat.df) <- paste0("TP_",timen)
        # if (nrow(heat.df) > 0){
        #   heat.df <- heat.df[seq(nrow(heat.df),1,-1),]
        # }
        
        # order the connections to the heatmap
        if (nrow(connect.df)>1){
          connect.df <- connect.df[hm_names,hm_names]
        }
        
        # lowest you can have is one connection
        all_width <- colSums(connect.df)
        all_width[all_width==0] <- .1
        tot_width <- sum(all_width)
        heat.df$Perc_Wid <- all_width/tot_width
        
        # push everything i need to the serverValues
        # stacked_heatmap.js
        serverValues[["time_points"]] <- paste0("TP_",timen)
        serverValues[["genes"]] <- rownames(connect.df)
        serverValues[["heights"]] <- 0:length(timen)
        serverValues[["heat.df"]] <- heat.df
        serverValues[["hm_hours_shifted"]] <- hm_hours_shifted
        
        # circ.db.js
        serverValues[["from_fc_cat"]] <- gene_map.df$fc_cat_genes
        serverValues[["from_genes"]] <- gene_map.df$int_genes_orig
        
        temp_tab <-data.frame(table(gene_map.df$fc_cat_genes[gene_map.df$int_genes_orig %in% rownames(connect.df)]))
        if (nrow(temp_tab) > 0){
          ind <- match(serverValues$fc_groups, temp_tab$Var1, nomatch = 0)
          not_ind <- setdiff(c(1:nrow(temp_tab)),ind)
          temp_tab[,] <- temp_tab[c(ind,not_ind),]
        }
        
        serverValues[["tot_fc_cats"]] <- as.numeric(temp_tab$Freq)
        serverValues[["chord_dat"]] <- connect.df
        
        serverValues[["gene.df"]] <- data.frame("Gene.Name"=rownames(connect.df),
                                                "Osc.Type"=hm_fc,
                                                "Hours.Shifted"=hm_hours_shifted)
        
        
        serverValues[["fc_rep"]] <- temp_tab$Var1[temp_tab$Freq!=0]
        serverValues[["color_pal_chord"]] <- color_pal[as.character(temp_tab$Var1[temp_tab$Freq!=0])]
        serverValues[["color_pal_dark_chord"]] <- color_pal_dark[as.character(temp_tab$Var1[temp_tab$Freq!=0])]
        
        if (sum(serverValues$tot_fc_cats) <= 1){
          serverValues[["group_comp_below_text"]] <- 
            HTML(
              paste0(
                "<b>NOT ENOUGH GENES REMAINING TO CREATE GROUP COMPARISON (need at least 2)</b><br><br>",
                "<b>GO Term:</b> ", path_trace[length(path_trace)], "<br>",
                "<b>Genes sig. annotated for GO Term with no STRING identifier: </b>", paste(no_string, collapse = " "), "<br>",
                "<b>Genes removed for having no connections or not having a high enough score (for the selected maximum number of genes): </b>", paste(kicked_out, collapse = " "), "<br>",
                "<b>Minimum combined STRING score of represented connections:</b>"
              )
            )
        } else {
          serverValues[["group_comp_below_text"]] <- 
            HTML(
              paste0(
                "<b>GO Term:</b> ", path_trace[length(path_trace)], "<br>",
                "<b>Genes sig. annotated for GO Term with no STRING identifier: </b>", paste(no_string, collapse = " "), "<br>",
                "<b>Genes removed for having no connections or not having a high enough score (for the selected maximum number of genes): </b>", paste(kicked_out, collapse = " "), "<br>",
                "<b>Minimum combined STRING score of represented connections:</b>", thresh
              )
            )
        }
        
        if (nrow(serverValues$chord_dat) & !grepl("\\(", path_trace[length(path_trace)])){
          path_trace[length(path_trace)] <<- paste0(path_trace[length(path_trace)]," (",nrow(connect.df)," genes, ",sum(serverValues$chord_dat)/2," chords)")
        }
      }
      incProgress(1/3, detail = paste("Finished! Started on:",Sys.time()))
    }
  }))
  
  observeEvent(input$restart,{
    # as.environment(globalenv())
    if (exists("fc_results") & !is.null(links)){
      serverValues[["darken"]] <- FALSE
      serverValues[["darken_except"]] <- 0
      
      for (inputId in names(input)) {
        serverValues[[inputId]] <- input[[inputId]]
      }
      if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){ # overwrites all other groups
        serverValues$fc_groups <- c("All.Circ.wo.OE.RE")
        
        # serverValues$fc_groups <- unique(c("Damped", "Harmonic", "Forced", serverValues$fc_groups))
        # serverValues$fc_groups <- serverValues$fc_groups[serverValues$fc_groups %in% c("Repressed","Damped","Forced","Harmonic","Overexpressed")]
      } else if ("All.Circ" %in% serverValues$fc_groups) { # overwrites all other groups except wo
        serverValues$fc_groups <- c("All.Circ")
        
        # serverValues$fc_groups <- c("Repressed", "Damped", "Harmonic", "Forced","Overexpressed")
      }
      
      # remove if in the "no significant"
      serverValues$fc_groups <- base::setdiff(serverValues$fc_groups,no_sig[[serverValues$ont]])
      
      serverValues[["sig"]] <- "Significant"
      
      # make stuff for the floor
      # groups of interest, provided by shiny app
      # currently fixed, will be updated
      if (length(serverValues$fc_groups) > 0){
        # get the table and such
        serverValues <- get_child_table("",serverValues, serverValues$sig)
        
        serverValues[["color_pal"]] <- color_pal[serverValues$fc_groups]
        serverValues[["color_pal_dark"]] <- color_pal_dark[serverValues$fc_groups]
        serverValues[["go_name"]] <- serverValues$ont
        # update the record global variables
        prev_parents <<- c("")
        prev_sig <<- c("")
        curr <<- ""
        prev_sig[1] <<- serverValues$sig
        if (is.null(input$dat_select$name)){
          good_name <- user_input_ont$orig_filen
        } else {
          good_name <- input$dat_select$name
        }
        
        path_trace <<- c(good_name, unname(c(ont_map[serverValues$ont])))
        
        if (serverValues$sig == "Significant"){
          path_trace[length(path_trace)] <<- paste(path_trace[(length(path_trace))],"*")
        }
        
        # we don't render anything for the chord diagram -- you need to have chosen something for that
      }
    }
  })
  
  observeEvent(input$num_chord,withProgress(message = "Computing!",detail = paste("Started on:",Sys.time()),value = 0,{
    if (exists("fc_results") & !is.null(links)){
      incProgress(1/2, detail = paste("Getting information for Chord Diagram. Started on:",Sys.time()))
      
      for (inputId in names(input)[names(input) != "fc_groups"]) {
        serverValues[[inputId]] <- input[[inputId]]
      }
      serverValues$fc_groups <- base::setdiff(serverValues$fc_groups,no_sig[[serverValues$ont]])
      
      if (curr != ""){
        # then we also have to get everything for the chord diagram and such
        go_term <- curr
        go_name <- character(0)
        
        for (f in serverValues$fc_groups){
          if (length(go_name) == 0){
            go_name <- fc_results[[f]][[serverValues$ont]][["go_results"]][fc_results[[f]][[serverValues$ont]][["go_results"]]$GO.ID == go_term,"Term"]
          }
        }
        
        int_genes <- c()
        fc_cat_genes <- c()
        for (f in serverValues$fc_groups){
          go_data <- fc_results[[f]][[serverValues$ont]][["go_data"]]
          go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
          
          # grab the genes connected to the term
          if (go_term %in% go_results$GO.ID) {
            go_genes <- genesInTerm(go_data, go_term)
            int_genes <- c(int_genes,go_genes[[go_term]])
          }
        }
        
        # int_genes <- int_genes[!duplicated(int_genes)]
        
        orig_sig <- serverValues$sig
        
        # string db ---
        
        # grab gene stringdb ids, original names
        
        gene_map.df <- data.frame("int_genes"=int_genes, stringsAsFactors = F)#, "fc_cat_genes" = fc_cat_genes)
        rownames(map_sub.df) <- map_sub.df$entrezgene
        gene_map.df$fc_cat_genes <- map_sub.df[gene_map.df$int_genes, "Osc.Type"]
        gene_map.df$STRING_id <- map_sub.df[gene_map.df$int_genes,"STRING_id"]
        gene_map.df$int_genes_orig <- map_sub.df[gene_map.df$int_genes,"query"]
        gene_map.df$period <- map_sub.df[gene_map.df$int_genes,"Period"]
        gene_map.df$pval <- map_sub.df[gene_map.df$int_genes,user_input_ont$pval_cat]
        # no factors!
        i <- sapply(gene_map.df, is.factor)
        gene_map.df[i] <- lapply(gene_map.df[i], as.character)
        
        # remove any mappings not found
        gene_map.df <- gene_map.df[!is.na(gene_map.df$int_genes),]
        gene_map.df <- gene_map.df[!is.na(gene_map.df$fc_cat_genes),]
        
        # remove not selected groups
        # need to change names for if we're looking at all circadian
        if ("All.Circ" %in% serverValues$fc_groups){
          gene_map.df$fc_cat_genes[gene_map.df$fc_cat_genes %in% c("Overexpressed","Repressed","Harmonic","Damped","Forced")] <- "All.Circ"
        } else if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){
          gene_map.df$fc_cat_genes[gene_map.df$fc_cat_genes %in% c("Harmonic","Damped","Forced")] <- "All.Circ.wo.OE.RE"
        }
        gene_map.df <- gene_map.df[gene_map.df$fc_cat_genes %in% serverValues$fc_groups,]
        # remove not selected periods if not all range
        if (!is_all_range){
          gene_map.df <- gene_map.df[gene_map.df$period <= high_range & gene_map.df$period >= low_range,]
        }
        gene_map.df <- gene_map.df[gene_map.df$pval < user_input_ont$sig_level & gene_map.df$int_genes %in% user_input_ont$gene_focus,]
        
        # if we toggle that we only want to see specified signifiance of group, remove others
        if (input$togg_sig_groups){
          keep <- c()
          for (f in serverValues$fc_groups){
            go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
            pval <- go_results[go_results$GO.ID == go_term,"classicFisher"]
            if (length(pval) > 0){
              if (orig_sig == "Significant" & pval < user_input_ont$ont_sig_level){
                keep <- c(keep,f)
              } else if (orig_sig != "Significant" & pval >= user_input_ont$ont_sig_level){
                keep <- c(keep,f)
              }
            }
          }
          
          gene_map.df <- gene_map.df[gene_map.df$fc_cat_genes %in% keep,]
        }
        
        # KEEP GENES THAT DO NOT HAVE STRING_ID MAPPING
        no_string <- gene_map.df[is.na(gene_map.df$STRING_id), "int_genes_orig"]
        
        gene_map.df <- gene_map.df[!is.na(gene_map.df$STRING_id),]
        gene_map.df <- gene_map.df[order(gene_map.df$fc_cat_genes),]
        gene_map.df <- gene_map.df[!duplicated(gene_map.df$STRING_id),]
        
        # get connections - hmmmm
        # these are connections that only appear in the gene ontology
        temp <- links[links$protein1 %in% gene_map.df$STRING_id & links$protein2 %in% gene_map.df$STRING_id,]
        # sort by highest to lowest combined score
        temp <- temp[order(-combined_score),]
        # get rid of reverse duplicates, since doing birectional
        if (nrow(temp) > 0){
          temp <- temp[!duplicated(t(apply(temp[,1:2], 1, sort))), ]
        }
        
        # at the moment, the threshold is 0
        thresh <- as.numeric(temp[nrow(temp),"combined_score"])
        num_low <- 0
        if (nrow(temp) > as.numeric(serverValues$num_chord)){
          # count the number of links that were kicked out for not being above the threshold
          num_low <- nrow(temp) - as.numeric(serverValues$num_chord)
          thresh <- as.numeric(temp[as.numeric(serverValues$num_chord),"combined_score"])
          
          temp <- temp[1:as.numeric(serverValues$num_chord),]
        }
        interacts <- as.data.frame(temp)
        colnames(interacts)[1:2] <- c("from","to") # the square matrix for chords
        #interacts <- interacts[!is.na(interacts$from) & !is.na(interacts$to),]
        
        # map interacts back
        # remove duplicates
        rownames(gene_map.df) <- gene_map.df$STRING_id
        interacts$from <- as.character(gene_map.df[interacts$from,"int_genes_orig"])
        interacts$to <- as.character(gene_map.df[interacts$to,"int_genes_orig"])
        
        # now add the type of gene, based on the from
        rownames(gene_map.df) <- gene_map.df$int_genes_orig
        interacts$fc_type <- gene_map.df[interacts$from, "fc_cat_genes"]
        # sort interacts by type
        interacts <- interacts[order(interacts$fc_type),]
        
        # now make a matrix
        connect.df <- data.frame(matrix(0,length(gene_map.df$int_genes_orig),length(gene_map.df$int_genes_orig)))
        colnames(connect.df) <- rownames(connect.df) <- gene_map.df$int_genes_orig
        
        if (nrow(interacts) > 0){
          for (i in 1:nrow(interacts)){
            connect.df[interacts$from[i] , interacts$to[i]] <- 1
          }
          # remove rows and columns with no connections
          # keep their names
          connect.df <- connect.df[rowSums(connect.df) > 0 | colSums(connect.df) > 0,
                                   rowSums(connect.df) > 0 | colSums(connect.df) > 0]
        }
        # make a symmetric matrix
        connect.df <- (t(connect.df)>0 & connect.df==0)+connect.df
        
        # now get the ones that didn't make the cut
        kicked_out <- gene_map.df$int_genes_orig[!gene_map.df$int_genes_orig %in% rownames(connect.df)]
        
        # heat maps ---
        
        int_genes_care <- rownames(connect.df)
        
        int_tr <- total_results[total_results$`Gene Name` %in% int_genes_care,]
        # need to change names for if we're looking at all circadian
        if ("All.Circ" %in% serverValues$fc_groups){
          int_tr$`Oscillation Type`[int_tr$`Oscillation Type` %in% c("Overexpressed","Repressed","Harmonic","Damped","Forced")] <- "All.Circ"
        } else if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){
          int_tr$`Oscillation Type`[int_tr$`Oscillation Type` %in% c("Harmonic","Damped","Forced")] <- "All.Circ.wo.OE.RE"
        }
        
        hm_total <- matrix(0,nrow(int_tr),length((end_num+1):(end_num+length(timen)*1)))
        hm_order <- rep(0,nrow(hm_total))
        hm_names <- rep("",nrow(hm_total))
        hm_hours_shifted <- rep(0,nrow(hm_total))
        hm_fc <- rep("",nrow(hm_total))
        count_next <- 0
        
        for (f in serverValues$fc_groups){
          int_sub_tr <- int_tr[int_tr$`Oscillation Type` == f & !is.na(int_tr$`Oscillation Type`),]
          if (nrow(int_sub_tr) > 0){
            # there should be no na rows in this data, by default
            # adjust phase
            int_sub_tr$`Phase Shift`[int_sub_tr$Initial.Amplitude < 0] <- int_sub_tr$`Phase Shift`[int_sub_tr$Initial.Amplitude < 0]+pi
            int_sub_tr$Initial.Amplitude[int_sub_tr$Initial.Amplitude < 0] <- -1*int_sub_tr$Initial.Amplitude[int_sub_tr$Initial.Amplitude < 0]
            
            # fixing the phase shift
            int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` > 2*pi] <- int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` > 2*pi]-2*pi
            int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` <0] <- int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` < 0]+2*pi
            
            # HEAT MAP STUFF
            if (!input$togg_fit){
              #get matrix of just the relative expression over time
              hm_mat <- as.matrix(int_sub_tr[,(end_num+1):(end_num+length(timen)*num_reps)])
              
              #if there are replicates, average the relative expression for each replicate
              mtx_reps <- list() # to store actual matrix
              mtx_count <- list() # to store how many are NA
              for (i in 1:num_reps){
                mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
                mtx_count[[i]] <- is.na(mtx_reps[[i]])
                mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
              }
              repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
              hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
              for (i in 1:num_reps){
                hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
                repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
              }
              repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
              hm_mat <- hm_mat/repmtx
              
              # center rows around mean
              # vector of row means
              all_row_mean <- rowMeans(hm_mat, na.rm = TRUE)
              hm_mat <- hm_mat - all_row_mean
              
              
            } else {
              # get fitted values for heatmap instead
              hm_mat <- as.matrix(int_sub_tr[,(end_num+length(timen)*num_reps+1):ncol(int_sub_tr)])
              
            }
            
            #normalize each row to be between -1 and 1
            for (i in 1:nrow(int_sub_tr)){
              gene_max <- max(abs((hm_mat[i,])),na.rm = TRUE)
              hm_mat[i,] <- hm_mat[i,]/gene_max 
            }
            #sort by phase shift, if more than 1 gene
            if (nrow(hm_mat) > 1){
              ord <- order(int_sub_tr$`Phase Shift`)
              hm_mat <- hm_mat[ord,]
            } else {
              ord <- 1
            }
            
            hm_total[(count_next+1):(count_next+nrow(hm_mat)),] <- hm_mat
            hm_order[(count_next+1):(count_next+nrow(hm_mat))] <- which(int_tr$`Oscillation Type` == f)[ord]
            hm_names[(count_next+1):(count_next+nrow(hm_mat))] <- int_sub_tr$`Gene Name`[ord]
            hm_hours_shifted[(count_next+1):(count_next+nrow(hm_mat))] <- int_sub_tr$`Hours Shifted`[ord]
            hm_fc[(count_next+1):(count_next+nrow(hm_mat))] <- f
            count_next <- count_next+nrow(hm_mat)
          }
        }
        
        # making the heat map a data frame
        heat.df <- as.data.frame(hm_total)
        colnames(heat.df) <- paste0("TP_",timen)
        # if (nrow(heat.df) > 0){
        #   heat.df <- heat.df[seq(nrow(heat.df),1,-1),]
        # }
        
        # order the connections to the heatmap
        if (nrow(connect.df)>1){
          connect.df <- connect.df[hm_names,hm_names]
        }
        
        # lowest you can have is one connection
        all_width <- colSums(connect.df)
        all_width[all_width==0] <- .1
        tot_width <- sum(all_width)
        heat.df$Perc_Wid <- all_width/tot_width
        
        # push everything i need to the serverValues
        # stacked_heatmap.js
        serverValues[["time_points"]] <- paste0("TP_",timen)
        serverValues[["genes"]] <- rownames(connect.df)
        serverValues[["heights"]] <- 0:length(timen)
        serverValues[["heat.df"]] <- heat.df
        serverValues[["hm_hours_shifted"]] <- hm_hours_shifted
        
        # circ.db.js
        serverValues[["from_fc_cat"]] <- gene_map.df$fc_cat_genes
        serverValues[["from_genes"]] <- gene_map.df$int_genes_orig
        
        temp_tab <-data.frame(table(gene_map.df$fc_cat_genes[gene_map.df$int_genes_orig %in% rownames(connect.df)]))
        if (nrow(temp_tab) > 0){
          ind <- match(serverValues$fc_groups, temp_tab$Var1, nomatch = 0)
          not_ind <- setdiff(c(1:nrow(temp_tab)),ind)
          temp_tab[,] <- temp_tab[c(ind,not_ind),]
        }
        
        serverValues[["tot_fc_cats"]] <- as.numeric(temp_tab$Freq)
        serverValues[["chord_dat"]] <- connect.df
        
        serverValues[["gene.df"]] <- data.frame("Gene.Name"=rownames(connect.df),
                                                "Osc.Type"=hm_fc,
                                                "Hours.Shifted"=hm_hours_shifted)
        
        
        serverValues[["fc_rep"]] <- temp_tab$Var1[temp_tab$Freq!=0]
        serverValues[["color_pal_chord"]] <- color_pal[as.character(temp_tab$Var1[temp_tab$Freq!=0])]
        serverValues[["color_pal_dark_chord"]] <- color_pal_dark[as.character(temp_tab$Var1[temp_tab$Freq!=0])]
        
        # remove the old listing of genes and chords
        path_trace[length(path_trace)] <<- strsplit(path_trace[length(path_trace)], split = " \\(")[[1]][1]
        
        if (sum(serverValues$tot_fc_cats) <= 1){
          serverValues[["group_comp_below_text"]] <- 
            HTML(
              paste0(
                "<b>NOT ENOUGH GENES REMAINING TO CREATE GROUP COMPARISON (need at least 2)</b><br><br>",
                "<b>GO Term:</b> ", path_trace[length(path_trace)], "<br>",
                "<b>Genes sig. annotated for GO Term with no STRING identifier: </b>", paste(no_string, collapse = " "), "<br>",
                "<b>Genes removed for having no connections or not having a high enough score (for the selected maximum number of genes): </b>", paste(kicked_out, collapse = " "), "<br>",
                "<b>Minimum combined STRING score of represented connections:</b>"
              )
            )
        } else {
          serverValues[["group_comp_below_text"]] <- 
            HTML(
              paste0(
                "<b>GO Term:</b> ", path_trace[length(path_trace)], "<br>",
                "<b>Genes sig. annotated for GO Term with no STRING identifier: </b>", paste(no_string, collapse = " "), "<br>",
                "<b>Genes removed for having no connections or not having a high enough score (for the selected maximum number of genes): </b>", paste(kicked_out, collapse = " "), "<br>",
                "<b>Minimum combined STRING score of represented connections:</b>", thresh
              )
            )
        }
      }
      
      if (nrow(serverValues$chord_dat) & !grepl("\\(", path_trace[length(path_trace)])){
        path_trace[length(path_trace)] <<- paste0(path_trace[length(path_trace)]," (",nrow(connect.df)," genes, ",sum(serverValues$chord_dat)/2," chords)")
      }
      
      incProgress(1/2, detail = paste("Finished! Started on:",Sys.time()))
    }
  }))
  
  # clicking forward and backwards through the ontology
  observeEvent(input$pie_forward,withProgress(message = "Computing!",detail = paste("Started on:",Sys.time()),value = 0, {
    serverValues$fc_groups <- base::setdiff(serverValues$fc_groups,no_sig[[serverValues$ont]])
    
    if (input$pie_forward != curr){
      incProgress(1/3, detail = paste("Getting information for Ontology Explorer. Started on:",Sys.time()))
      if (input$pie_forward == "B"){ # go back up
        if (length(prev_parents) > 1){
          curr <<- prev_parents[length(prev_parents)]
          prev_parents <<- prev_parents[-length(prev_parents)]
          path_trace <<- path_trace[-length(path_trace)]
          if (length(path_trace)>1){
            path_trace <<- path_trace[-length(path_trace)]
          }
          prev_sig <<- prev_sig[-length(prev_sig)]
          serverValues$sig <- prev_sig[length(prev_sig)]
          
          go_term <- curr
          go_name <- character(0)
          for (f in serverValues$fc_groups){
            if (length(go_name) == 0){
              go_name <- fc_results[[f]][[serverValues$ont]][["go_results"]][fc_results[[f]][[serverValues$ont]][["go_results"]]$GO.ID == go_term,"Term"]
            }
          }
          serverValues[["go_name"]] <- go_name
        }
      } else { # go forward down
        hasChild <- F
        for (f in serverValues$fc_groups){
          go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
          if (sum(go_results$GO.ID == input$pie_forward) > 0 && 
              go_results$hasChild[go_results$GO.ID == input$pie_forward]){
            hasChild <- T
          }
        }
        
        # if (hasChild){
          prev_parents <<- c(prev_parents, curr)
          # if (curr != ""){
          #   path_trace <<- c(path_trace, curr)
          # }
          
          prev_sig <<- c(prev_sig, serverValues$sig)
          curr <<- input$pie_forward
        # }
      }
      
      orig_sig <- serverValues$sig
      # serverValues$update_val <- F
      serverValues <- get_child_table(curr,serverValues, serverValues$sig)
      prev_sig[length(prev_sig)] <- serverValues$sig
      
      if (curr == ""){
        serverValues[["go_name"]] <- serverValues$ont
      } else {
        go_term <- curr
        go_name <- character(0)
        for (f in serverValues$fc_groups){
          if (length(go_name) == 0){
            go_name <- fc_results[[f]][[serverValues$ont]][["go_results"]][fc_results[[f]][[serverValues$ont]][["go_results"]]$GO.ID == go_term,"Term"]
          }
        }
        
        serverValues[["go_name"]] <- go_name
      }
      
      incProgress(1/3, detail = paste("Getting information for Chord Diagram. Started on:",Sys.time()))
      
      if (curr != ""){
        # then we also have to get everything for the chord diagram and such
        go_term <- curr
        go_name <- character(0)
        
        for (f in serverValues$fc_groups){
          if (length(go_name) == 0){
            go_name <- fc_results[[f]][[serverValues$ont]][["go_results"]][fc_results[[f]][[serverValues$ont]][["go_results"]]$GO.ID == go_term,"Term"]
          }
        }
        
        hasChild <- T
        
        if (input$pie_forward == "B"|| hasChild){
          path_trace <<- c(path_trace, go_name)
          
          if (orig_sig == "Significant"){
            path_trace[length(path_trace)] <<- paste(path_trace[(length(path_trace))],"*")
          }
          
        }
        
        int_genes <- c()
        fc_cat_genes <- c()
        for (f in serverValues$fc_groups){
          go_data <- fc_results[[f]][[serverValues$ont]][["go_data"]]
          go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
          
          # grab the genes connected to the term
          if (go_term %in% go_results$GO.ID) {
            go_genes <- genesInTerm(go_data, go_term)
            int_genes <- c(int_genes,go_genes[[go_term]])
          }
        }
        
        # int_genes <- int_genes[!duplicated(int_genes)]
        
        # string db ---
        
        # grab gene stringdb ids, original names
        
        gene_map.df <- data.frame("int_genes"=int_genes, stringsAsFactors = F)#, "fc_cat_genes" = fc_cat_genes)
        rownames(map_sub.df) <- map_sub.df$entrezgene
        gene_map.df$fc_cat_genes <- map_sub.df[gene_map.df$int_genes, "Osc.Type"]
        gene_map.df$STRING_id <- map_sub.df[gene_map.df$int_genes,"STRING_id"]
        gene_map.df$int_genes_orig <- map_sub.df[gene_map.df$int_genes,"query"]
        gene_map.df$period <- map_sub.df[gene_map.df$int_genes,"Period"]
        gene_map.df$pval <- map_sub.df[gene_map.df$int_genes,user_input_ont$pval_cat]
        # no factors!
        i <- sapply(gene_map.df, is.factor)
        gene_map.df[i] <- lapply(gene_map.df[i], as.character)
        
        # remove any mappings not found
        gene_map.df <- gene_map.df[!is.na(gene_map.df$int_genes),]
        gene_map.df <- gene_map.df[!is.na(gene_map.df$fc_cat_genes),]
        
        # remove not selected groups
        # need to change names for if we're looking at all circadian
        if ("All.Circ" %in% serverValues$fc_groups){
          gene_map.df$fc_cat_genes[gene_map.df$fc_cat_genes %in% c("Overexpressed","Repressed","Harmonic","Damped","Forced")] <- "All.Circ"
        } else if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){
          gene_map.df$fc_cat_genes[gene_map.df$fc_cat_genes %in% c("Harmonic","Damped","Forced")] <- "All.Circ.wo.OE.RE"
        }
        gene_map.df <- gene_map.df[gene_map.df$fc_cat_genes %in% serverValues$fc_groups,]
        # remove not selected periods if not all range
        if (!is_all_range){
          gene_map.df <- gene_map.df[gene_map.df$period <= high_range & gene_map.df$period >= low_range,]
        }
        gene_map.df <- gene_map.df[gene_map.df$pval < user_input_ont$sig_level & gene_map.df$int_genes %in% user_input_ont$gene_focus,]
        
        # if we toggle that we only want to see specified signifiance of group, remove others
        if (input$togg_sig_groups){
          keep <- c()
          for (f in serverValues$fc_groups){
            go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
            pval <- go_results[go_results$GO.ID == go_term,"classicFisher"]
            if (length(pval) > 0){
              if (orig_sig == "Significant" & pval < user_input_ont$ont_sig_level){
                keep <- c(keep,f)
              } else if (orig_sig != "Significant" & pval >= user_input_ont$ont_sig_level){
                keep <- c(keep,f)
              }
            }
          }
          
          gene_map.df <- gene_map.df[gene_map.df$fc_cat_genes %in% keep,]
        }
        
        # KEEP GENES THAT DO NOT HAVE STRING_ID MAPPING
        no_string <- gene_map.df[is.na(gene_map.df$STRING_id), "int_genes_orig"]
        
        gene_map.df <- gene_map.df[!is.na(gene_map.df$STRING_id),]
        gene_map.df <- gene_map.df[order(gene_map.df$fc_cat_genes),]
        gene_map.df <- gene_map.df[!duplicated(gene_map.df$STRING_id),]
        
        # get connections - hmmmm
        # these are connections that only appear in the gene ontology
        temp <- links[links$protein1 %in% gene_map.df$STRING_id & links$protein2 %in% gene_map.df$STRING_id,]
        # sort by highest to lowest combined score
        temp <- temp[order(-combined_score),]
        # get rid of reverse duplicates, since doing birectional
        if (nrow(temp) > 0){
          temp <- temp[!duplicated(t(apply(temp[,1:2], 1, sort))), ]
        }
        
        # at the moment, the threshold is 0
        thresh <- as.numeric(temp[nrow(temp),"combined_score"])
        num_low <- 0
        if (nrow(temp) > as.numeric(serverValues$num_chord)){
          # count the number of links that were kicked out for not being above the threshold
          num_low <- nrow(temp) - as.numeric(serverValues$num_chord)
          thresh <- as.numeric(temp[as.numeric(serverValues$num_chord),"combined_score"])
          
          temp <- temp[1:as.numeric(serverValues$num_chord),]
        }
        interacts <- as.data.frame(temp)
        colnames(interacts)[1:2] <- c("from","to") # the square matrix for chords
        #interacts <- interacts[!is.na(interacts$from) & !is.na(interacts$to),]
        
        # map interacts back
        # remove duplicates
        rownames(gene_map.df) <- gene_map.df$STRING_id
        interacts$from <- as.character(gene_map.df[interacts$from,"int_genes_orig"])
        interacts$to <- as.character(gene_map.df[interacts$to,"int_genes_orig"])
        
        # now add the type of gene, based on the from
        rownames(gene_map.df) <- gene_map.df$int_genes_orig
        interacts$fc_type <- gene_map.df[interacts$from, "fc_cat_genes"]
        # sort interacts by type
        interacts <- interacts[order(interacts$fc_type),]
        
        # now make a matrix
        connect.df <- data.frame(matrix(0,length(gene_map.df$int_genes_orig),length(gene_map.df$int_genes_orig)))
        colnames(connect.df) <- rownames(connect.df) <- gene_map.df$int_genes_orig
        
        if (nrow(interacts) > 0){
          for (i in 1:nrow(interacts)){
            connect.df[interacts$from[i] , interacts$to[i]] <- 1
          }
          # remove rows and columns with no connections
          # keep their names
          connect.df <- connect.df[rowSums(connect.df) > 0 | colSums(connect.df) > 0,
                                   rowSums(connect.df) > 0 | colSums(connect.df) > 0]
        }
        # make a symmetric matrix
        connect.df <- (t(connect.df)>0 & connect.df==0)+connect.df
        
        # now get the ones that didn't make the cut
        kicked_out <- gene_map.df$int_genes_orig[!gene_map.df$int_genes_orig %in% rownames(connect.df)]
        
        # heat maps ---
        
        int_genes_care <- rownames(connect.df)
        
        int_tr <- total_results[total_results$`Gene Name` %in% int_genes_care,]
        # need to change names for if we're looking at all circadian
        if ("All.Circ" %in% serverValues$fc_groups){
          int_tr$`Oscillation Type`[int_tr$`Oscillation Type` %in% c("Overexpressed","Repressed","Harmonic","Damped","Forced")] <- "All.Circ"
        } else if ("All.Circ.wo.OE.RE" %in% serverValues$fc_groups){
          int_tr$`Oscillation Type`[int_tr$`Oscillation Type` %in% c("Harmonic","Damped","Forced")] <- "All.Circ.wo.OE.RE"
        }
        
        hm_total <- matrix(0,nrow(int_tr),length((end_num+1):(end_num+length(timen)*1)))
        hm_order <- rep(0,nrow(hm_total))
        hm_names <- rep("",nrow(hm_total))
        hm_hours_shifted <- rep(0,nrow(hm_total))
        hm_fc <- rep("",nrow(hm_total))
        count_next <- 0
        
        for (f in serverValues$fc_groups){
          int_sub_tr <- int_tr[int_tr$`Oscillation Type` == f & !is.na(int_tr$`Oscillation Type`),]
          if (nrow(int_sub_tr) > 0){
            # there should be no na rows in this data, by default
            # adjust phase
            int_sub_tr$`Phase Shift`[int_sub_tr$Initial.Amplitude < 0] <- int_sub_tr$`Phase Shift`[int_sub_tr$Initial.Amplitude < 0]+pi
            int_sub_tr$Initial.Amplitude[int_sub_tr$Initial.Amplitude < 0] <- -1*int_sub_tr$Initial.Amplitude[int_sub_tr$Initial.Amplitude < 0]
            
            # fixing the phase shift
            int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` > 2*pi] <- int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` > 2*pi]-2*pi
            int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` <0] <- int_sub_tr$`Phase Shift`[int_sub_tr$`Phase Shift` < 0]+2*pi
            
            # HEAT MAP STUFF
            if (!input$togg_fit){
              #get matrix of just the relative expression over time
              hm_mat <- as.matrix(int_sub_tr[,(end_num+1):(end_num+length(timen)*num_reps)])
              
              #if there are replicates, average the relative expression for each replicate
              mtx_reps <- list() # to store actual matrix
              mtx_count <- list() # to store how many are NA
              for (i in 1:num_reps){
                mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
                mtx_count[[i]] <- is.na(mtx_reps[[i]])
                mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
              }
              repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
              hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
              for (i in 1:num_reps){
                hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
                repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
              }
              repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
              hm_mat <- hm_mat/repmtx
              
              # center rows around mean
              # vector of row means
              all_row_mean <- rowMeans(hm_mat, na.rm = TRUE)
              hm_mat <- hm_mat - all_row_mean
              
              
            } else {
              # get fitted values for heatmap instead
              hm_mat <- as.matrix(int_sub_tr[,(end_num+length(timen)*num_reps+1):ncol(int_sub_tr)])
              
            }
            
            #normalize each row to be between -1 and 1
            for (i in 1:nrow(int_sub_tr)){
              gene_max <- max(abs((hm_mat[i,])),na.rm = TRUE)
              hm_mat[i,] <- hm_mat[i,]/gene_max 
            }
            #sort by phase shift, if more than 1 gene
            if (nrow(hm_mat) > 1){
              ord <- order(int_sub_tr$`Phase Shift`)
              hm_mat <- hm_mat[ord,]
            } else {
              ord <- 1
            }
            
            hm_total[(count_next+1):(count_next+nrow(hm_mat)),] <- hm_mat
            hm_order[(count_next+1):(count_next+nrow(hm_mat))] <- which(int_tr$`Oscillation Type` == f)[ord]
            hm_names[(count_next+1):(count_next+nrow(hm_mat))] <- int_sub_tr$`Gene Name`[ord]
            hm_hours_shifted[(count_next+1):(count_next+nrow(hm_mat))] <- int_sub_tr$`Hours Shifted`[ord]
            hm_fc[(count_next+1):(count_next+nrow(hm_mat))] <- f
            count_next <- count_next+nrow(hm_mat)
          }
        }
        
        # making the heat map a data frame
        heat.df <- as.data.frame(hm_total)
        colnames(heat.df) <- paste0("TP_",timen)
        # if (nrow(heat.df) > 0){
        #   heat.df <- heat.df[seq(nrow(heat.df),1,-1),]
        # }
        
        # order the connections to the heatmap
        if (nrow(connect.df)>1){
          connect.df <- connect.df[hm_names,hm_names]
        }
        
        # lowest you can have is one connection
        all_width <- colSums(connect.df)
        all_width[all_width==0] <- .1
        tot_width <- sum(all_width)
        heat.df$Perc_Wid <- all_width/tot_width
        
        # push everything i need to the serverValues
        # stacked_heatmap.js
        serverValues[["time_points"]] <- paste0("TP_",timen)
        serverValues[["genes"]] <- rownames(connect.df)
        serverValues[["heights"]] <- 0:length(timen)
        serverValues[["heat.df"]] <- heat.df
        serverValues[["hm_hours_shifted"]] <- hm_hours_shifted
        
        # circ.db.js
        serverValues[["from_fc_cat"]] <- gene_map.df$fc_cat_genes
        serverValues[["from_genes"]] <- gene_map.df$int_genes_orig
        
        temp_tab <-data.frame(table(gene_map.df$fc_cat_genes[gene_map.df$int_genes_orig %in% rownames(connect.df)]))
        if (nrow(temp_tab) > 0){
          ind <- match(serverValues$fc_groups, temp_tab$Var1, nomatch = 0)
          not_ind <- setdiff(c(1:nrow(temp_tab)),ind)
          temp_tab[,] <- temp_tab[c(ind,not_ind),]
        }
        
        serverValues[["tot_fc_cats"]] <- as.numeric(temp_tab$Freq)
        serverValues[["chord_dat"]] <- connect.df
        
        serverValues[["gene.df"]] <- data.frame("Gene.Name"=rownames(connect.df),
                                                "Osc.Type"=hm_fc,
                                                "Hours.Shifted"=hm_hours_shifted)
        
        
        serverValues[["fc_rep"]] <- temp_tab$Var1[temp_tab$Freq!=0]
        serverValues[["color_pal_chord"]] <- color_pal[as.character(temp_tab$Var1[temp_tab$Freq!=0])]
        serverValues[["color_pal_dark_chord"]] <- color_pal_dark[as.character(temp_tab$Var1[temp_tab$Freq!=0])]
        
        if (sum(serverValues$tot_fc_cats) <= 1){
          serverValues[["group_comp_below_text"]] <- 
            HTML(
              paste0(
                "<b>NOT ENOUGH GENES REMAINING TO CREATE GROUP COMPARISON (need at least 2)</b><br><br>",
                "<b>GO Term:</b> ", path_trace[length(path_trace)], "<br>",
                "<b>Genes sig. annotated for GO Term with no STRING identifier: </b>", paste(no_string, collapse = " "), "<br>",
                "<b>Genes removed for having no connections or not having a high enough score (for the selected maximum number of genes): </b>", paste(kicked_out, collapse = " "), "<br>",
                "<b>Minimum combined STRING score of represented connections:</b>"
              )
            )
        } else {
          serverValues[["group_comp_below_text"]] <- 
            HTML(
              paste0(
                "<b>GO Term:</b> ", path_trace[length(path_trace)], "<br>",
                "<b>Genes sig. annotated for GO Term with no STRING identifier: </b>", paste(no_string, collapse = " "), "<br>",
                "<b>Genes removed for having no connections or not having a high enough score (for the selected maximum number of genes): </b>", paste(kicked_out, collapse = " "), "<br>",
                "<b>Minimum combined STRING score of represented connections:</b>", thresh
              )
            )
        }
        
        if (nrow(serverValues$chord_dat) & !grepl("\\(", path_trace[length(path_trace)])){
          path_trace[length(path_trace)] <<- paste0(path_trace[length(path_trace)]," (",nrow(connect.df)," genes, ",sum(serverValues$chord_dat)/2," chords)")
        }
      }
      
      incProgress(1/3, detail = paste("Finished! Started on:",Sys.time()))
    }
  }))
  
  # switching between significant and not significant view
  observeEvent(input$sig_in, {
    serverValues$fc_groups <- base::setdiff(serverValues$fc_groups,no_sig[[serverValues$ont]])
    
    if (input$sig_in != serverValues$sig){
      criteria <- F
      for (f in serverValues$fc_groups){ # go through listed categories
        go_results <- fc_results[[f]][[serverValues$ont]][["go_results"]]
        if (sum(go_results$GO.ID == curr) > 0 | curr == ""){
          criteria <- (curr == "" || go_results$classicFisher[go_results$GO.ID == curr] > user_input_ont["ont_sig_level"] || go_results$hasChild[go_results$GO.ID == curr])
        }
        if(criteria){ break }
      }
      
      if (criteria){
        orig_sig <- serverValues$sig
        # serverValues$sig <- input$sig_in
        # prev_sig[length(prev_sig)] <<- serverValues$sig
        serverValues <- get_child_table(curr,serverValues, input$sig_in)
        
        if (orig_sig != serverValues$sig){
          if (serverValues$sig == "Significant" & !grepl("*",path_trace[length(path_trace)], fixed = T)){
            top_str <- strsplit(path_trace[length(path_trace)],split = " (",fixed = TRUE)
            path_trace[length(path_trace)] <<- paste0(top_str[[1]][1]," *",
                                                      " (",top_str[[1]][2])
          } else {
            path_trace[length(path_trace)] <<- gsub("( \\*)","", path_trace[length(path_trace)])
          }
        }
      }
    }
    
  })
  
  # clicking on gene chord to open/display information
  observeEvent(input$url, {
    # need to map to uniprot
    serverValues[["term_of_interest"]] <- input$url
    gene_care <- bg[which(bg[,user_input_ont$id_type]==input$url),"UNIPROT"]
    
    serverValues[["url"]] <- paste0("https://www.uniprot.org/uniprot/",gene_care)
    
    # also make the ribbon dataframe for the gene expression
    rep_genes <- total_results[total_results$'Gene Name'==input$url,(end_num+1):(end_num+(length(timen)*num_reps))]
    
    ribbon.df <- data.frame(matrix(ncol = 4+num_reps, nrow = length(timen)))
    colnames(ribbon.df) <- c("Times","Fit","Min","Max", paste(rep("Rep",num_reps),c(1:num_reps), sep=".")) # assigning column names
    ribbon.df$Times <- timen
    ribbon.df$Fit <- t(total_results[total_results$'Gene Name'==input$url,c((end_num+1+(length(timen)*num_reps)):ncol(total_results))]) # assigning the fit
    ribbon.df$Min <- sapply(seq(1,ncol(rep_genes), by = num_reps), function(x) min(unlist(rep_genes[,c(x:(num_reps-1+x))]), na.rm = TRUE)) # getting min values of replicates
    ribbon.df$Max <- sapply(seq(1,ncol(rep_genes), by = num_reps), function(x) max(unlist(rep_genes[,c(x:(num_reps-1+x))]), na.rm = TRUE)) # getting max values of replicates
    for (i in 1:num_reps){ # assign each of the replicates
      ribbon.df[,4+i] <- t(rep_genes[,seq(i,ncol(rep_genes),by=num_reps)])
    }
    colnames(ribbon.df) <- c("Times","Fit","Min","Max", paste(rep("Rep",num_reps),c(1:num_reps), sep=".")) # assigning column names
    
    serverValues[["ribbon.df"]] <- ribbon.df 
    serverValues[["tr_sub"]] <-  total_results[total_results$'Gene Name'==input$url,] # for printing relevant parameters
    serverValues[["is_go"]] <- FALSE
    
  })
  
  # observing clicking on a go term on the ontology map
  observeEvent(input$go_url_map,{
    f <- input$fc_map
    go_data <- fc_results[[f]][[input$ont]][["go_data"]]
    care_group <- fc_results[[f]][[input$ont]][["go_results"]]
    go_genes <- genesInTerm(go_data, input$go_url_map)[[input$go_url_map]]
    map_intr <- map_sub.df[map_sub.df$entrezgene %in% go_genes,]
    # subset to only our category and significance
    # have to change f if all circadian w/wo oe/re
    if (f == "All.Circ"){
      f <- c("Overexpressed","Repressed","Harmonic","Forced","Damped")
    } else if (f == "All.Circ.wo.OE.RE"){
      f <- c("Harmonic","Forced","Damped")
    }
    
    if (!is_all_range){
      map_intr_sub <- map_intr[map_intr$Osc.Type %in% f & map_intr[[user_input_ont$pval_cat]] < user_input_ont$sig_level &
                                 map_intr$Period <= high_range & map_intr$Period >= low_range & map_intr$entrezgene %in% user_input_ont$gene_focus,]
    } else {
      map_intr_sub <- map_intr[map_intr$Osc.Type %in% f & map_intr[[user_input_ont$pval_cat]] < user_input_ont$sig_level & 
                                 map_intr$entrezgene %in% user_input_ont$gene_focus,]
    }
    go_genes_sig <- map_intr_sub$query[!is.na(map_intr_sub$query)]
    
    # get the term name
    child <- care_group[which(input$go_url_map==care_group$GO.ID)[1],"Term"]
    
    serverValues[["url"]] <- paste0("https://www.ebi.ac.uk/QuickGO/term/",input$go_url_map)
    serverValues[["term_of_interest"]] <- child 
    serverValues[["is_go"]] <- TRUE
    
    start_string <- paste0("<b>AC Category: </b>", input$fc_map, "<br>")
    if (input$fc_map == "All.Circ"){
      start_string <- paste0(
        "<b>AC Category: </b>", input$fc_map, "<br>",
        "<b>&emsp; - Total no. of sig. genes: </b>", length(go_genes_sig), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Repressed: </b>", sum(map_intr_sub$Osc.Type == "Repressed", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Damped: </b>", sum(map_intr_sub$Osc.Type == "Damped", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Harmonic: </b>", sum(map_intr_sub$Osc.Type == "Harmonic", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Forced: </b>", sum(map_intr_sub$Osc.Type == "Forced", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Overexpressed: </b>", sum(map_intr_sub$Osc.Type == "Overexpressed", na.rm = T), "<br>"
      )
    } else if (input$fc_map == "All.Circ.wo.OE.RE"){
      start_string <- paste0(
        "<b>AC Category: </b>", input$fc_map, "<br>",
        "<b>&emsp; - Total no. of sig. genes: </b>", length(go_genes_sig), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Damped: </b>", sum(map_intr_sub$Osc.Type == "Damped", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Harmonic: </b>", sum(map_intr_sub$Osc.Type == "Harmonic", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Forced: </b>", sum(map_intr_sub$Osc.Type == "Forced", na.rm = T), "<br>"
      )
    }
    
    serverValues[["ont_map_below_text"]] <-
      HTML(
        paste0(
          start_string,
          "<b>GO Term: </b>", serverValues$term_of_interest, "<br>",
          "<b>GO ID: </b>", input$go_url_map, "<br>",
          "<b>GO Term Level: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Level"], "<br>",
          "<b>More information about this GO Term: </b>", "<a href='", serverValues$url,"' target='_blank'>QuickGO</a>", "<br>",
          "<b>Enrichment Information: </b>", "<br>",
          "<b>&emsp; - Fraction Annotated: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Significant"]/care_group[which(input$go_url_map==care_group$GO.ID)[1],"Annotated"], "<br>",
          "<b>&emsp;&emsp;   - Number of sig. genes in AC Category annotated for this GO Term: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Significant"], "<br>",
          "<b>&emsp;&emsp;   - Total number of genes annotated for this GO Term: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Annotated"], "<br>",
          "<b>&emsp; - P-Value: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"classicFisher"], "<br>",
          "<b>&emsp; - P-Value Type: </b>", user_input_ont$ont_pval_cat, "<br>",
          "<b>&emsp; - Fold Enrichment: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Fold.Enrichment"], "<br>",
          "<b>&emsp;&emsp;   - Number of sig. genes in AC Category annotated for the GO Term: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Significant"], "<br>",
          "<b>&emsp;&emsp;   - Number of genes Expected for the GO Term: </b>", care_group[which(input$go_url_map==care_group$GO.ID)[1],"Expected"], "<br>",
          "<b>List of sig. genes in AC Category annotated : </b>", "<br>",
          paste(go_genes_sig, collapse = " ")
          
        )
      )
  })
  
  # observing click on a go term on the ontology navigator
  observeEvent(input$go_url_nav,{
    for (f in serverValues$fc_groups){
      if (!f %in% no_sig[[serverValues$ont]]){
        go_data <- fc_results[[f]][[input$ont]][["go_data"]]
        care_group <- fc_results[[f]][[serverValues$ont]][["go_results"]]
        if (input$go_url_nav %in% care_group$GO.ID){
          go_row <- which(care_group$GO.ID == input$go_url_nav)[1]
          break
        }
      }
    }
    # get the term name
    child <- care_group[go_row,"Term"]
    
    serverValues[["url"]] <- paste0("https://www.ebi.ac.uk/QuickGO/term/",input$go_url_nav)
    serverValues[["term_of_interest"]] <- child 
    serverValues[["is_go"]] <- TRUE
    
    # calculate the total amount of significant genes
    
    serverValues[["ont_nav_below_text"]] <-
      HTML(
        paste0(
          "<b>GO Term: </b>", serverValues$term_of_interest, "<br>",
          "<b>GO ID: </b>", input$go_url_nav, "<br>",
          "<b>GO Term Level: </b>", care_group[go_row,"Level"], "<br>",
          "<b>More information about this GO Term: </b>", "<a href='", serverValues$url,"' target='_blank'>QuickGO</a>", "<br>",
          "<b>Total Fraction Annotated: </b>", serverValues$go.df[serverValues$go.df$GO_ID == input$go_url_nav,"Tot_Genes"], "<br>",
          "<b>&emsp;&emsp;   - Number of sig. genes in all AC Categories annotated for this GO Term: </b>", care_group[go_row,"Annotated"]*serverValues$go.df[serverValues$go.df$GO_ID == input$go_url_nav,"Tot_Genes"], "<br>",
          "<b>&emsp;&emsp;   - Total number of genes annotated for this GO Term: </b>", care_group[go_row,"Annotated"], "<br>"
        )
      )
  })
  
  # observing click on an ac cat bar on the ontology map
  observeEvent(input$fc_nav, {
    # split fc nav into the go id and the category
    f <- unlist(strsplit(input$fc_nav, " "))[1]
    go_id <- unlist(strsplit(input$fc_nav, " "))[2]
    
    go_data <- fc_results[[f]][[input$ont]][["go_data"]]
    go_genes <- genesInTerm(go_data, go_id)[[go_id]]
    map_intr <- map_sub.df[map_sub.df$entrezgene %in% go_genes,]
    # subset to only our category and significance
    f_repl <- f
    if (f == "All.Circ"){
      f_repl <- c("Overexpressed","Repressed","Harmonic","Forced","Damped")
    } else if (f == "All.Circ.wo.OE.RE"){
      f_repl <- c("Harmonic","Forced","Damped")
    }
    
    if (!is_all_range){
      map_intr_sub <- map_intr[map_intr$Osc.Type %in% f_repl & map_intr[[user_input_ont$pval_cat]] < user_input_ont$sig_level &
                                 map_intr$Period <= high_range & map_intr$Period >= low_range & map_intr$entrezgene %in% user_input_ont$gene_focus,]
    } else {
      map_intr_sub <- map_intr[map_intr$Osc.Type %in% f_repl & map_intr[[user_input_ont$pval_cat]] < user_input_ont$sig_level & 
                                 map_intr$entrezgene %in% user_input_ont$gene_focus,]
    }
    go_genes_sig <- map_intr_sub$query[!is.na(map_intr_sub$query)]
    
    # get the term name
    care_group <- fc_results[[f]][[input$ont]][["go_results"]]
    child <- care_group[which(go_id==care_group$GO.ID)[1],"Term"]
    
    serverValues[["url"]] <- paste0("https://www.ebi.ac.uk/QuickGO/term/",go_id)
    serverValues[["term_of_interest"]] <- child 
    serverValues[["is_go"]] <- TRUE
    
    start_string <- paste0("<b>AC Category: </b>", f, "<br>")
    if (f == "All.Circ"){
      start_string <- paste0(
        "<b>AC Category: </b>", f, "<br>",
        "<b>&emsp; - Total no. of sig. genes: </b>", length(go_genes_sig), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Repressed: </b>", sum(map_intr_sub$Osc.Type == "Repressed", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Damped: </b>", sum(map_intr_sub$Osc.Type == "Damped", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Harmonic: </b>", sum(map_intr_sub$Osc.Type == "Harmonic", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Forced: </b>", sum(map_intr_sub$Osc.Type == "Forced", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Overexpressed: </b>", sum(map_intr_sub$Osc.Type == "Overexpressed", na.rm = T), "<br>"
      )
    } else if (f == "All.Circ.wo.OE.RE"){
      start_string <- paste0(
        "<b>AC Category: </b>", f, "<br>",
        "<b>&emsp; - Total no. of sig. genes: </b>", length(go_genes_sig), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Damped: </b>", sum(map_intr_sub$Osc.Type == "Damped", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Harmonic: </b>", sum(map_intr_sub$Osc.Type == "Harmonic", na.rm = T), "<br>",
        "<b>&emsp;&emsp; - No. of sig. genes Forced: </b>", sum(map_intr_sub$Osc.Type == "Forced", na.rm = T), "<br>"
      )
    }
    
    serverValues[["ont_nav_below_text"]] <-
      HTML(
        paste0(
          start_string,
          "<b>GO Term: </b>", serverValues$term_of_interest, "<br>",
          "<b>GO ID: </b>", go_id, "<br>",
          "<b>GO Term Level: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Level"], "<br>",
          "<b>More information about this GO Term: </b>", "<a href='", serverValues$url,"' target='_blank'>QuickGO</a>", "<br>",
          "<b>Enrichment Information: </b>", "<br>",
          "<b>&emsp; - Fraction Annotated: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Significant"]/care_group[which(go_id==care_group$GO.ID)[1],"Annotated"], "<br>",
          "<b>&emsp;&emsp;   - Number of sig. genes in AC Category annotated for this GO Term: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Significant"], "<br>",
          "<b>&emsp;&emsp;   - Total number of genes annotated for this GO Term: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Annotated"], "<br>",
          "<b>&emsp; - P-Value: </b>", care_group[which(go_id==care_group$GO.ID)[1],"classicFisher"], "<br>",
          "<b>&emsp; - P-Value Type: </b>", user_input_ont$ont_pval_cat, "<br>",
          "<b>&emsp; - Fold Enrichment: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Fold.Enrichment"], "<br>",
          "<b>&emsp;&emsp;   - Number of sig. genes in AC Category annotated for the GO Term: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Significant"], "<br>",
          "<b>&emsp;&emsp;   - Number of genes Expected for the GO Term: </b>", care_group[which(go_id==care_group$GO.ID)[1],"Expected"], "<br>",
          "<b>List of sig. genes in AC Category annotated : </b>", "<br>",
          paste(go_genes_sig, collapse = " ")
          
        )
      )
  })
  
  observeEvent(input$togg_path, {
    serverValues[["togg_path"]] <- input$togg_path
  })
  
  observeEvent(input$togg_dark, {
    if (input$togg_dark){
      serverValues[["textcol_bg"]] <- "#FFFFFF"
      # serverValues[["textcol_fore"]] <- "#000000"
      serverValues[["col_bg"]] <- "#000000"
    } else {
      serverValues[["textcol_bg"]] <- "#000000"
      # serverValues[["textcol_fore"]] <- "#000000"
      serverValues[["col_bg"]] <- "#FFFFFF"
    }
  })
  
  # render output ontology explorer ----
  
  output$ont_map_ui <- renderUI({
    if (input$viz_wid == "" | grepl("\\D", input$viz_wid)){
      wid <- "100%"
    } else {
      wid <- paste0(input$viz_wid,"px")
    }
    
    if (input$viz_height == "" | grepl("\\D", input$viz_height)){
      height <- "840px"
    } else {
      height <- paste0(input$viz_height,"px")
    }
    d3Output("ont_map",width = wid, height = height)
  })
  
  output$ont_nav_ui <- renderUI({
    if (input$viz_wid == "" | grepl("\\D", input$viz_wid)){
      wid <- "100%"
    } else {
      wid <- paste0(input$viz_wid,"px")
    }
    
    if (input$viz_height == "" | grepl("\\D", input$viz_height)){
      height <- "840px"
    } else {
      height <- paste0(input$viz_height,"px")
    }
    d3Output("ont_nav",width = wid, height = height)
  })
  
  output$group_comp_ui <- renderUI({
    if (input$viz_wid == "" | grepl("\\D", input$viz_wid)){
      wid <- "100%"
    } else {
      wid <- paste0(input$viz_wid,"px")
    }
    
    if (input$viz_height == "" | grepl("\\D", input$viz_height)){
      height <- "840px"
    } else {
      height <- paste0(input$viz_height,"px")
    }
    d3Output("group_comp",width = wid, height = height)
  })
  
  # now send it off to d3!
  
  output$ont_nav <- renderD3({ 
    if ((!is.null(serverValues$restart) | !is.null(input$update_map)) &&
        length(serverValues$fc_groups) > 0 &&
        (serverValues$restart > 0 | (!is.null(input$update_map)  && input$update_map > 0))){
      
      r2d3(data=serverValues$heat.df, script = "js_scripts/all_tabs.js", d3_version = 5,
             options(r2d3.theme = list(
               background = serverValues$col_bg,
               foreground = serverValues$textcol_bg)),
             options = list(color_pal = unname(serverValues$color_pal),
                            color_pal_dark = unname(serverValues$color_pal_dark),
                            color_pal_chord = unname(serverValues$color_pal_chord),
                            color_pal_dark_chord = unname(serverValues$color_pal_dark_chord),
                            
                            time_points = serverValues$time_points,
                            genes = serverValues$genes,
                            heights = serverValues$heights,
                            chord_dat = serverValues$chord_dat,
                            fc_cats = serverValues$fc_groups,
                            fc_rep = serverValues$fc_rep,
                            from_fc_cat = serverValues$from_fc_cat,
                            from_genes = serverValues$from_genes,
                            tot_fc_cats = serverValues$tot_fc_cats,
                            go_name = serverValues$go_name,
                            
                            go_df = data_to_json(serverValues$go.df),
                            fc_cats_go = paste0("Tot_",serverValues$fc_groups),
                            fc_full_names = serverValues$fc_groups,
                            fc_names = serverValues$fc_abbrev,
                            sig = serverValues$sig,
                            
                            textcol_bg = serverValues$textcol_bg,
                            curr = curr,
                            font_size = input$font_size,
                            which_camp = "ont_nav",
                            path_trace = c(rbind(rev(path_trace),rep("^",length(path_trace))))[-(length(path_trace)*2)],
                            togg_path = serverValues$togg_path,
                            togg_border = input$togg_border,
                            togg_label = input$togg_label
                           )
      )
    }
  })
  
  output$group_comp <- renderD3({ 
    if ((!is.null(serverValues$restart) | (!is.null(input$update_map)  && input$update_map > 0)) &&
        length(serverValues$fc_groups) > 0 && 
        (serverValues$restart > 0 | (!is.null(input$update_map)  && input$update_map > 0)) && 
        sum(serverValues$tot_fc_cats) > 1){
      # print(sum(serverValues$chord_dat))
      
      r2d3(data=serverValues$heat.df, script = "js_scripts/all_tabs.js", d3_version = 5,
           options(r2d3.theme = list(
             background = serverValues$col_bg,
             foreground = serverValues$textcol_bg)),
           options = list(color_pal = unname(serverValues$color_pal),
                          color_pal_dark = unname(serverValues$color_pal_dark),
                          color_pal_chord = unname(serverValues$color_pal_chord),
                          color_pal_dark_chord = unname(serverValues$color_pal_dark_chord),
                          
                          time_points = serverValues$time_points,
                          genes = serverValues$genes,
                          heights = serverValues$heights,
                          chord_dat = serverValues$chord_dat,
                          fc_cats = serverValues$fc_groups,
                          fc_rep = serverValues$fc_rep,
                          from_fc_cat = serverValues$from_fc_cat,
                          from_genes = serverValues$from_genes,
                          tot_fc_cats = serverValues$tot_fc_cats,
                          go_name = serverValues$go_name,
                          fc_hours_shifted = serverValues$hm_hours_shifted,
                          tot_chords = sum(serverValues$chord_dat)/2,
                          
                          go_df = data_to_json(serverValues$go.df),
                          fc_cats_go = paste0("Tot_",serverValues$fc_groups),
                          fc_full_names = serverValues$fc_groups,
                          fc_names = serverValues$fc_abbrev,
                          sig = serverValues$sig,
                          
                          textcol_bg = serverValues$textcol_bg,
                          curr = curr,
                          font_size = input$font_size,
                          which_camp = "group_comp",
                          path_trace = c(rbind(rev(path_trace),rep("^",length(path_trace))))[-(length(path_trace)*2)],
                          togg_path = serverValues$togg_path,
                          togg_border = input$togg_border,
                          togg_label = input$togg_label
           )
      )
    }
  })
  
  output$ont_map <- renderD3({ 
    if (!is.null(input$update_map) && input$update_map > 0 & input$ont_in_map!= "NA" & !input$ont_in_map %in% c("biological_process", "molecular_function", "cellular_component")){
      r2d3(data = data_to_json(serverValues$dat_sankey), script = "js_scripts/all_tabs.js", 
           d3_version = 5, 
           dependencies = c("js_scripts/sankey/d3-sankey.js",
                            "js_scripts/sankey/d3-sankey.min.js"),
           options = list("all_sig"=serverValues$all_sig,
                          "sig_color"=unname(serverValues$sig_color),
                          "no_sig_color"=serverValues$no_sig_color,
                          "font_size" = input$font_size,
                          "which_camp" = "ont_map",
                          textcol_bg = serverValues$textcol_bg,
                          togg_border = input$togg_border,
                          togg_label = input$togg_label
                          ),
           options(r2d3.theme = list(
             background = serverValues$col_bg,
             foreground = serverValues$textcol_bg)))
    }
  })
  
  output$ont_map_below <- renderUI({
    if(!is.null(serverValues$ont_map_below_text)){
      if (serverValues$is_go){
        serverValues$ont_map_below_text
      } else {
        HTML("Click on a GO Term for a link to more information about that GO Term, as well as enrichment information for that term for the selected AC category.")
      }
    } else {
      HTML("Click on a GO Term for a link to more information about that GO Term, as well as enrichment information for that term for the selected AC category.")
    }
  })
  
  output$ont_nav_below <- renderUI({
    if(!is.null(serverValues$ont_nav_below_text)){
      if (serverValues$is_go){
        serverValues$ont_nav_below_text
      } else {
        HTML("Click on a GO Term for a link to more information about that GO Term. Click on an AC group's bar text to get enrichment information for that term for the selected AC category.")
      }
    } else  {
      HTML("Click on a GO Term for a link to more information about that GO Term. Click on an AC group's bar text to get enrichment information for that term for the selected AC category.")
    }
  })
  
  output$group_comp_below <- renderUI({
    if(!is.null(serverValues$group_comp_below_text)){
      serverValues$group_comp_below_text
    } else  {
      HTML("")
    }
  })
  
  output$frame <- renderUI({
    if(!is.null(serverValues$url)) {
      
      if (!serverValues$is_go){
        HTML(paste0("<center><h3>For more information on ",serverValues$term_of_interest,", please see: <a href='", serverValues$url,"' target='_blank'>UNIPROT</a>.</center></h3>"))
      } else {
        HTML(paste0("<center><h3>For more information on ",serverValues$term_of_interest,", please see: <a href='", serverValues$url,"' target='_blank'>QuickGO</a>.</center></h3>"))
      }
      # tags$iframe(src=serverValues$url, height = 1077, width = 700)
    } 
  })
  
  output$expr <- renderPlot({
    if(!is.null(serverValues$url) && !serverValues$is_go) {
      if (num_reps == 1){ # single replicate
        # generate a graph of gene expression with original data and fitted data
        
        # getting the total results: original and fitted values
        data.m <- data.frame(matrix(0,length(timen),3))
        colnames(data.m) <- c("Original","Fit","Times")
        data.m$Original <- as.numeric(serverValues$tr_sub[,c((end_num+1):(length(timen)+end_num))])
        data.m$Fit <- as.numeric(serverValues$tr_sub[,-c(1:(length(timen)+end_num))])
        data.m$Times <- timen
        
        # create gene expression plot
        col_vect <- c("Original"="black","Fit"="green")
        plot_viz <- ggplot(data = data.m,aes(x=Times))+
          geom_line(aes(y=Original,colour="Original"))+
          geom_line(aes(y=Fit,colour="Fit"))+
          scale_color_manual("",values=col_vect)+
          ggtitle(paste(input$gene_name))+
          theme(text= element_text(size = 20),plot.title = element_text(hjust = .5),
                legend.position = "bottom",legend.direction = "horizontal")+
          labs(x="Hours",y="Expression")
        
      } else{
        plot_viz<-ggplot(data = serverValues$ribbon.df,aes(x=Times))+ # declare the dataframe and main variables
          geom_ribbon(aes(x=Times, ymax=Max, ymin=Min, colour="Original"),
                      fill = "gray", alpha = 0.5)+ # create shading
          geom_line(aes(y=Fit,colour="Fit"))+ # fitted values
          ggtitle(paste(serverValues$tr_sub$`Gene Name`))+ # gene name is title
          scale_color_manual("",values=color_bar)+
          scale_fill_manual("",values=color_bar)+
          theme(text= element_text(size = 20),plot.title = element_text(hjust = .5),
                legend.position = "bottom",legend.direction = "horizontal")+
          labs(x="Hours", y="Expression") #Label for axes
        
        # add specific replicate lines 
        for (i in 1:num_reps){
          plot_viz <- plot_viz + geom_line(data = serverValues$ribbon.df,aes_string(x="Times",y=paste("Rep",i,sep = ".")), colour=color_bar[i], alpha=0.6)
        }
      }
      
      return(plot_viz)
    }
  })
  
  output$expr_text <- renderPrint({
    if(!is.null(serverValues$url) && !serverValues$is_go) {
      cat(paste("Gene Name:",serverValues$tr_sub$`Gene Name`,"\n"))
      cat(paste("Convergence:", serverValues$tr_sub$Convergence,"\n"))
      cat(paste("Iterations:",serverValues$tr_sub$Iterations,"\n"))
      cat(paste("Amplitude Change Coefficient:", serverValues$tr_sub$Amplitude.Change.Coefficient,"\n"))
      cat(paste("Oscillation Type:",serverValues$tr_sub$`Oscillation Type`,"\n"))
      cat(paste("Initial.Amplitude:", serverValues$tr_sub$Initial.Amplitude,"\n"))
      cat(paste("Radian.Frequency:",serverValues$tr_sub$Radian.Frequency,"\n"))
      cat(paste("Period:",serverValues$tr_sub$Period,"\n"))
      cat(paste("Phase Shift:",serverValues$tr_sub$`Phase Shift`,"\n"))
      cat(paste("Hours Shifted:",serverValues$tr_sub$`Hours Shifted`,"\n"))
      cat(paste("Slope:",serverValues$tr_sub$`Slope`,"\n"))
      cat(paste("P-Value:",serverValues$tr_sub$`P.Value`,"\n"))
      cat(paste("BH Adj P-Value:",serverValues$tr_sub$`BH.Adj.P.Value`,"\n"))
      cat(paste("BY Adj P-Value:",serverValues$tr_sub$`BY.Adj.P.Value`,"\n"))
    }
  })
  
  output$gene_list <- renderDataTable({
    if (!is.null(serverValues$restart) | (!is.null(input$update_map)  && input$update_map > 0) &&
         length(serverValues$fc_groups) > 0 &&
        (serverValues$restart > 0 | (!is.null(input$update_map)  && input$update_map > 0))){
      serverValues$gene.df
    }
    
  })
  
  output$org_tbl <- renderTable(all_org_info,  
                            striped = TRUE,  
                            width = '100%')
  
  output$user_data_info <- renderPrint({
    cat("")
    cat("User ENCORE Inputs:\n")
    cat(paste("ENCORE End Date and Time: ",user_input_ont$save_date,"\n"))
    cat(paste("Original File Name: ",user_input_ont$orig_filen,"\n"))
    cat(paste("Save File Name: ",user_input_ont$save_filen,"\n"))
    cat(paste("Organism: ",id_to_common[user_input_ont$org_tax],"\n")) # YOU SHOULD TAKE THE NAME OF THIS
    cat(paste("Gene Name ID Type: ",user_input_ont$id_type,"\n")) # YOU SHOULD TAKE THE NAME OF THIS
    cat(paste("Genes to Focus Enrichment On: ", paste(user_input_ont$gene_focus, collapse = " "),"\n"))
    cat(paste("Ontology Types to Compute: ",user_input_ont$ont_group,"\n"))
    cat(paste("ECHO P-Value Category: ",user_input_ont$pval_cat,"\n"))
    cat(paste("ECHO Significance Cutoff: ",user_input_ont$sig_level,"\n"))
    cat(paste("Ontology P-Value Category: ",user_input_ont$ont_pval_cat,"\n"))
    cat(paste("Ontology Significance Cutoff: ",user_input_ont$ont_sig_level,"\n"))
    cat(paste("Period Restriction, Low Value: ",user_input_ont$low_range,"\n"))
    cat(paste("Period Restriction, High Value: ",user_input_ont$high_range,"\n"))
  })
  
  #function to download png of plot
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0(serverValues$go_name, '.png') },
    content = function(file) 
    {
      png(file)
      print(monitor_viz)
      dev.off()
    },
    contentType='image/png'
  )# now send it off to d3!
  
  # backend for calculating ontology ----
  
  observeEvent(input$run_ont, withProgress(message = paste("Creating ENCORE File. Started on:",Sys.time()), value = 0, {
    # progress bar 
    incProgress(0/26, detail = paste("Loading selected data. Started on:",Sys.time()))
    
    missing_genes <- c()
    
    # get gene ids
    
    gene_focus <- strsplit(input$gene_focus,"\n", fixed=T)[[1]]
    gene_focus <- gene_focus[!duplicated(gene_focus)]
    
    # first, save all user inputs
    user_input_ont <<- list(
      "save_date" = Sys.time(),
      "orig_filen" = input$dat_ont$name,
      "save_filen" = input$save_filen,
      "org_tax" = input$org_tax,
      "id_type" = input$id_type,
      "pval_cat" = input$pval_cat,
      "sig_level" = input$sig_level,
      "ont_pval_cat" = input$ont_pval_cat,
      "ont_sig_level" = input$ont_sig_level,
      "low_range" = input$low_range,
      "high_range" = input$high_range,
      "ont_group" = input$ont_group,
      "gene_focus" = gene_focus
    )
    
    
    # save and load orignal file
    orig_filen <<- input$dat_ont$name 
    save_filen <<- input$save_filen
    load(input$dat_ont$datapath, envir = globalenv())
    
    # get subsets of certain range
    
    if (input$low_range == "" | input$high_range == ""){
      is_all_range <<- TRUE
      low_range <<- NA
      high_range <<- NA
    } else {
      is_all_range <<- FALSE
      low_range <<- as.numeric(sapply(input$low_range, function(x) eval(parse(text=x))))
      high_range <<- as.numeric(sapply(input$high_range, function(x) eval(parse(text=x))))
    }
    
    # specifics based on organism ----
    
    all_genes <- read.csv(paste0("data/",input$org_tax,"_background.csv"),header = T,stringsAsFactors = F)
    
    # remove duplicates
    total_results <- total_results[order(total_results$`BH Adj P-Value`),] # order by pvalue
    total_results <- total_results[!duplicated(total_results$`Gene Name`),] # remove duplicates
    # rename columns of total_results for ease
    colnames(total_results)[14:16] <- c("P.Value","BH.Adj.P.Value","BY.Adj.P.Value")
    
    orig_total_results <- total_results
    
    # add background genes not in list to total_results
    
    existing_bg <- all_genes[,input$id_type][all_genes[,input$id_type]!=""]
    background <- setdiff(existing_bg,total_results$`Gene Name`)
    if (length(background) > 0){
      orig_len <- nrow(total_results)
      emp <- data.frame(matrix(NA,
                               length(background),
                               ncol(total_results)))
      colnames(emp) <- colnames(total_results)
      total_results <- rbind(total_results, emp)
      total_results$`Gene Name`[-(1:orig_len)] <- background
      # since they were put at the end, we don't have to order
      total_results <- total_results[!duplicated(tolower(total_results$`Gene Name`)),] # remove duplicates
    }
    # SOME MAPPING
    
    incProgress(1/26, detail = paste("Mapping gene names. Started on:",Sys.time()))
    
    if (input$org_tax == "5141"){ # neurospora has no mapping
      map_sub.df <- data.frame(matrix(0,nrow(total_results),2))
      colnames(map_sub.df) <- c("query","entrezgene")
      # entrez gene is a misnomer for neurospora
      map_sub.df$query <- map_sub.df$entrezgene  <- total_results$`Gene Name`
    } else { 
      # switch statement for different organism types
      # get the entrez genes - named vector
      if (input$id_type != 'ENTREZID'){
        switch(input$org_tax,
               "10090" = { # mouse
                 entrezgene <- mapIds(org.Mm.eg.db,
                                      total_results$`Gene Name`,
                                      'ENTREZID', input$id_type,
                                      multiVals = "first")
               },
               "9606" = { # human
                 entrezgene <- mapIds(org.Hs.eg.db,
                                      total_results$`Gene Name`,
                                      'ENTREZID', input$id_type,
                                      multiVals = "first")
               },
               "7227" = { # drosophila
                 entrezgene <- mapIds(org.Dm.eg.db,
                                      total_results$`Gene Name`,
                                      'ENTREZID', input$id_type,
                                      multiVals = "first")
               },
               "7165" = { # anopheles
                 entrezgene <- mapIds(org.Ag.eg.db,
                                      total_results$`Gene Name`,
                                      'ENTREZID', input$id_type,
                                      multiVals = "first")
               },
               "4932" = { # yeast
                 entrezgene <- mapIds(org.Sc.sgd.db,
                                      total_results$`Gene Name`,
                                      'ENTREZID', input$id_type,
                                      multiVals = "first")
               },
               "511145" = { # e coli
                 entrezgene <- mapIds(org.Hs.eg.db,
                                      total_results$`Gene Name`,
                                      'ENTREZID', input$id_type,
                                      multiVals = "first")
               })
      } else {
        entrezgene <- total_results$`Gene Name`
        names(entrezgene) <- entrezgene
      }
      
      # put into data frame form
      map_sub.df <- data.frame(matrix(0,length(entrezgene),2)) # query and entrezid
      colnames(map_sub.df) <- c("query","entrezgene")
      map_sub.df$query <- names(entrezgene)
      map_sub.df$entrezgene <- as.character(lapply(entrezgene, `[[`, 1))
      
      # first sort the data by numerical ids - why?
      map_sub.df <- map_sub.df[order(map_sub.df$entrezgene),]
      map_sub.df <- map_sub.df[!(duplicated(map_sub.df$query)),] # now remove duplicates
      # order both total results and map sub, in order to add relevant columns
      total_results <- total_results[order(total_results$`Gene Name`),]
      map_sub.df <- map_sub.df[order(map_sub.df$query),]
      
      # remove genes with no entrez id mapping, but save them first
      missing_genes <- c(missing_genes, total_results[is.na(map_sub.df$entrezgene),"Gene Name"])
      total_results <- total_results[!is.na(map_sub.df$entrezgene),]
      map_sub.df <- map_sub.df[!is.na(map_sub.df$entrezgene),]
    }
    
    map_sub.df <- cbind(map_sub.df, 
                        "Osc.Type" = total_results$`Oscillation Type`, 
                        "Period" = total_results$Period,
                        "BH Adj P-Value" = total_results$`BH.Adj.P.Value`,
                        "BY Adj P-Value" = total_results$`BY.Adj.P.Value`,
                        "P-Value" = total_results$`P.Value`)
    
    incProgress(1/26, detail = paste("Mapping STRING names. Started on:",Sys.time()))
    
    # STRINGDB
    string_db <- STRINGdb$new(version="10", # most recent available
                              species=as.numeric(input$org_tax), # except the e coli???
                              score_threshold=0, 
                              input_directory="")
    map_sub.df <- string_db$map(map_sub.df, "entrezgene", removeUnmappedRows = F)
    # remove missing and duplicated genes 
    missing_genes <- c(missing_genes, total_results[is.na(map_sub.df$entrezgene),"Gene Name"])
    total_results <- total_results[!is.na(map_sub.df$entrezgene),]
    map_sub.df <- map_sub.df[!is.na(map_sub.df$entrezgene),]
    
    # not removing missing string ids -- hurts ontologies a lot
    # missing_genes <- c(missing_genes, total_results[is.na(map_sub.df$STRING_id),"Gene Name"])
    # total_results <- total_results[!is.na(map_sub.df$STRING_id),]
    # map_sub.df <- map_sub.df[!is.na(map_sub.df$STRING_id),]
    
    total_results <- total_results[order(total_results$`P.Value`),]
    map_sub.df <- map_sub.df[order(map_sub.df$`P.Value`),]
    total_results <- total_results[!duplicated(map_sub.df$entrezgene),]
    map_sub.df <- map_sub.df[!duplicated(map_sub.df$entrezgene),]
    
    # push map_sub.df to global to save
    map_sub.df <<- map_sub.df
    
    # update gene focus if none specified
    if (sum(nchar(gene_focus)) == 0){
      gene_focus <- map_sub.df$query
    }
    
    # get gene ontology enrichments for each fc category
    
    ont_list <- get_fc_results(input$org_tax,
                               input$pval_cat,
                               input$sig_level,
                               input$ont_pval_cat,
                               input$ont_sig_level,
                               gene_focus,
                               input$ont_group)
    
    fc_results <<- ont_list[[1]]
    no_sig <<- ont_list[[2]]
    total_results <<- orig_total_results # just because of name changes
    
    # now put missing genes for output upon download
    missing_genes <<- intersect(orig_total_results$`Gene Name`, missing_genes)
    
    output$dropped_genes <- renderText(c(paste0("Unmapped genes (",length(missing_genes),"):\n"),
                                         missing_genes))
    
    incProgress(1/26, detail = paste("Finished! Started on:",Sys.time()))
  }))
  
  # rendering output for ENCORE file ----
  
  observeEvent(input$ont_file,withProgress(message = "Saving ENCORE file!",detail = paste("Started on:",Sys.time()),value = 0, {
    incProgress(1/2)
    
    output$download_ont_file <- renderText({
      "Files will appear in the 'downloads' subfolder of the ENCORE folder."
    })
    
    if (is_all_range){
      filen <- (paste0("downloads//",input$save_filen,'_ENCORE.RData'))
    } else {
      filen <- (paste0("downloads//",input$save_filen,"_",low_range,"_to_",high_range,"_ENCORE.RData"))
    }
    
    save(file=filen,
         list=c("total_results","fc_results","orig_filen", "map_sub.df", "timen", "num_reps", "no_sig","is_all_range","low_range","high_range","user_input_ont", "user_input","missing_genes"),
         file)
    
    incProgress(1/2)
    
    
  }))
  
  # output$ont_file <- downloadHandler( # Visualization results
  #   filename = function() { 
  #     if (is_all_range){
  #       return(paste0(input$save_filen,'_ENCORE.RData'))
  #     } else {
  #       return(paste0(input$save_filen,"_",low_range,"_to_",high_range,"_ENCORE.RData"))
  #     }
  #     },
  #   content = function(file) {
  #     
  #       print(paste("Saving file... Started on:",Sys.time()))
  #       save(file=file,
  #            list=c("total_results","fc_results","orig_filen", "map_sub.df", "timen", "num_reps", "no_sig","is_all_range","low_range","high_range","user_input_ont", "missing_genes"),
  #            file)
  #       print(paste("Finished!:",Sys.time()))
  #     
  #   }
  # )
  
  output$org_id_tbl <- renderTable(avail_id_types,  
                                     striped = TRUE,  
                                     width = '100%')
  
  output$org_id_ex_tbl <- renderTable(org_id_ex,  
                                      striped = TRUE,  
                                      width = 'auto')
  
  
}

# run ----

shinyApp(ui, server)