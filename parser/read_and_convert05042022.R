##################################################
## Project: COVID-19 Disease Map
## Script purpose: Translate COVID-19 Disease Map diagrams to HiPathia 
## using a dedicated converter
## Date: 05.04.2022
## Author: Marek Ostaszewski
##################################################
library(httr)
source("https://gitlab.lcsb.uni.lu/marek.ostaszewski/systems_biology_translators/-/raw/master/hipathia/convert_to_hipathia.R")
source("https://gitlab.lcsb.uni.lu/marek.ostaszewski/systems_biology_translators/-/raw/master/hipathia/reduce_hipathia_graph.R")

setwd("/home/m3m/INFO_PROYECTO/Coronavirus/converter_2HiPathia/SIF_05042022/")

### A wrapper function for conversion/reduction calls
### Uses MINERVA project url, writes results to target directory,
### 'use_model_name' allows to use diagram name instead of id when writing
### 'plot_graph' plots the outcome of the conversion/reduction
convert_minerva_to_hipathia <- function(minerva_url, project = NULL, target_directory, 
                                        use_model_name = F, plot_graph = F) {
  if(!dir.exists(target_directory)) { dir.create(target_directory, recursive = T) }
  message("Converting ", minerva_url)
  components <- get_map_components(minerva_url, project)

  if(is.null(project)) { project <- get_default_project(minerva_url) }
  
  for(mid in components$models$idObject) {
    message("Processing the diagram id: ", mid)
    
    ### Translate the diagram
    message(paste0("Creating casq sif file ..."))
    raw_sif <- get_casq_sif(minerva_url, project, mid, target_directory)
    if(is.null(raw_sif)) {
      ### Skip if the diagram was empty after running casq
      message("casq produced an empty file, skipping...")
      next
    }
    hid <- paste0("hsa", mid) ### Hipathia prefix, add species here
    
    ### Create the hipathia diagram
    message(paste0("Creating HiPathia sif ..."))
    raw_graph <- get_hipathia_files(raw_sif = raw_sif,
                                    node_prefix = hid,
                                    model_elements = components$map_elements[[as.character(mid)]],
                                    target_directory)
    ### Reduce the diagram
    message(paste0("Cleaning HiPathia sif ..."))
    hipathia_graph <- read_hipathia_graph(hid, target_directory)
    reduced_hipathia_graph <- process_hipathia_graph(hipathia_graph, 
                                                     verbose = F)
    ### Plot if needed
    if(plot_graph) { plot_hipathia(reduced_hipathia_graph) }
  
    if(is_hipathia_graph(reduced_hipathia_graph)) {
      ### Write to target path
      write_hipathia_graph(reduced_hipathia_graph, split_effectors = T,
                           paste0(target_directory,ifelse(use_model_name, 
                                                          with(components$models, name[idObject == mid]), hid)))
      
      name.pathways_hsa_file <- data.frame(code = hid, name = components$models$name[components$models$idObject == mid])
      write.table(name.pathways_hsa_file, file = paste0(target_directory,"name.pathways_hsa.txt"), append = T, sep = "\t", col.names = F, row.names = F, quote = F )
      
    } else {
      message("The produced graph has a structure incompatible with HiPathia.\nNo output written.")
    }
  }
}


write_hipathia_graph <- function(graph, file_name, split_effectors = FALSE) {
  att_file <- data.frame(get.vertex.attribute(graph))
  sif_file <- data.frame(V1 = get.edgelist(graph)[,1], V2 = get.edge.attribute(graph, "sign"), V3 = get.edgelist(graph)[,2])

  if(split_effectors) { ### If requested, create two separate tables for the core network and functions
    ### Split the annotation file
    functions_att_file <- dplyr::filter(att_file, type == "function")
    att_file <- dplyr::filter(att_file, type != "function")
    ### Split the sif file
    functions_sif_file <- dplyr::filter(sif_file, V3 %in% functions_att_file$name)
    sif_file <- dplyr::filter(sif_file, !(V3 %in% functions_att_file$name))

    write.table(functions_att_file, file = paste0(file_name,"_functions.att"),
                sep = "\t", quote = F, row.names = F)
    message("write_hipathia_files - function attribute table write: ", paste0(file_name,"_functions.att"))

    cat(paste0("\n", paste(colnames(functions_sif_file), collapse = "\t")),
        file = paste0(file_name,"_functions.sif"))
    write.table(functions_sif_file, file = paste0(file_name,"_functions.sif"),
                sep = "\t", quote = F, row.names = F, col.names = F, append = T)
    message("write_hipathia_files - sif file write: ", paste0(file_name,"_functions.sif"))
  }

  write.table(att_file, file = paste0(file_name,".att"),
              sep = "\t", quote = F, row.names = F)
  message("write_hipathia_files - attribute table write: ", paste0(file_name,".att"))

  cat(paste0("\n", paste(colnames(sif_file), collapse = "\t"),"\n"),file = paste0(file_name,".sif"))
  write.table(sif_file, file = paste0(file_name,".sif"),
              sep = "\t", quote = F, row.names = F, col.names = F, append = T)
  message("write_hipathia_files - sif file write: ", paste0(file_name,".sif"))
  return(TRUE)
}


### The address of the COVID-19 map in MINERVA
c19dmap <- "https://covid19map.elixir-luxembourg.org/minerva/api/"

### Individual diagrams
convert_minerva_to_hipathia(c19dmap, target_directory = "./minerva/",
                            use_model_name = F)

### WikiPathways diagrams, internal representation
convert_minerva_to_hipathia(c19dmap, project = "internal_wikipathways_selected", target_directory = "./_notgit/wp/", 
                            use_model_name = T)

### Reactome diagrams, internal representation
convert_minerva_to_hipathia(c19dmap, project = "internal_reactome_selected", target_directory = "./_notgit/rc/",
                            use_model_name = T)
