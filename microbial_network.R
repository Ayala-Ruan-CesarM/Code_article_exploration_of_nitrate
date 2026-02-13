library(tidyverse)
library(NetCoMi)
library(igraph)
library(phyloseq)

file_paths <- list.files(path = "Simulations_genera_level/",
                         pattern = "PWIW_simulated_neg_binom_05_rsed57_n200.tsv",
                         full.names = TRUE)
for (file in file_paths) {
  
  # 2a: Read the TSV file
  cat("Processing file:", basename(file), "\n")
  otu_data <- read_tsv(file, show_col_types = FALSE)
  
  # 2b: Process the OTU data into a phyloseq-compatible matrix
  otu_matrix <- otu_data %>%
    column_to_rownames("#OTU_ID")
  
  # 2c: Create the phyloseq OTU table
  otu_table <- otu_table(as.matrix(otu_matrix), taxa_are_rows = TRUE)
  
  # 2d: Dynamically generate sample metadata from the column names
  # The column names, excluding '#OTU_ID', are your sample IDs
  sample_names <- colnames(otu_matrix)
  
  # 2e: Extract the 'Treatment' (PW or IW) from the sample names
  treatment <- ifelse(grepl("PW_", sample_names), "PW", "IW")
  
  # 2f: Create the sample data frame
  sample_info <- data.frame(
    SampleID = sample_names,
    Treatment = treatment,
    row.names = sample_names
  )
  
  # 2g: Create the sample data object for phyloseq
  sample_data <- sample_data(sample_info)
  
  # Extract the OTU IDs from the original data frame
  otu_ids <- otu_data$`#OTU_ID`
  
  # Create a new data frame for the taxonomic information
  taxa_df <- data.frame(
    Rank6 = otu_ids,
    Rank7 = NA, # You can fill this with other data or leave it as NA
    stringsAsFactors = FALSE
  )
  
  # Set the row names of the tax_table to match the OTU IDs
  rownames(taxa_df) <- otu_ids
  
  # Create the phyloseq 'tax_table' object
  tax_table_object <- tax_table(as.matrix(taxa_df))
  
  # 2h: Create the final phyloseq object
  phyloseq_object_test <- phyloseq(otu_table, sample_data, tax_table_object)
  
  # You can now save or analyze this object
  # For example, to save the object for later use:
  object_name <- paste0("phyloseq_object_", gsub("\\.tsv$", "", basename(file)))
  assign(object_name, phyloseq_object_test)

}
#print(phyloseq_object_test)

# Subset the phyloseq object to only include the PW samples without QC filterd
physeq_subset_PW_200 <- subset_samples(phyloseq_object_test, Treatment == "PW")
#saveRDS(physeq_subset_PW_200, "Simulations_genera_level/physeq_subset_PW_200")
physeq_filtered_PW <- filter_taxa(physeq_subset_PW_200, function(x) sum(x > 0) > (0.2 * length(x)), TRUE)

# Construct the network using the SPRING method
net_spring_high_820freq_rep100  <- netConstruct(
  data = physeq_subset_PW_200,
  taxRank = "Rank6",
  cores = 12,
  filtTax = "highestFreq",
  filtTaxPar = list(highestFreq = 820), # 275 > 10000, 415>5000, 820>1000
  measure = "spring",
  measurePar = list(nlambda=30, 
                    rep.num=100,
                    Rmethod = "aprox",
                    verboseR = TRUE,
                    quantitative = TRUE),
  normMethod = "clr", 
  zeroMethod = "pseudo",
  sparsMethod = "none", 
  dissFunc = "signed",
  verbose = 3,
  seed = 57)

net_analyzed_820_f2 <- netAnalyze(net_spring_high_820freq_rep100,
                                  entrLCC = TRUE,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector",
                                  weightDeg = FALSE, 
                                  normDeg = FALSE) 

################ All Taxa plot. Annex Fig. ################
                                  
png("FigAnnex_network_820.png", 
    width = 2000, 
    height = 1600, 
    res = 300)
                                  
 p_F2_820 <- plot(net_analyzed_820_f2,
                 edgeWidth = 0.05,
                 edgeTranspLow = 50,
                 edgeTranspHigh = 50, # THe higher the more transparent it is
                 labels = FALSE,
                 rmSingles = TRUE)
dev.off()
################ Fig. 5 on the manuscript ################
png("Fig5_network_820.png", 
    width = 2000, 
    height = 1600, 
    res = 300)
                                  
p_F2_820 <- plot(net_analyzed_820_f2,
                 nodeFilter = "highestDegree",
                 nodeFilterPar = 75,
                 nodeSize = "degree",
                 edgeWidth = 1,
                 edgeTranspLow = 50,
                 edgeTranspHigh = 20, # THe higher the more transparent it is
                 cexLabels = 2.0,
                 cexNodes = 1.5,
                 labels = TRUE,
                 rmSingles = TRUE,
                 mar = c(1, 3, 4, 8)
                ) 

legend(0.8, 0.9, cex = 0.5 , title = "Estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)

#1. Get the cluster IDs and node colors
cluster_data <- net_analyzed_820_f2$clustering$clust1
plot_colors <- p_F2_820$nodecolor$nodecol1

# 2. Get the names of the nodes that were plotted
plotted_node_names <- names(plot_colors)

# 3. Filter the cluster data to include only the nodes that were plotted
aligned_cluster_data <- cluster_data[match(plotted_node_names, names(cluster_data))]

# 4. Create a data frame to combine the aligned vectors
node_info <- data.frame(
  cluster_id = aligned_cluster_data,
  color = plot_colors
)

# 5. Get the unique cluster IDs and their corresponding colors
legend_data <- node_info[!duplicated(node_info$cluster_id), ]

# Order the legend data by cluster ID
legend_data <- legend_data[order(legend_data$cluster_id), ]

# Get the final lists for the legend function
cluster_names <- paste("Cluster", legend_data$cluster_id)
cluster_colors <- legend_data$color

legend(0.9, 0.15,
       title = "Clusters",
       legend = cluster_names,
       col = cluster_colors,
       pch = 19,
       bty = "n",
       cex = 0.6)

dev.off()
