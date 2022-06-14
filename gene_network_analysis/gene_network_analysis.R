# Loading required packages

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(igraph)
library(corrplot)

# We extracted all genes involved on the siRNA screening
genes <- read.table("genes.txt")
genes <- genes$V1

### Enrichment analysis of genes in study

GO_BP <- enrichGO(gene = genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

bp_bp <- barplot(GO_BP) + ggtitle("Biological Process")
bp_bp
dotplot(GO_BP)
cnetplot(GO_BP)

# We have a great enrichment of biological processes associated with actin cytoskeleton
# organization as we expected

GO_MF <- enrichGO(gene = genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = 'SYMBOL',
                  ont = "MF",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

bp_mf <- barplot(GO_MF) + ggtitle("Molecular Funcion")
bp_mf
dotplot(GO_MF)
cnetplot(GO_MF)

# We can see two principal molecular functions highly enriched (actin binding and
# GTPase regulator activity). These two clusters are related since actin cytoskeleton
# regulation is mostly mediated by Rho-GTPases family.

GO_CC <- enrichGO(gene = genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = 'SYMBOL',
                  ont = "CC",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

bp_cc <- barplot(GO_CC) + ggtitle("Cellular Component")
bp_cc
dotplot(GO_CC)
cnetplot(GO_CC)

# We can see that cellular compartments enriched are related to cytoskeletal functions
# as cell cortex, lamellipodium, focal adhesions...


### Building a gene network

# We use our list of genes involved on the study to build an interaction network
# using STRING. As all genes are involved in cytoskeletal function and our network
# is highly connected (507 nodes and 7625 edges) we decided to be very restrictive
# with the interactions and we build a network showing only interactions with the
# highest confidence (0.9) and only experiments and database validated interactions.
# We downloaded STRING network and loaded it on CytoScape for further analysis and
# representation of the network.
# After this, we load the network on R for statistical analysis using igraph package


### Reading graph file
gene_network <- read.graph("gene_network_reduced.graphml", format = "graphml")
gene_network
# Our network is composed by 135 nodes and 395 edges

### Network topology analysis

diameter(gene_network, directed = FALSE)
radius(gene_network, mode = "all")

path_length <- path.length.hist(graph = gene_network,directed = FALSE)
hist(path_length$res, freq = TRUE, main = "Shortest path length distribution")

network_degrees <- degree(gene_network)
network_degrees
degree_histogram <- hist(network_degrees,freq = F,col="blue",
                         xlab="Node degree", ylab="Probability",
                         main="Degree distribution")
# Network degree distribution appears to fit with a scale-free network, with
# a huge amount of nodes with low degree and some nodes with high degree (hubs).

### Scale-free network test
# On a scale-free network the network degree distribution follows a power-law
# distribution
# To check if our network is a scale-free network we base on linear regresion and 
# Kolmogorov-Smirnov test

### Linear regresion
degree.frequencies <- table(network_degrees)
degrees <- as.numeric(names(degree.frequencies))

# Log transformation
log10.degrees.frequencies <- log10(degree.frequencies)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies)))

lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)

n <- lm.r$coefficients[[1]]
m <- lm.r$coefficients[[2]]

# Representation
plot(log10.node.degrees,log10.degrees.frequencies,pch=1,cex=1.25,lwd=3,col="blue",
     xlab="log10(k)",ylab="log10(p(k))",cex.lab=1.5,cex.axis=1.25,font.axis=2)
lines(log10.node.degrees,m*log10.node.degrees+n,col="red",lwd=5)

summary(lm.r)
# H0: There is no relationship between response and predictor variables
# H1: There exists a relationship between response and predictor variables

# R^2 = 0.737 and p-value <<< 0.05. Our gene network fits to a scale-free network.


### Kolmogorov-Smirnov test
# H0: Node degree distribution fits a power-law distribution
# H1: Node degree distribution doesn't fit a power-law distribution
network.degree.distribution <- degree.distribution(gene_network)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]
# p-value >>> 0.05. We don't reject null hypothesis

# Based on these two tests, we can affirm that our gene network is a scale-free network

### Hubs identification

p95 <- quantile(network_degrees, probs=0.95)
hubs <- which(network_degrees > p95)
length(hubs)
hubs_gene_names <- names(hubs)
hubs_gene_names

# We have 7 hubs on our network. Hub genes are ACTA1, PXN, PAK1, PAK2, PTK2, PAK3
# & CDC42

hubs.enrichGO <- enrichGO(gene = hubs_gene_names, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

barplot(hubs.enrichGO)
dotplot(hubs.enrichGO)
cnetplot(hubs.enrichGO)

### Closeness or centrality

closeness <- closeness(gene_network, mode = "all")
p95 <- quantile(closeness, 0.95)
closeness.p95 <- which(closeness > p95)
length(closeness.p95)
closeness.gene.names <- names(closeness.p95)
closeness.gene.names

# We have 9 genes with high closeness: ACTA1, PXN, PAK1, PAK2, PTK2, PAK3, ARHGEF7,
# CDC42 & RAC1

closeness.enrichGO <- enrichGO(gene = closeness.gene.names, OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL",ont = "BP", pAdjustMethod = "BH",
                               qvalueCutoff = 0.05)

barplot(closeness.enrichGO)
dotplot(closeness.enrichGO)
cnetplot(closeness.enrichGO)

### Betweenness

node.betweenness <- betweenness(gene_network, normalized = TRUE)
p95 <- quantile(node.betweenness, 0.95)
betweenness.p95 <- which(node.betweenness > p95)
length(betweenness.p95)
betweenness.gene.names <- names(betweenness.p95)
betweenness.gene.names

# We have 8 genes with high betweenness: PXN, PAK1, PAK2, PTK2, PAK3, ARHGEF7,
# CDC42 & RAC1

betweenness.enrichGO <- enrichGO(gene = betweenness.gene.names, OrgDb = org.Hs.eg.db,
                                 keyType = "SYMBOL",
                                 ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

barplot(betweenness.enrichGO)
dotplot(betweenness.enrichGO)
cnetplot(betweenness.enrichGO)



# Plotting network data extracted from CytoScape analysis

network_data <- read.csv("network_data_wo_low_connectivity.csv")

str(network_data)

# Delete repeated and useless variables for our analysis
del_variables <- c(6,8,11,13,15,16,17)
network_data <- network_data[-del_variables]
summary(network_data)

# Correlation between variables without names column
cor_matrix <- round(cor(network_data[-7]),3)

cor_test <- cor.mtest(network_data[-7], conf.level = 0.95)


corrplot(cor_matrix, p.mat = cor_test$p, diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = "AOE",
         tl.col="black", tl.srt=45)

# Plot Closeness centrality and Stress

ggplot(network_data, aes(x=ClosenessCentrality,y=Stress, size=Degree, color=Degree)) +
  geom_point() +
  geom_text(data=subset(network_data, Stress > 20000 & ClosenessCentrality > 0.2),
            aes(ClosenessCentrality,Stress,label=name), hjust=1.1,vjust=1.1) +
  xlab("Closeness Centrality") + ylab("Stress") +
  guides(size = FALSE) + scale_color_gradient(low="blue", high="red")




