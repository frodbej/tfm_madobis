# Loading required packages

library(ggplot2)

# Loading data

hit_table_classifier <- read.csv("hit_table_classifier.csv")
rownames(hit_table_classifier) <- hit_table_classifier$Image_Metadata_Gene
hit_table_classifier$Image_Metadata_Gene <- NULL

sum(hit_table_classifier$Long_projections..Cell.Count)
sum(hit_table_classifier$Total..Cell.Count)

# We have 6671 cells with long projections out of 57718 cells in our screening

summary(hit_table_classifier$Long_projections..Cell.Count)
sd(hit_table_classifier$Long_projections..Cell.Count)

summary(hit_table_classifier$Enriched.Score.long_projections)
# Enrichment score ranges from -1.46 to 1.36

# To study in more detail gene conditions with a high enrichment score in long projections
# cells we filter by the 95th percentile

p95 <- quantile(hit_table_classifier$Enriched.Score.long_projections, 0.95)

ggplot(hit_table_classifier, aes(x=Enriched.Score.long_projections)) + 
  geom_histogram(color="black", fill="white", binwidth = 0.1) + ylab("Count") +
  xlab("Enrichment Score") + 
  geom_vline(aes(xintercept=p95),
             color="blue", linetype="dashed", size=1)


high_enrichment <- hit_table_classifier[hit_table_classifier$Enriched.Score.long_projections > p95,]

rownames(high_enrichment)[which.max(high_enrichment$Enriched.Score.long_projections)]

# MYO5B condition has the highest enrichment in long projections cells

high_enrichment <- high_enrichment[order(-high_enrichment$Enriched.Score.long_projections),]
high_enrichment["Genes"] <- row.names(high_enrichment)

ggplot(data=high_enrichment, aes(x=Genes,y=Enriched.Score.long_projections)) +
  geom_bar(stat="identity") + scale_x_discrete(limits = high_enrichment$Genes) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Long Projections Enrichment Score")


# Save gene names to enter in STRING
genes_long_proj <- high_enrichment$Genes
write.table(genes_long_proj,"genes_long_projections.csv", quote = F, sep = "\t",
            row.names = F, col.names = F)
