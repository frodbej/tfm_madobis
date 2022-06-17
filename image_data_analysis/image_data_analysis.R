### Loading required packages

library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(datawizard)
library(gplots)

### Loading data

metadata <- read.csv("Metadata_ID_Exp_ABD_Final.csv")
per_image_data <- read.csv("per_image.csv")

# Distribution of genes in the plates

sup_metadata <- metadata[metadata$ï..Experiment=="A" & metadata$Channel=="DAPI",]
sup_p_1 <- sup_metadata[sup_metadata$Plate==1,]
sup_p_2 <- sup_metadata[sup_metadata$Plate==2,]

sort(table(sup_p_1$Gene)/2, decreasing = TRUE)
sort(table(sup_p_2$Gene)/2, decreasing = TRUE)

# The genes are not evenly distributed on the plates


# As control, we check we have the same number of cells than nuclei
sum(per_image_data$Image_Count_Cell)
sum(per_image_data$Image_Count_FilteredNuclei)

###PCA

# We have to remove logical and character data but maintaining metadata as
# qualitative supplementary variables

df_for_pca <- per_image_data[,-which(sapply(per_image_data, class) == "logical")]
metadata_variables <- df_for_pca[,grepl("Image_Metadata",names(df_for_pca))]
df_for_pca <- df_for_pca[,-which(sapply(df_for_pca, class) == "character")]
df_for_pca <- df_for_pca[,-which(grepl("Image_Metadata",names(df_for_pca)))]
df_for_pca <- data.frame(df_for_pca,metadata_variables)
quali.sup <- which(grepl("Image_Metadata",names(df_for_pca)))

res.pca <- PCA(df_for_pca, scale.unit=TRUE, quali.sup = quali.sup, graph = FALSE)
fviz_pca_ind(res.pca, col.ind = df_for_pca$Image_Metadata_Condition, addEllipses = T)
fviz_pca_ind(res.pca, col.ind = factor(df_for_pca$Image_Metadata_Plate), addEllipses = T)
fviz_pca_ind(res.pca, col.ind = df_for_pca$Image_Metadata_Experiment, addEllipses = T)

# We see there are no differences between experiments and between plates 1 and 2
# because confidence ellipses overlap, so we discard batch effect


### Data exploration

metadata_genes <- unique(metadata$Gene)
length(metadata_genes)
# We are studying 521 genes according to metadata
write.table(metadata_genes,"genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)

genes <- unique(per_image_data$Image_Metadata_Gene)
n_genes <- length(genes)
n_genes
# We have at least one image of each gene condition in our three experiments

total_genes <- table(metadata$Gene)
# Divide by 3 because each gene is repeated 3 times on metadata, one time for each
# of the three channels
total_genes <- total_genes/3
available_genes <- table(per_image_data$Image_Metadata_Gene)
available_genes

table(per_image_data$Image_Metadata_Condition)

# Sort genes to see which ones we have more information about
sort(total_genes, decreasing = TRUE)[1:10]
# We take 798 images of Mock wells
sort(available_genes, decreasing = TRUE)[1:10]
# We get 703 images of Mock wells with at least one cell

# We build a data frame with images expected and observed for each gene

genes_in_study <- data.frame(total_genes,row.names = names(total_genes))[-1]
colnames(genes_in_study) <- "Total"

# For loop iterates on all expected genes, if we have images of wells with this gene
# we insert the number of available images. If we don't have images for
# that gene we insert zero.

for (i in 1:length(metadata_genes)){
  if (metadata_genes[i] %in% names(available_genes)){
    genes_in_study[metadata_genes[i],"Studied"] = available_genes[metadata_genes[i]]
  }else{
    genes_in_study[metadata_genes[i],"Studied"] = 0
  }
}

# Calculate the ratio of observed and expected images

genes_in_study <- genes_in_study %>% 
  mutate(Proportion = genes_in_study$Studied/genes_in_study$Total)


# Plot sorted proportions to see the information available for each of the genes

ggplot(data=genes_in_study, aes(x=1:n_genes,y=sort(Proportion,decreasing = T))) +
  geom_line() + geom_point() + xlab("Genes") + ylab("Proportion of available information")


sum(genes_in_study$Proportion == 1)/n_genes

# We have all the information available for more than half of the genes

sum(genes_in_study$Proportion < 0.5)/n_genes

# We have less than 2% of the genes with less than half of expected information

### Plot control conditions

# siTOX siRNA targets an antiapoptotic protein so if our screening is OK we should
# see less cells on these wells than on control wells.

control_condition <- per_image_data[per_image_data$Image_Metadata_Condition == "Control",]
nuclei_count_control <- control_condition$Image_Count_FilteredNuclei
sum(nuclei_count_control)

siTOX_condition <- per_image_data[per_image_data$Image_Metadata_Gene == "siTOX",]
nuclei_count_siTOX <- siTOX_condition$Image_Count_FilteredNuclei
sum(nuclei_count_siTOX)

boxplot(nuclei_count_control,nuclei_count_siTOX, ylab = "Cell count per image",
        names = c("Control", "siTOX"))

sup_df_control <- data.frame("Image_Count_FilteredNuclei"=nuclei_count_control,
                             "Condition"=rep("Control",length(nuclei_count_control)))
sup_df_siTOX <- data.frame("Image_Count_FilteredNuclei"=nuclei_count_siTOX,
                           "Condition"=rep("siTOX",length(nuclei_count_siTOX)))
sup_df <- rbind(sup_df_control,sup_df_siTOX)
sup_df$Condition <- as.factor(sup_df$Condition)


ggplot(sup_df, aes(x=Condition, y=Image_Count_FilteredNuclei, fill=Condition)) + 
  geom_violin() + geom_boxplot(width=0.1) + scale_fill_brewer(palette="Dark2") +
  ylab("Cells per image")


t.test(nuclei_count_control,nuclei_count_siTOX)
t.test(nuclei_count_control,nuclei_count_siTOX, alternative = "greater")

# There are significant differences between the means


# We check KIF11 and INCENP controls according to the area of nuclei

control_nuclei_area <- control_condition$Mean_FilteredNuclei_AreaShape_Area

KIF11_condition <- per_image_data[per_image_data$Image_Metadata_Gene == "KIF11",]
KIF11_nuclei_area <- KIF11_condition$Mean_FilteredNuclei_AreaShape_Area

INCENP_condition <- per_image_data[per_image_data$Image_Metadata_Gene == "INCENP",]
INCENP_nuclei_area <- INCENP_condition$Mean_FilteredNuclei_AreaShape_Area


boxplot(control_nuclei_area,KIF11_nuclei_area,INCENP_nuclei_area)

sup_df_control <- data.frame("Mean_FilteredNuclei_AreaShape_Area"=control_nuclei_area,
                             "Condition"=rep("Control",length(control_nuclei_area)))
sup_df_KIF11 <- data.frame("Mean_FilteredNuclei_AreaShape_Area"=KIF11_nuclei_area,
                             "Condition"=rep("KIF11",length(KIF11_nuclei_area)))
sup_df_INCENP <- data.frame("Mean_FilteredNuclei_AreaShape_Area"=INCENP_nuclei_area,
                             "Condition"=rep("INCENP",length(INCENP_nuclei_area)))
sup_df <- rbind(sup_df_control,sup_df_KIF11,sup_df_INCENP)
sup_df$Condition <- as.factor(sup_df$Condition)


ggplot(sup_df, aes(x=Condition, y=Mean_FilteredNuclei_AreaShape_Area, fill=Condition)) + 
  geom_violin() + geom_boxplot(width=0.1) + scale_fill_brewer(palette="Dark2") +
  ylab("Nuclei Area")


t.test(control_nuclei_area,KIF11_nuclei_area)
t.test(control_nuclei_area,INCENP_nuclei_area)
t.test(control_nuclei_area,KIF11_nuclei_area, alternative = "less")
t.test(control_nuclei_area,INCENP_nuclei_area, alternative = "less")

# there are significant differences between the area of the nuclei in the
# KIF11 and INCENP conditions and the control condition.


# We check if we have the same number of cells on plates 1 and 2.

plate1 <- per_image_data[per_image_data$Image_Metadata_Plate == 1,]
plate2 <- per_image_data[per_image_data$Image_Metadata_Plate == 2,]

cell_count_p1 <- plate1$Image_Count_Cell
cell_count_p2 <- plate2$Image_Count_Cell

boxplot(cell_count_p1,cell_count_p2, names=c("Plate 1", "Plate 2"),
        ylab="Cell count per image")

sup_df_p1 <- data.frame("Image_Count_Cell"=cell_count_p1,
                             "Plate"=rep("Plate 1",length(cell_count_p1)))
sup_df_p2 <- data.frame("Image_Count_Cell"=cell_count_p2,
                           "Plate"=rep("Plate 2",length(cell_count_p2)))
sup_df_plate <- rbind(sup_df_p1,sup_df_p2)
sup_df_plate$Plate <- as.factor(sup_df_plate$Plate)

ggplot(sup_df_plate, aes(x=Plate, y=Image_Count_Cell, fill=Plate)) + 
  geom_violin() + geom_boxplot(width=0.1) + scale_fill_brewer(palette="Dark2") +
  ylab("Cells per image")


t.test(cell_count_p1,cell_count_p2)
t.test(cell_count_p1,cell_count_p2, alternative = "less")

# There are significant differences between the number of cells in plates 1 and 2

# We check if we have the same number of cells on experiments A, B and D

exp_a <- per_image_data[per_image_data$Image_Metadata_Experiment == "A",]
exp_b <- per_image_data[per_image_data$Image_Metadata_Experiment == "B",]
exp_d <- per_image_data[per_image_data$Image_Metadata_Experiment == "D",]

cell_count_exp_a <- exp_a$Image_Count_Cell
cell_count_exp_b <- exp_b$Image_Count_Cell
cell_count_exp_d <- exp_d$Image_Count_Cell

boxplot(cell_count_exp_a,cell_count_exp_b,cell_count_exp_d,
        names = c("Experiment A","Experiment B","Experiment D"),
        ylab = "Cell count per image")

sup_df_a <- data.frame("Image_Count_Cell"=cell_count_exp_a,
                        "Experiment"=rep("A",length(cell_count_exp_a)))
sup_df_b <- data.frame("Image_Count_Cell"=cell_count_exp_b,
                        "Experiment"=rep("B",length(cell_count_exp_b)))
sup_df_d <- data.frame("Image_Count_Cell"=cell_count_exp_d,
                       "Experiment"=rep("D",length(cell_count_exp_d)))

sup_df_experiment <- rbind(sup_df_a,sup_df_b,sup_df_d)
sup_df_experiment$Experiment <- as.factor(sup_df_experiment$Experiment)

ggplot(sup_df_experiment, aes(x=Experiment, y=Image_Count_Cell, fill=Experiment)) + 
  geom_violin() + geom_boxplot(width=0.1) + scale_fill_brewer(palette="Dark2") +
  ylab("Cells per image")


t.test(cell_count_exp_a,cell_count_exp_b)
t.test(cell_count_exp_a,cell_count_exp_d)
t.test(cell_count_exp_b,cell_count_exp_d)

t.test(cell_count_exp_a,cell_count_exp_b, alternative = "less")
t.test(cell_count_exp_a,cell_count_exp_d, alternative = "less")
# There is a significantly lower number of cells in experiment A than in experiments B and D


# Loading per object data

per_object_data <- read.csv("per_object.csv")

# Data exploration

dim(per_object_data)
head(per_object_data[1:10])
str(per_object_data)
summary(per_object_data[1:5])


# We have high dimensional data composed of 599 variables and 57718 cells.
# To handle such a huge dataset we'll use dimensionality reduction (PCA) to compress
# data in a smaller number of dimensions

# Per object data don't have metadata information, so we merge metadata from per
# image data using the column ImageNumber which is common to both data

sup_table <- select(per_image_data,ImageNumber,Image_Metadata_Experiment,Image_Metadata_Condition,
                    Image_Metadata_Plate,Image_Metadata_Gene,Image_Metadata_Replica)

per_object_data <- merge(per_object_data,sup_table,by = "ImageNumber")


# We can see that our variables have very different ranges, for this reason we need to
# scale our data.

# We first remove logical data for PCA

df_pca <- per_object_data[,-which(sapply(per_object_data, class) == "logical")]

sum(is.na(df_pca))

# We have 127 NA values in our dataset
which(colSums(is.na(df_pca)) > 0 )

# We can see that the columns with NA values are related with Neighbors variables, so
# NA values should be for cells that don't have neighbors. We'll remove this cells.

df_pca <- df_pca[-which(rowSums(is.na(df_pca))>0),]

# Check that now we have no NA values
sum(is.na(df_pca))


quali.sup <- 1:3
metadata_cols <- grep("Metadata",colnames(df_pca))
quali.sup <- c(quali.sup,metadata_cols)

# Scale data using scale.unit = TRUE
res_pca <- PCA(df_pca, scale.unit=TRUE, quali.sup = quali.sup, ncp = ncol(df_pca)%/%3,
               graph = FALSE)

fviz_screeplot(res_pca, ncp=20)

head(res_pca$eig)
# Our first 6 components retain 38.5% of variance

# We can plot variables and his relation with the components that explain more variance of data

fviz_pca_var(res_pca) # We see nothing because there are a lot of variables

# Plot individuals

fviz_pca_ind(res_pca, label="none", habillage = as.factor(df_pca$Image_Metadata_Experiment),
             addEllipses=TRUE, ellipse.level=0.95)
fviz_pca_ind(res_pca, label="none", habillage = as.factor(df_pca$Image_Metadata_Condition),
             addEllipses=TRUE, ellipse.level=0.95)
fviz_pca_ind(res_pca, label="none", habillage = as.factor(df_pca$Image_Metadata_Plate),
             addEllipses=TRUE, ellipse.level=0.95)

# We can see that confidence ellipses overlap between cells from different experiments,
# conditions and plates so we don't detect signs of batch effect


# We're retaining the number of components needed to explain the 90% of variance for
# k-Means clustering

eigenvalues <- as.data.frame(res_pca$eig)
eigenvalues[eigenvalues[,3]>90,][1,]
# We retain the first 126 principal components

df_for_km <- as.data.frame(res_pca$ind$coord)[,1:126]

# We estimate the optimal number of clusters for k-Means using Elbow method

set.seed(123)

wss <- 0

# Compute and plot wss for k = 1 to k = 15

for(k in 1:15){
  km <- kmeans(df_for_km,centers = k,nstart=5)
  wss[k] <- km$tot.withinss
}

k.values <- 1:15

plot(k.values, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# It is not very clear the optimal number of clusters. We're going to take 2 and 5 clusters

km_2 <- kmeans(df_for_km, centers = 2, nstart = 5)

fviz_cluster(km_2, data = df_for_km)

km_5 <- kmeans(df_for_km, centers = 5, nstart = 5)

fviz_cluster(km_5, data = df_for_km)

# We can see that clusters have poor quality, it is possible that this is due to the high
# dimensionality of data and the existence of some noisy variables.


# After initial exploration of the complete per object (cells) data with a large number
# of variables, we decided to explore the complete data using machine learning algorithms
# in CellProfiler Analyst, and perform statistical analysis in R of a manageable set of
# variables.

# For detailed analysis of features we select a set features based on three criteria:
# i) fidelity
# ii) biological relevance
# iii) interpretability

df_reduced <- select(per_object_data,ImageNumber,ObjectNumber,
                     Image_Metadata_Experiment,Image_Metadata_Condition,
                     Image_Metadata_Plate,Image_Metadata_Gene,Image_Metadata_Replica,
                     Cell_Neighbors_NumberOfNeighbors_Adjacent,
                     Cell_AreaShape_Area,Cell_AreaShape_Perimeter,
                     Cell_AreaShape_FormFactor,Cell_AreaShape_Eccentricity,
                     Cell_AreaShape_MajorAxisLength,
                     Cell_AreaShape_MinorAxisLength,
                     Cell_AreaShape_Compactness,
                     Cell_Intensity_MeanIntensity_act,
                     Cell_Intensity_MeanIntensity_tub,
                     FilteredNuclei_AreaShape_Area,
                     FilteredNuclei_AreaShape_Perimeter,
                     FilteredNuclei_AreaShape_FormFactor,
                     FilteredNuclei_AreaShape_Eccentricity,
                     FilteredNuclei_Intensity_MeanIntensity_DNA)

sup_variables <- colnames(df_reduced[1:8])

### Data preprocessing

# Check if we have missing values

sum(is.na(df_reduced))

# We transform categorical data in factors

cat_variables <- c("Image_Metadata_Experiment","Image_Metadata_Condition",
                  "Image_Metadata_Plate","Image_Metadata_Gene","Image_Metadata_Replica")

df_reduced[cat_variables] <- lapply(df_reduced[cat_variables], as.factor)


boxplot(df_reduced[!(colnames(df_reduced) %in% sup_variables)])
# We can see that Cell Area is the most variable feature and has a huge range, masking
# the other features boxplots
summary(df_reduced)

# Correlation between variables

cor_matrix <- round(cor(df_reduced[9:ncol(df_reduced)]),3)

corrplot(cor_matrix, type="upper", diag = F, 
         tl.col="black", tl.srt=45)


# Distribution of variables

hist.data.frame(df_reduced[9:ncol(df_reduced)])

# Our features have a different range of values so we need to standarize data.
# We are using robust standarization that centers around the median using MAD, and are robust
# to outliers


# Data are standarized independently by experiment and by plate and then merged

st_a_1 <- df_reduced[df_reduced$Image_Metadata_Experiment=="A" & df_reduced$Image_Metadata_Plate==1,] 
st_a_1 <- standardize(st_a_1[9:ncol(st_a_1)],robust = TRUE)
st_a_2 <- df_reduced[df_reduced$Image_Metadata_Experiment=="A" & df_reduced$Image_Metadata_Plate==2,]
st_a_2 <- standardize(st_a_2[9:ncol(st_a_2)],robust = TRUE)
st_b_1 <- df_reduced[df_reduced$Image_Metadata_Experiment=="B" & df_reduced$Image_Metadata_Plate==1,]
st_b_1 <- standardize(st_b_1[9:ncol(st_b_1)],robust = TRUE)
st_b_2 <- df_reduced[df_reduced$Image_Metadata_Experiment=="B" & df_reduced$Image_Metadata_Plate==2,]
st_b_2 <- standardize(st_b_2[9:ncol(st_b_2)],robust = TRUE)
st_d_1 <- df_reduced[df_reduced$Image_Metadata_Experiment=="D" & df_reduced$Image_Metadata_Plate==1,]
st_d_1 <- standardize(st_d_1[9:ncol(st_d_1)],robust = TRUE)
st_d_2 <- df_reduced[df_reduced$Image_Metadata_Experiment=="D" & df_reduced$Image_Metadata_Plate==2,]
st_d_2 <- standardize(st_d_2[9:ncol(st_d_2)],robust = TRUE)

df_stand <- rbind(st_a_1,st_a_2,st_b_1,st_b_2,st_d_1,st_d_2)

# Check that data are centered around the median
mean(st_a_1$Cell_AreaShape_Area)
median(st_a_1$Cell_AreaShape_Area)

summary(df_stand)

hist.data.frame(df_stand)


# Distribution of features is more homogeneous now

# We can't check if data fit normal distribution using Shapiro Test because
# p-value calculation can't be trusted for sample sizes above 5000

shapiro.test(df_stand$Cell_AreaShape_FormFactor)


# Group by gene using mean

df_stand <- cbind(df_stand, select(df_reduced,Image_Metadata_Gene))

gene_group_df <- df_stand %>% 
                     group_by(Image_Metadata_Gene) %>%
                     summarise_all("mean")
rownames(gene_group_df) <- gene_group_df$Image_Metadata_Gene

# Heatmap

gene_group_mat <- as.matrix(gene_group_df[-1])
row.names(gene_group_mat) <- gene_group_df[["Image_Metadata_Gene"]]


heatmap(gene_group_mat)
heatmap.2(gene_group_mat, trace = "none", scale = "column")


# Plotting Z-scores of different features to identify possible hits

# Since we have standardized the data from different plates and experiments independently,
# we have to subtract the median from the total set to center the z-scores around the median.

medians <- apply(gene_group_df[-1], 2, median)

for (i in 2:ncol(gene_group_df)) {
  gene_group_df[i] <- gene_group_df[i]-medians[i-1]
}


plot_hits <- function(feature){
  ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=feature)) + 
    geom_point() + ylim(-5,5) +
    geom_text(aes(label=ifelse(feature > 1 | feature < -1,
                               as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2)
  
}

### Initial representation
# Cell Area

plot_hits(gene_group_df$Cell_AreaShape_Area)

# Cell Perimeter

plot_hits(gene_group_df$Cell_AreaShape_Perimeter)

# Cell Form Factor

plot_hits(gene_group_df$Cell_AreaShape_FormFactor)

# Eccentricity

plot_hits(gene_group_df$Cell_AreaShape_Eccentricity)

# Cell Major Axis Length

plot_hits(gene_group_df$Cell_AreaShape_MajorAxisLength)

# Cell Minor Axis Length

plot_hits(gene_group_df$Cell_AreaShape_MinorAxisLength)

# Cell Compactness

plot_hits(gene_group_df$Cell_AreaShape_Compactness)

# Tubulin Intensity

plot_hits(gene_group_df$Cell_Intensity_MeanIntensity_tub)

# Actin intensity

plot_hits(gene_group_df$Cell_Intensity_MeanIntensity_act)

# DAPI Intensity

plot_hits(gene_group_df$FilteredNuclei_Intensity_MeanIntensity_DNA)

# Nuclei Area

plot_hits(gene_group_df$FilteredNuclei_AreaShape_Area)

# Nuclei Perimeter

plot_hits(gene_group_df$FilteredNuclei_AreaShape_Perimeter)

# Nuclei Form Factor

plot_hits(gene_group_df$FilteredNuclei_AreaShape_FormFactor)

# Nuclei Eccentricity

plot_hits(gene_group_df$FilteredNuclei_AreaShape_Eccentricity)


### Detailed representation for manuscript

# Cell Area
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=Cell_AreaShape_Area)) + 
  geom_point() + ylim(-3,3) +
  geom_text(aes(label=ifelse(Cell_AreaShape_Area > 1 | Cell_AreaShape_Area < -1,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2, size=4.5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Cell Area Z-Scores") +
  geom_hline(yintercept=1,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-1,linetype="dashed", 
             color = "red", size=1)



# Cell Perimeter
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=Cell_AreaShape_Perimeter)) + 
  geom_point() + ylim(-3,3) +
  geom_text(aes(label=ifelse(Cell_AreaShape_Perimeter > 1 | Cell_AreaShape_Perimeter < -1,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2, size=4.5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Cell Perimeter Z-Scores") + 
  geom_hline(yintercept=1,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-1,linetype="dashed", 
             color = "red", size=1)

# Nuclei Area
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=FilteredNuclei_AreaShape_Area)) + 
  geom_point() + ylim(-3,3) +
  geom_text(aes(label=ifelse(FilteredNuclei_AreaShape_Area > 1 | FilteredNuclei_AreaShape_Area < -1,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2, size=4.5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Nuclei Area Z-Scores") +
  geom_hline(yintercept=1,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-1,linetype="dashed", 
             color = "red", size=1)  

# Nuclei Perimeter
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=FilteredNuclei_AreaShape_Perimeter)) + 
  geom_point() + ylim(-3,3) +
  geom_text(aes(label=ifelse(FilteredNuclei_AreaShape_Perimeter > 1 | FilteredNuclei_AreaShape_Perimeter < -1,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2,size=4.5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Nuclei Perimeter Z-Scores") +
  geom_hline(yintercept=1,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-1,linetype="dashed", 
             color = "red", size=1) 


# Cell Major Axis Length
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=Cell_AreaShape_MajorAxisLength)) + 
  geom_point() + ylim(-3,3) +
  geom_text(aes(label=ifelse(Cell_AreaShape_MajorAxisLength > 1 | Cell_AreaShape_MajorAxisLength < -1,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2, size=4.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Cell Major Axis Length Z-Scores") +
  geom_hline(yintercept=1,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-1,linetype="dashed", 
             color = "red", size=1)

# Cell Minor Axis Length
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=Cell_AreaShape_MinorAxisLength)) + 
  geom_point() + ylim(-3,3) +
  geom_text(aes(label=ifelse(Cell_AreaShape_MinorAxisLength > 1 | Cell_AreaShape_MinorAxisLength < -1,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2, size=4.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Cell Minor Axis Length Z-Scores") +
  geom_hline(yintercept=1,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-1,linetype="dashed", 
             color = "red", size=1)



# Cell FormFactor

ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=Cell_AreaShape_FormFactor)) + 
  geom_point() + ylim(-2,2) +
  geom_text(aes(label=ifelse(Cell_AreaShape_FormFactor > 0.75 | Cell_AreaShape_FormFactor < -0.75,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.2,vjust=0.2, size=4.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Cell Form Factor Z-Scores") +
  geom_hline(yintercept=0.75,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-0.75,linetype="dashed", 
             color = "red", size=1)

# Cell Compactness
ggplot(gene_group_df, aes(x=Image_Metadata_Gene, y=Cell_AreaShape_Compactness)) + 
  geom_point() + ylim(-2,2) +
  geom_text(aes(label=ifelse(Cell_AreaShape_Compactness > 0.75 | Cell_AreaShape_Compactness < -0.75,
                             as.character(Image_Metadata_Gene),'')),hjust=-0.1,vjust=0.2, size=4.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Cell Compactness Z-Scores") +
  geom_hline(yintercept=0.75,linetype="dashed", 
             color = "red", size=1) +
  geom_hline(yintercept=-0.75,linetype="dashed", 
             color = "red", size=1)


