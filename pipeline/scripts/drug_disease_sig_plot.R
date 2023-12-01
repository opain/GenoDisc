# Make a plot showing disease and drug signitures
set.seed(1)
n_gene<-1000
n_drug<-5
genes<-paste0('gene_',1:n_gene)
drugs<-paste0('drug_',1:n_drug)
drug <- matrix(
  rbinom(n_gene * n_drug, 2, 0.5) - 1,
  nrow = n_drug,
  ncol = n_gene)

disease<-rnorm(n_gene)

signitures<-data.frame(
  name=c('disease', drugs)
)

signitures<-data.table(cbind(
  signitures,
  rbind(
    disease,
    drug
  )
))

names(signitures)<-c('name',genes)

# Retain genes that are strongly associated with disease
signitures<-signitures[,c('name',genes[abs(disease) > 1.96]), with=F]

signitures_melt<-melt(signitures)

genes[order(disease)]

signitures_melt$variable<-factor(signitures_melt$variable, levels=genes[order(disease)])

library(ggplot2)
library(cowplot)

ggplot(signitures_melt, aes(x=variable, y=value, group=name, colour=name)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_line() +
  facet_grid(name ~ .) +
  theme_half_open() +
  theme(axis.text.x = element_blank()) + 
  panel_border() +
  labs(x = 'Genes', y = 'Z-score')
  
############################################
# Simulate drug gene 
############################################

set.seed(123)  # For reproducibility

# Simulate data
n_genes = 100  # Number of genes
disease_z_scores = rnorm(n_genes)  # Disease-associated Z scores
tissue_z_scores = rnorm(n_genes)  # Tissue-specific Z scores
noise = rnorm(n_genes, sd = 0.5)  # Random noise

# Simulate drug Z scores with negative association
drug_z_scores = -disease_z_scores * tissue_z_scores + noise

# Combine into a data frame
data = data.frame(
  Gene = paste("Gene", 1:n_genes),
  Disease_Z = disease_z_scores,
  Drug_Z = drug_z_scores,
  Tissue_Z = tissue_z_scores
)

# View the first few rows of the data
data_melt<-melt(data)

ggplot(data_melt, aes(x=Gene, y=value, group=variable, colour=variable)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_line() +
  facet_grid(variable ~ .) +
  theme_half_open() +
  theme(axis.text.x = element_blank()) + 
  panel_border() +
  labs(x = 'Genes', y = 'Z-score')

