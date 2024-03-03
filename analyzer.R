data <- read.csv('C:/Users/Valentin/PycharmProjects/Projet-antibio/table.csv')

nombre_antibiotiques <- ncol(data) - 3
#167 antibio test??

nombre_especes_uniques <- length(unique(data$taxon_id))
#3525

comptage_genomes <- aggregate(genome_id ~ taxon_id, data = data, FUN = function(x) length(unique(x)))

names(comptage_genomes) <- c("taxon_id", "nombre_genomes")

library(ggplot2)

comptage_genomes_filtre <- comptage_genomes[comptage_genomes$nombre_genomes > 1,]

ggplot(comptage_genomes_filtre, aes(x = taxon_id, y = nombre_genomes)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Nombre de G??nomes par Taxon ID (Filtr??)",
       x = "Taxon ID",
       y = "Nombre de G??nomes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 


library(plotly)


# Convertir le ggplot en graphique plotly interactif
p <- ggplot(comptage_genomes, aes(x = taxon_id, y = nombre_genomes)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Nombre de G??nomes par Taxon ID (Interactif)",
       x = "Taxon ID",
       y = "Nombre de G??nomes")

ggplotly(p)

comptage_genomes_filtre$taxon_id <- factor(_filtre$taxon_id)

ggplot(comptage_genomes_filtre, aes(x = taxon_id, y = nombre_genomes)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Nombre de G??nomes par Taxon ID (Filtr??)",
       x = "Taxon ID",
       y = "Nombre de G??nomes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
