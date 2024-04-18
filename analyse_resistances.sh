
#!/bin/bash

# Répertoire contenant les fichiers genome.fasta à traiter
repertoire="/genomes/"

# Liste des bases de données à utiliser
databases=("ncbi" "resfinder" "card" "argannot")

# Redirection de la sortie vers test.txt
> results_all_genomes_all_db.txt

# Parcourir chaque fichier genome.fasta dans le répertoire
for file in "$repertoire"*.fasta; do
    echo "Traitement du fichier : $file"
    # Parcourir chaque base de données
    for db in "${databases[@]}"; do
        echo "Analyse avec la base de données $db :"
        /home/anais/miniconda3/bin/abricate --db "$db" "$file" >> results_all_genomes_all_db.txt
    done
    echo "--------------------------------------------------" >> results_all_genomes_all_db.txt
done

