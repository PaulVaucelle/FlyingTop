#!/bin/bash

# Vérifie qu'au moins deux fichiers sont passés en argument
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fichier_liste_repertoires> <fichier_noms_output>"
    exit 1
fi

# Fichiers d'entrée
file_list=$1
output_list=$2

# Vérifie que les fichiers existent
if [ ! -f "$file_list" ] || [ ! -f "$output_list" ]; then
    echo "Erreur : Les fichiers $file_list et/ou $output_list n'existent pas."
    exit 1
fi

# Lire les fichiers ligne par ligne et associer les répertoires aux fichiers de sortie
paste "$file_list" "$output_list" | while IFS=$'\t' read -r dir output_file; do
    echo "Traitement du répertoire : $dir"
    echo "  - Fichier de sortie : $output_file"
    
    # Vérifie que le répertoire existe
    if [ ! -d "$dir" ]; then
        echo "Erreur : Le répertoire $dir n'existe pas."
        continue
    fi

    # Se déplacer dans le répertoire
    cd "$dir" || { echo "Impossible de se déplacer dans $dir"; continue; }
    find . -type f -empty -delete
    # Cherche tous les fichiers ROOT
    root_files=(*.root)
    if [ ${#root_files[@]} -eq 0 ]; then
        echo "  - Aucun fichier ROOT trouvé dans $dir."
        cd - > /dev/null
        continue
    fi

    # Fusionner les fichiers ROOT avec hadd
    echo "  - Fusion des fichiers ROOT dans $output_file"
    hadd -ff "$output_file" *.root


    # Si la fusion a réussi, supprimer tous les fichiers sauf l'output
    if [ $? -eq 0 ]; then
        echo "  - Suppression des fichiers ROOT sauf $output_file."
        for file in *.root; do
            if [ "$file" != "$output_file" ]; then
                rm -f "$file"
            fi
        done
    else
        echo "  - La fusion a échoué, les fichiers ROOT ne seront pas supprimés."
    fi
    
    # Retour au répertoire précédent
    cd - > /dev/null
done

echo "Traitement terminé."
