#!/bin/bash

# Vérifie qu'un fichier contenant la liste des répertoires est passé en argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <fichier_liste_repertoires>"
    # exit 1
fi

# Lire la liste des répertoires
file_list=$1

# Boucle sur chaque répertoire mentionné dans le fichier
while IFS= read -r base_dir; do
    echo "Traitement du répertoire : $base_dir"
    
    # Vérifie que le répertoire existe
    if [ ! -d "$base_dir" ]; then
        echo "Erreur : Le répertoire $base_dir n'existe pas."
        continue
    fi

    # Liste les sous-répertoires triés par ordre alphabétique
    subdirs=( "$base_dir"/000* )
    subdir_count=${#subdirs[@]}
    
    # Si plus d'un sous-répertoire existe
    if [ "$subdir_count" -gt 1 ]; then
        echo "  - $subdir_count sous-répertoires trouvés."

        # Le répertoire cible est toujours "0000"
        target_dir="$base_dir/0000"
        if [ ! -d "$target_dir" ]; then
            echo "  - Le répertoire 0000 est manquant. Création en cours."
            mkdir -p "$target_dir"
        fi

        # Déplace les fichiers des répertoires "0001" et suivants vers "0000"
        for subdir in "${subdirs[@]:1}"; do
            echo "  - Déplacement des fichiers depuis $subdir vers $target_dir"
            mv "$subdir"/* "$target_dir/" 2>/dev/null
        done

        echo "  - Suppression des répertoires vides restants."
        for subdir in "${subdirs[@]:1}"; do
            rmdir "$subdir" 2>/dev/null || echo "    Impossible de supprimer $subdir (non vide ?)"
        done
    else
        echo "  - Aucun traitement nécessaire (moins de 2 sous-répertoires)."
    fi
done < "$file_list"

echo "Traitement terminé."
