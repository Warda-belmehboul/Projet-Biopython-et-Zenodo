#!/usr/bin/env python3
"""
SCRIPT PRINCIPAL D'EXÉCUTION
Université de Jijel - Master I Sciences de la Nature et de la Vie
Travail Pratique : Logiciels Libres et Open Source
"""

import sys
import os
from datetime import datetime

def main():
    """Fonction principale d'exécution"""
    
    print("\n" + "=" * 80)
    print("UNIVERSITÉ DE JIJEL - FACULTÉ DES SCIENCES DE LA NATURE ET DE LA VIE")
    print("MASTER I - TRAVAIL PRATIQUE : LOGICIELS LIBRES ET OPEN SOURCE")
    print("=" * 80)
    print(f"Date d'exécution : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Étudiant : [Votre Nom]")
    print(f"Groupe : [Numéro du Groupe]")
    print("=" * 80)
    
    # Vérification des prérequis
    print("\n[VÉRIFICATION DES PRÉREQUIS]")
    print("-" * 40)
    
    try:
        import Bio
        print(f"✓ Biopython version {Bio.__version__} installé")
    except ImportError:
        print("✗ Biopython non installé. Installation : pip install biopython")
        sys.exit(1)
    
    try:
        import pandas as pd
        print(f"✓ Pandas version {pd.__version__} installé")
    except ImportError:
        print("✗ Pandas non installé")
    
    try:
        import matplotlib
        print(f"✓ Matplotlib version {matplotlib.__version__} installé")
    except ImportError:
        print("✗ Matplotlib non installé")
    
    # Vérification du fichier XML
    xml_file = "4967697.xml"
    if os.path.exists(xml_file):
        print(f"✓ Fichier XML trouvé : {xml_file}")
        file_size = os.path.getsize(xml_file) / 1024
        print(f"  Taille : {file_size:.1f} KB")
    else:
        print(f"✗ Fichier XML non trouvé : {xml_file}")
        print("  Veuillez placer le fichier 4967697.xml dans le même répertoire")
        sys.exit(1)
    
    # Import et exécution de l'analyse
    print("\n" + "=" * 80)
    print("LANCEMENT DE L'ANALYSE")
    print("=" * 80)
    
    try:
        # Import dynamique des fonctions
        from analyse_metadonnees_zenodo import (
            analyser_metadonnees_dublin_core,
            analyser_sequences_transposons,
            analyser_alignements_transposons,
            rechercher_sequences_ncbi,
            simulation_analyse_blast,
            generer_rapport_complet
        )
        
        # Exécution complète
        generer_rapport_complet()
        
        print("\n" + "=" * 80)
        print("RÉSUMÉ DE L'EXÉCUTION")
        print("=" * 80)
        
        # Statistiques d'exécution
        fichiers_generes = [
            "pgbd5_transposon.fasta",
            "similarite_transposons.png",
            "rapport_analyse_zenodo_biopython.txt"
        ]
        
        print("\nFichiers générés :")
        for fichier in fichiers_generes:
            if os.path.exists(fichier):
                size_kb = os.path.getsize(fichier) / 1024
                print(f"  ✓ {fichier} ({size_kb:.1f} KB)")
            else:
                print(f"  ✗ {fichier} (non généré)")
        
        print("\nCompétences démontrées :")
        competences = [
            "Analyse XML Dublin Core",
            "Manipulation de séquences ADN",
            "Alignements et similarités",
            "Recherche dans bases de données",
            "Visualisation scientifique",
            "Génération de rapports"
        ]
        
        for i, competence in enumerate(competences, 1):
            print(f"  {i}. {competence}")
        
        print("\n" + "=" * 80)
        print("TRAVAIL TERMINÉ AVEC SUCCÈS")
        print("=" * 80)
        
    except Exception as e:
        print(f"\n[ERREUR] Problème lors de l'exécution : {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
