# Instructions pour le développement

Depuis la racine du dépôt, lire [les contrats d’architecture](docs/architecture.md) et [la feuille de route de Luna](docs/LUNA.md). Le [cahier des charges](CdC_noyau_trivarie_sans_trim_corrige.md) définit les garanties scientifiques. Les notes de discussion historiques conservent le contexte des choix initiaux.

Lorsque l’outil de graphe est disponible, privilégier son index du dépôt, la recherche de symboles et le suivi des appels. En cas de résultat insuffisant, utiliser la recherche de texte et les fichiers sources. Les commandes du projet s’exécutent depuis un clone standard. Sur les postes équipés de RTK, appliquer aussi les conventions RTK locales.

Le noyau, les constructeurs et les audits de production sont développés en Rust. Python sert aux oracles indépendants et aux références historiques. Chaque lot géométrique livre du Rust accessible depuis la crate et testé avec Cargo.

1. Découper par capacité réutilisable. Pour chaque lot, consigner la révision de départ, les fichiers concernés, les préconditions, les seuils, l’oracle indépendant, les entrées à refuser et les critères de sortie avant de développer.
2. Écrire le code, les docstrings, les commentaires techniques et les commits en anglais. Préserver les recettes et les formats versionnés existants. Une extension de famille possède un identifiant et des conventions explicites.
3. Garder la bibliothèque indépendante de l’interface utilisateur. Le CLI et les futurs adaptateurs passent par le même constructeur. Toute adoption d’un modèle passe par l’audit et un remplacement atomique. Tester la conservation de l’ancien état après refus.
4. Distinguer fidélité de représentation, régularité locale, injectivité, compatibilité des traces, couverture et aptitude à l’analyse. Un échantillonnage fournit un diagnostic. Un certificat nomme son domaine et sa méthode. Une preuve incomplète donne un résultat indéterminé.
5. Développer les oracles indépendamment du constructeur. Ranger les références analytiques dans tests/references, les entrées invalides dans tests/counterexamples et chaque défaut numérique reproduit dans tests/regressions avec sa recette et l’assertion qui échoue sur l’ancienne révision.
6. Rejouer `make check` à l’intégration. Rapporter révision, recette, environnement, commandes, mesures, couverture et limites. Compter un refus sur un cas obligatoire comme échec de couverture. Garder build/ pour les résultats régénérables.
7. upstream/ contient des références en lecture seule. Documenter origine, révision, licence, avis à conserver et contrat numérique avant de reprendre du code. Le choix de licence du code original appartient au porteur du projet avant diffusion (NF-8).
8. En cas de conflit entre une construction et son oracle, conserver le contre-exemple et identifier l’hypothèse en cause avant de modifier les seuils ou le périmètre.
