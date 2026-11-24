# Instructions pour le développement

Sur cette machine, suivre @/home/jeff/.codex/RTK.md et préfixer les commandes shell par `rtk`. Pour découvrir le code, privilégier le graphe MCP codebase-memory après indexation du projet. Utiliser search_graph, trace_path et get_code_snippet, puis query_graph. Réserver la recherche texte aux documents, aux littéraux et aux résultats insuffisants.

Lire docs/architecture.md et docs/LUNA.md. Le CdC corrigé définit les garanties scientifiques. Les deux notes d’adaptation précisent la progression du démonstrateur.

Le noyau, les constructeurs et les audits de production sont développés en Rust. Python sert aux oracles indépendants et aux références historiques. Chaque lot géométrique livre du Rust accessible depuis la crate et testé avec Cargo.

1. Découper par capacité réutilisable. Pour chaque lot, consigner la révision de départ, les fichiers concernés, les préconditions, les seuils, l’oracle indépendant, les entrées à refuser et les critères de sortie avant de développer.
2. Écrire le code, les docstrings, les commentaires techniques et les commits en anglais. Préserver les recettes et les formats versionnés existants. Une extension de famille possède un identifiant et des conventions explicites.
3. Garder la bibliothèque indépendante de l’interface utilisateur. Le CLI et les futurs adaptateurs passent par le même constructeur. Toute adoption d’un modèle passe par l’audit et un remplacement atomique. Tester la conservation de l’ancien état après refus.
4. Distinguer fidélité de représentation, régularité locale, injectivité, compatibilité des traces, couverture et aptitude à l’analyse. Un échantillonnage fournit un diagnostic. Un certificat nomme son domaine et sa méthode. Une preuve incomplète donne un résultat indéterminé.
5. Développer les oracles indépendamment du constructeur. Ranger les références analytiques dans tests/references, les entrées invalides dans tests/counterexamples et chaque défaut numérique reproduit dans tests/regressions avec sa recette et l’assertion qui échoue sur l’ancienne révision.
6. Rejouer `make check` à l’intégration. Rapporter révision, recette, environnement, commandes, mesures, couverture et limites. Compter un refus sur un cas obligatoire comme échec de couverture. Garder build/ pour les résultats régénérables.
7. upstream/ contient des références en lecture seule. Documenter origine, révision, licence, avis à conserver et contrat numérique avant de reprendre du code. Le choix de licence du code original appartient au porteur du projet avant diffusion (NF-8).
8. En cas de conflit entre une construction et son oracle, conserver le contre-exemple et identifier l’hypothèse en cause avant de modifier les seuils ou le périmètre.
