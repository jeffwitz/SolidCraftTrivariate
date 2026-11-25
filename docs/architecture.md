# Architecture et contrats

## Noyau natif Rust

La crate `solidcraft-trivariate` expose `RationalBezierMap::new` et `evaluate` depuis `src/lib.rs`. Toutes les nouvelles fonctionnalités de production passent par cette crate. Les modèles et les tests Python acquis pendant les premiers lots fournissent des références pour qualifier les portages.

## Planificateur

Le chemin polygonal `sketch_layer_partition_v1` prend un contour de stock, des contours de coupe, les intervalles axiaux matière/vide et les guides explicites de l’esquisse. Dans la stratégie automatique, le constructeur compare des propositions de séparatrices. Il construit les cartes de chaque proposition admise et vérifie l’aire matérielle exacte de chaque couche, la positivité du jacobien sur le fermé et les incidences des faces entières. Il sélectionne ensuite le plus petit nombre de blocs. À égalité, la partition historique conserve la priorité. Le plan sérialisé permet de retrouver la méthode, les segments avec leur provenance et les effectifs des couches. [Construction](../src/sketch_layers.rs) · [Lecture et rejeu du plan](../src/sketch_layers_json.rs)

Une esquisse déjà quadrangulée peut garder ses faces convexes si toutes leurs arêtes coïncident entre voisines. Le plan de la poche relie les quatre coins de l’empreinte à ceux du stock. Sur certains contours simples à nombre pair de sommets, un deuxième procédé détache des quadrilatères convexes le long de diagonales visibles. Il rattache ensuite chaque trou quadrangulaire aux quatre coins de sa cellule d’accueil. Les séparatrices s’arrêtent sur les sommets correspondants. Le procédé est déterministe et peut échouer à trouver une décomposition existante. La grille orthogonale et la triangulation suivie de quadrangulation restent disponibles pour les autres configurations. [Arrangement](../src/sketch_arrangement.rs) · [Conformité des quadrilatères](../src/sketch_quads.rs) · [Diagonales locales](../src/sketch_quad_separators.rs)

Les identifiants de la poche décrivent la paroi et la couche par rapport au fond. Dans les autres partitions guidées, une cellule est identifiée par ses coins rationnels canoniques et son intervalle axial. Ces clés protègent les cellules inchangées contre une renumérotation globale. Une modification de coordonnées peut changer les clés géométriques des cellules incidentes. L’optimisation de la correspondance entre deux plans de topologies différentes et la sélection générale du nombre de cellules incidentes à chaque coin restent ouvertes.

Les opérations à axes croisés et les contours NURBS suivent des chemins dédiés. `import-solidcraft-sketch` lit l’historique de features de l’équerre et construit le modèle depuis ses esquisses. `plan-sketch` explicite les événements d’intersection des profils NURBS pris en charge. Les gabarits spécialisés et les modèles mécaniques fournissent des références de forme, de coût et de conditionnement. La comparaison de la poche et du C utilise des raccords face entière et une partition commune aux couches. Le cas à cinq volumes proposé avec jonctions en T requiert un autre contrat de raccord pour cette famille. [Comparaison et bornes](partition-poche-comparaison.md) · [Chemins de démonstration](../README.md)

### Objectif de sélection et d’édition

Le contrat cible prend comme entrée les contours, le constructeur 3D (extrusion, révolution…) et les traitements prévus. Parmi les candidats qui satisfont la fidélité de forme, la couverture, les raccords et la régularité sur le domaine fermé, il privilégie le plus petit nombre de volumes logiques. Pour une édition, il essaie d’abord la connectivité existante. À nombre de volumes égal, il privilégie la stabilité des identités, puis compare les travées, les points de contrôle et le coût d’analyse. Une extraction de Bézier pour certifier ou intégrer une travée ne crée aucune nouvelle identité de volume de conception.

Chaque gabarit doit annoncer sa classe de topologies, ses hypothèses de régularité et ses règles de placement des nœuds. Pour la révolution d’un profil, le planificateur examine où les coins paramétriques des cartes aboutissent sur le contour après édition. Deux faces extérieures d’une même cellule tensorielle qui se rejoignent sur une peau devenue lisse imposent un déterminant nul sur leur arête commune. Cette obstruction écarte ce candidat précis. Dans ce cas le plan doit documenter la raison du refus avant d’envisager un changement de topologie. Le rapport distingue un candidat admissible, le meilleur candidat testé et un minimum démontré. Cette règle guide le développement du planificateur.

## Modèle historique des références Python

La recette versionnée est l’intention de construction. `trivariate.build.build(recipe)` produit un modèle candidat en mémoire. `trivariate.audit.audit(model)` prend uniquement ce modèle et reconstruit les coefficients attendus à partir des dimensions déclarées. Le CLI sérialise en JSON après validation, dans un fichier temporaire du même répertoire que la destination. Le renommage est atomique. En cas de deux publications concurrentes valides, la dernière remplace la précédente.

La famille fermée `affine_pair_v1` encode deux cartes sur [0,1]³, en millimètres. Chaque axe est de degré 1, avec nœuds [0,0,1,1]. Les huit coefficients (wx,wy,wz,w) sont rangés avec k variant le plus vite, puis j, puis i. Tous les poids valent 1. Les identités left et right proviennent des deux positions du gabarit. Le côté u=1 de left et le côté u=0 de right référencent une même interface orientée par (v,w).

Pour les largeurs entières positives a,b, les cartes sont (au,v,w) et (a+bu,v,w). L’égalité de tous les coefficients avec ces cartes garantit le déterminant positif sur le fermé, l’injectivité par bloc, l’égalité des traces et une partition de [0,a+b] × [0,1]². Le calcul entier du prototype rend cette vérification exacte dans cette famille. L’oracle renvoie valid, invalid ou indeterminate pour un schéma inconnu.

## Extension prévue

| Couche | Données et responsabilité |
| --- | --- |
| BlockPlan | Features résolues, paramètres avec unités, provenance et incidences attendues. |
| TrivariateComplex | Cartes rationnelles natives, orientations et frontière. BlockSide désigne un côté logique, TopologicalFace une incidence partagée, DesignFace l’identité sélectionnée par l’utilisateur. |
| Vérification | Preuves par bloc puis par interface et pour la partition complète. Consigner budget, motifs de refus et bornes. |
| AnalysisSpace | Bases, raffinement, poids de trace, assemblage H¹, contraintes et quadrature à partir de la géométrie acceptée. |
| Adaptateurs | SolidCraft fournit l’intention des features. FieldCraft fournit le corpus et la comparaison EF. Le rendu échantillonne les cartes. |

Le premier raccord courbe porte sur des faces entières. Les côtés ont une orientation explicite, les degrés, nœuds, multiplicités, poids et coefficients homogènes doivent coïncider après réorientation et normalisation commune. Les raccords emboîtés et la frontière partielle relèvent d’une qualification ultérieure.

## Références locales

- SolidCraft et son sous-arbre solidcore fournissent la façade de corps, l’évaluation des features et les règles de coordination. FieldCraft fournit `tests/corpus/filleted-bracket-blocks/{recipe,fixture}.json` et `crates/testkit/tests/bracket_blocks.rs`.
- AutoCraft, craft-contrib, Craft Anything et ArtCraft documentent les revues, les tests, les commandes, la provenance et les instructions communes. Les clones de Craft Anything et ArtCraft utilisent un checkout partiel pour les chemins étudiés. Les six révisions figurent dans references.lock.
- Les références de test utilisent la bibliothèque standard Python. Les applications Rust amont disposent de leurs propres dépendances et cycles de compilation. Leur compilation interviendra lors du travail sur les adaptateurs.


## Parcours effectivement implémenté en Rust

Le CLI lit des recettes JSON versionnées. Le point d’entrée `interchange::build` route chaque famille vers son constructeur puis son audit. L’audit réexamine les coefficients homogènes, le domaine des paramètres et les incidences propres à cette famille. L’export STEP assemble les faces orientées de la frontière et leurs composantes, y compris les coquilles internes des cavités. [Routage](../src/interchange.rs) · [Frontière](../src/boundary_shell.rs)

Les charges et appuis utilisent les identités des faces de conception. L’espace Q2 conserve sa propre numérotation et ses contraintes de trace. Le plan dyadique appartient au job de calcul. L’opérateur de relèvement élimine les équations dépendantes avant la résolution. Le résultat conserve le plan, les blocs, les charges nodales, les déplacements, les réactions et le bilan des travaux. [Jobs mécaniques](../src/mechanical_job.rs) · [Contraintes](../src/constraint_system.rs)

Le job v5 permet de contrôler la quadrature des faces rationnelles. Avec `derivative: "freeze_primal"`, la dérivée réemploie les subdivisions primales de chaque face et compte les panneaux supplémentaires dans le budget. Le job v6 accepte une force volumique spatialement affine et une règle `body_quadrature` avec ordres minimum et maximum, tolérance relative et plafond de points. L’estimation compare la charge nodale complète, le volume, la résultante et le moment entre deux ordres. Les comptes rendent visibles les ordres rejetés et le coût de sélection. Le seuil de convergence entre règles reste un indicateur numérique. [Surface](../src/adaptive_traction.rs) · [Volume](../src/volume_load.rs)

La session géométrique d’édition d’un intervalle conserve les couches dont la géométrie et le masque de matière restent compatibles. La session mécanique de l’équerre mémorise les contributions locales avec les contrôles homogènes, les paramètres d’intégration et les caractéristiques du matériau. Ces deux sessions publient un nouvel état après les contrôles d’admission. Le bilan indique pour chaque expérience le périmètre du cache, les vérifications froides et les mesures de temps. [Édition géométrique](../src/interval_bores_session.rs) · [Édition mécanique](../src/bracket_mechanics_session.rs) · [Bilan](bilan-campagne-20261125.md)
