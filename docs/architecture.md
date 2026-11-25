# Architecture et contrats

## Noyau natif Rust

La crate `solidcraft-trivariate` expose `RationalBezierMap::new` et `evaluate` depuis `src/lib.rs`. Toutes les nouvelles fonctionnalités de production passent par cette crate. Les modèles et les tests Python acquis pendant les premiers lots fournissent des références pour qualifier les portages.

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
