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
