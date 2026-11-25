# Partition locale des équerres FieldCraft

La référence de départ est la révision `aa94b1b`. Ses cinq équerres ont chacune 48 blocs et 74 interfaces. Sur chaque aile, la face extérieure percée porte 24 panneaux.

## Critères de qualification

- Reconstruire les cinq solides depuis les mêmes recettes. Conserver le contour en L, le congé R6, les positions et les diamètres des deux trous. Préserver la lecture et l’audit des anciens modèles.
- Sur un prisme rectangulaire percé, relier chaque côté extérieur directement à un quart de cercle. L’objectif est de quatre blocs traversants. Sur l’équerre, essayer quatre blocs dans chaque aile et deux dans le coude.
- Vérifier sur le cube fermé la positivité du déterminant de chaque carte. Vérifier les raccords en coordonnées homogènes, les parois cylindriques et les empreintes planes.
- Comparer chaque STEP avec le solide historique correspondant par intégration adaptative et par classification de matière, y compris à ±0,02 mm de chaque paroi. Les seuils géométriques sont ceux du rejeu précédent.
- Compter les panneaux des faces extérieures percées directement sur les cartes. La réduction du nombre de blocs ne constitue pas une preuve de minimalité sur l’ensemble des partitions 3D possibles.

## Résultat géométrique

| Mesure par équerre | Grille précédente | Partition locale |
|---|---:|---:|
| Blocs trivariés | 48 | 10 |
| Raccords entre blocs | 74 | 11 |
| Faces extérieures de blocs | 140 | 38 |
| Panneaux sur chacune des faces planes percées x = 0 et y = 0 | 24 | 5 |

Les cinq variantes partagent cette topologie. Les quatre cartes de chaque aile vont directement du cercle au bord de la partie rectangulaire de l’aile. Le cinquième panneau sur la grande face appartient au coude. Les deux cartes du coude sont communes aux deux ailes. Les onze raccords comparent les réseaux homogènes, après élévation de degré si nécessaire. Le prisme rectangulaire à un seul trou comporte quatre blocs et quatre raccords, pour chacun des axes d’extrusion X, Y et Z testés. Les profils en L réfléchis et translatés ont aussi été reconstruits dans les quatre orientations du plan.

Les huit faces cylindriques de perçage sont exportées sur des supports STEP `CYLINDRICAL_SURFACE`. La reconnaissance vérifie les coordonnées et les poids du quart de cercle, le rayon et la direction d’extrusion. Lors du premier essai, le STEP décrivait les mêmes surfaces par des B-splines rationnelles. Sur la variante Ø10/Ø5, FreeCAD/OCCT classait huit points de matière dans le vide sur un plan tangent au petit trou. Les coordonnées et les empreintes des fichiers du premier essai sont conservées dans [le contre-exemple](../../tests/counterexamples/fieldcraft_compact_x80_y40_cylindrical_step_supports.json).

FreeCAD/OCCT a relu cinq solides valides à une coquille, avec 38 faces chacun. Les 7 680 points des deux grilles, les 80 sondes radiales et les 40 sondes sur les plans tangents concordent avec les références. Sur les cinq paires de STEP, les plus grands écarts absolus mesurés par intégration adaptative sont 2,19 × 10⁻¹¹ mm³ en volume et 1,64 × 10⁻⁷ mm² en aire.

## Calcul sur la partition

Sur l’équerre d’origine, le calcul Q2 sans raffinement reproduit la résultante de 2 000 N sur 320 mm² avec 171 degrés de liberté scalaires, contre 705 auparavant. À quadrature d’ordre 16, l’énergie de déformation vaut 135,445 N·mm avec dix blocs et 144,124 N·mm avec l’ancienne partition. L’écart de 6,02 % porte sur deux espaces d’analyse différents. La validation géométrique laisse ouverte la convergence en contrainte et en énergie aux raffinements supérieurs.

## Fichiers et rejeu

[Télécharger les cinq équerres et le prisme témoin](fieldcraft-replay/fieldcraft-compact.zip) · [Ouvrir l’équerre dans la galerie](https://jeffwitz.github.io/SolidCraftTrivariate/?model=fieldcraft_x64_y64)

Depuis un clone du dépôt, les commandes suivantes reconstruisent le prisme témoin.

```sh
mkdir -p build
cargo run -- build tests/references/four_sector_plate.recipe.json build/plate.model.json
cargo run -- audit build/plate.model.json
cargo run -- export-step build/plate.model.json build/plate.step
```

Pour l’équerre, utiliser la recette `docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.recipe.json`. Le constructeur choisit la partition depuis les esquisses. Les modèles à 48 blocs restent auditables avec le plan enregistré dans leurs fichiers.
