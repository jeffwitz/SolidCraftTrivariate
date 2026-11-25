# Cycle Luna — esquisses NURBS multiboucles

Départ de `81cb61d`, branche `campagne-nurbs-multiboucles`. Ce cycle succède à la fenêtre de dix heures clôturée. Il porte sur les esquisses à plusieurs contours, leur partition automatique, les poches et l’espace IGA.

## M1 — contours en bandes, trou traversant et poche

La famille `nurbs_multiloop_v1` reçoit `bounds_xy_mm:[[xmin,xmax],[ymin,ymax]]`, `height_mm` et `loops`. Chaque entrée comporte `id`, `lower`, `upper` (NURBS du format nurbs_sketch_strip_v1) et `z_interval_mm:[lo,hi]`. Les deux extrémités de chaque paire sont reliées par des segments verticaux. Le contour extérieur est un rectangle. Les projections en x des contours intérieurs doivent être strictement séparées, les contours strictement à l’intérieur du rectangle. Ce premier domaine inclut des contours à flancs non coniques. Le planificateur trie les événements physiques et prolonge les séparatrices à travers les intervalles libres. Les cellules suivent les couches déduites de tous les débuts et fins d’enlèvement. Les noms des pièces n’interviennent pas.

Construction Rust dans src/nurbs_multiloop.rs, intégration et audit de persistance depuis le CLI. Oracle géométrique indépendant par fonctions B-spline sur les nœuds sources, intégrale des aires et soustraction des volumes par intervalle axial. Comparer le STEP relu à 1e−8 relatif, les positions de paroi à 1e−10 mm et les classements intérieur/extérieur/face. Prospects : même boucle en trou traversant puis en poche borgne, deuxième boucle indépendante de forme différente et avec un autre niveau, permutation des entrées et renommage. Vérifier l’invariance de matière, la provenance du fond de poche et des parois. Refuser un contour touchant le stock, des projections qui se chevauchent et un intervalle dégénéré. Vérifier la conservation de la cible publiée sur refus. Distinguer les décisions prévues, les réussites et tout refus d’un cas obligatoire.

L’oracle IGA vérifie l’identité des degrés de liberté des faces, le champ affine et son gradient, puis l’équilibre d’une contrainte constante au moyen des efforts de volume et de bord. Seuil du résidu relatif 1e−8 à l’ordre 16. Une estimation numérique de volume ou de conditionnement se rapporte à ses points et à son ordre. Rejouer make check et Clippy strict avant intégration. Le journal consigne les révisions et résultats mesurés.

## Flotte activée

Les métadonnées de session et le modèle effectif du premier tour identifient trois exécutions `gpt-6-luna` distinctes.

| Rôle | Session |
|---|---|
| Construction Rust | `01a21208-822a-7f02-8a5d-a90895f7a829` |
| Oracle géométrique et STEP | `01a21208-b7a7-79d3-93c9-bdc4262b1012` |
| Oracle IGA | `01a21208-e10c-78d1-a96e-9423ec00f071` |

## Premier essai prospectif

Le profil initial de la deuxième boucle présente trois pôles successifs à x = 70 mm. Au nœud intérieur t = 2, sa dérivée en x est nulle. Cette entrée échoue donc au critère de régularité du lot. La recette originale est conservée dans `tests/counterexamples/nurbs_multiloop_stationary.json`. Le prospect de remplacement utilise les abscisses 55, 62, 67, 70, 74 et 80 mm. L’intégrale indépendante donne 112940,85164966024 mm³ de matière pour la recette mixte de remplacement. Les seuils restent ceux du préenregistrement.

Le constructeur est intégré à `8025838`. Sur un instantané de l’index, `make check` a exécuté 395 tests Rust et 58 tests Python ; Clippy strict passe. Le modèle mixte compte 45 blocs, 75 interfaces et 120 côtés extérieurs. Trois altérations isolées (contrôle homogène, nom de boucle, étiquette de face) donnent toutes `multiloop_stored_partition_mismatch`.

### Qualification de M1

Le trou traversant, la poche, leur combinaison et les contours permutés/renommés produisent respectivement 12, 27, 45 et 45 blocs. Chaque STEP relu avec OCC 7.8.1 comporte un solide fermé valide. Le volume du cas mixte atteint 112940,85164966017 mm³, pour 112940,85164966024 mm³ selon l’intégrale des courbes de l’esquisse. Les 1515 classifications concordent ; l’écart maximal mesuré sur les neuf points contrôlés par face est 2,85e−14 mm. Les quatre entrées invalides ont conservé la cible octet pour octet. La règle de Gauss de la première tentative de contrôle volumique évaluait les blocs hors de leur cube unité. Les deux rapports permettent de rejouer le diagnostic et sa correction. [Qualification STEP](results/nurbs-multiloop-delivery-20261125.json) · [Premier diagnostic](results/nurbs-multiloop-delivery-initial-20261125.json) · [Référence des profils](results/nurbs-multiloop-geometry.json)

Le cas mixte possède 1800 degrés de liberté scalaires Q3, 1200 correspondances nodales d’interface et 2880 contrôles du champ affine et de son gradient. Le résidu relatif d’équilibre à q16 vaut 2,65e−12. Le test historique appliquait les poids deux fois aux coordonnées homogènes. Après correction, son résidu vaut 1,96e−12. Les deux oracles M1 ont utilisé le binaire de SHA-256 `72d0b6163bcd2eb5468b8783b9abaf52fc562cd87c5fffe7cc821889c9304040`. [Qualification IGA](results/nurbs-multiloop-iga.json) · [Rejeu historique](results/nurbs-patch-equilibrium-20261125.json) · [Mesure antérieure](results/nurbs-patch-equilibrium-20261125-double-weight-diagnostic.json)

## M2 — assemblage des contours à partir d’arêtes indépendantes

Départ du constructeur `8025838`. La recette `nurbs_edge_loops_v1` reçoit le stock et la hauteur de M1, un tableau non ordonné `edges:[{id,curve}]` et des opérations `cuts:[{id,seed_edge,z_interval_mm}]`. Chaque opération désigne une arête du contour. Le constructeur doit découvrir les autres arêtes, les deux flancs et leur sens depuis leurs extrémités homogènes exactes. Domaine initial : composantes simples de quatre arêtes, deux fermetures verticales linéaires et deux flancs NURBS synchronisables au sens M1. Les composantes sont disjointes en projection x. Les noms et l’ordre des arêtes ou des opérations ne déterminent pas le partitionnement.

Fichiers prévus : `src/nurbs_edge_loops.rs`, intégration de la famille dans les entrées Rust existantes, tests ciblés et recettes/oracles dans `tests/references`. Prospect : deux boucles libres, une poche ouverte par le haut et un trou traversant, arêtes permutées, orientations inversées, identifiants et jauges changés. Retrouver le volume de l’intégrale indépendante (écart relatif ≤ 1e−8), les positions de paroi à 1e−10 mm, les faces du fond et les traces conformes. Vérifier le champ affine et le résidu d’équilibre Q3 à q16 ≤ 1e−8. Refuser une chaîne ouverte, une bifurcation, une fermeture oblique et deux opérations sur la même composante. Vérifier la conservation de la cible sur chaque refus. Les arrangements d’arêtes qui se coupent en leur intérieur et les boucles à nombre quelconque d’arêtes dépassent ce lot.

M2 est suspendu à la demande de Jeff après le passage des trois tests Rust ciblés et la construction du prospect à 45 blocs. La qualification indépendante de M2 reste ouverte.

## S1 — plaque à trou rectangulaire vif et trou NURBS fermé G2

Révision de départ `58c4409`. Priorité donnée par Jeff au contrôle des raccords. Fichiers dédiés prévus dans `src/g2_plate.rs`, `tests/g2_plate.rs`, `tests/references/g2_plate*` et `docs/results/g2-plate*`. Les trois exécutions Luna du tableau précédent prennent la construction Rust, la comparaison B-Rep et l’oracle des jets de courbe et du champ IGA.

Le prospect est une plaque extrudée avec deux découpes traversantes. Les quatre coins vifs du rectangle doivent être détectés. L’autre profil est une B-spline rationnelle cubique périodique à nœuds simples et poids positifs, non conique, sans segment de fermeture ajouté. La recette et les positions des pôles doivent être communes à toutes les représentations. Vérifier à tous les raccords, fermeture incluse, la position (1e−10 mm), la tangente unitaire orientée (1e−10) et le vecteur de courbure géométrique (1e−10 mm⁻¹). Un témoin de raccord G1 non G2 doit franchir les deux premiers tests et échouer au troisième. Mesurer aussi les raccords entre faces extrudées : normale et courbure selon la section, courbure le long des génératrices. Les arêtes du plat sont des jonctions vives ; l’espace d’analyse assemblé est C0.

Essayer en premier le backend B-Rep réel de SolidCore sur cette esquisse et son extrusion. Consigner sa révision, la commande, la voie d’accès à la géométrie et le résultat. Un essai OCCT est admis après un échec constaté et tracé de SolidCore. Comparer le volume à l’intégrale indépendante des courbes (1e−8 relatif), la paroi aux courbes extrudées et les classifications dans les deux vides et la matière. Vérifier traces, reproduction affine et équilibre Q3 à q16 (1e−8 relatif). Si le certificat de régularité, la fermeture ou la comparaison échoue, conserver le cas avec le détail de l’échec.
