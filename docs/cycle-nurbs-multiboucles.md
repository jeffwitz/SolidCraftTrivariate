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
