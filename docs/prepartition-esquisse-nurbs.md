# Préparer la partition dans l’esquisse NURBS

Une opération de la plaque croisée spécifie un repère d’esquisse, une direction, une longueur d’extrusion et les boucles à soustraire. Le même format de courbe décrit le contour rectangulaire et le contour libre. Chaque courbe fournit son degré, ses nœuds et ses pôles homogènes `[wx, wy, w]`. [Recette des deux esquisses](../tests/references/g2_cross_sketches.recipe.json)

## Décider où placer les lignes de construction

1. Extraire les travées sur l’intervalle actif de la courbe. La fonction Rust `extract_active_bezier` évalue la forme polaire de la B-spline dans les coordonnées homogènes, en arithmétique rationnelle vérifiée. Elle accepte les degrés 1 à 3, les nœuds entiers à une échelle commune, les poids positifs et les séquences périodiques étendues comme les courbes serrées aux extrémités.
2. Examiner les deux côtés de chaque nœud et la fermeture. À un nœud intérieur de multiplicité m, une base de degré p garantit ordinairement C^(p−m) pour m ≤ p. Les pôles et les poids peuvent donner davantage de régularité. `plan_loop` expose le nœud, sa multiplicité, la travée sortante et la continuité géométrique mesurée exactement sur les jets.
3. Conserver la provenance des jonctions. Une discontinuité de position ou une vitesse nulle demande un traitement explicite. Un changement de tangente marque un coin. Une tangente commune avec un saut de courbure marque une jonction G1. Un nœud introduit pour raffiner un côté rectiligne peut conserver la continuité G2.
4. Ajouter les événements imposés par **l’autre extrusion**. Les parois du rectangle transverse rencontrent le profil périodique à des paramètres qui peuvent se trouver à l’intérieur des travées. Chaque événement géométrique doit porter les identités des deux esquisses. Ces événements décident où subdiviser les travées et où prolonger les lignes de construction en faces. Un nœud de représentation et une intersection entre opérations ont des provenances distinctes.

## Critère géométrique implémenté

Pour les deux jets réguliers d’une courbe plane orientée, le test commence par l’égalité des positions. Il vérifie ensuite que les vitesses sont colinéaires et de même sens. Si la vitesse à droite est λ fois celle de gauche avec λ > 0, la jonction est G2 si la différence entre l’accélération à droite et λ² fois l’accélération à gauche est parallèle à la tangente. Ce test autorise un changement de vitesse ainsi qu’une accélération tangentielle différente. Il utilise des opérations rationnelles et signale les dépassements arithmétiques. [Implémentation Rust](../src/sketch_junction.rs)

Les essais portent sur un rectangle cubique à quatre coins et une subdivision supplémentaire, un cercle et une ellipse rationnels à nœuds doubles, une boucle cubique périodique non uniforme et deux travées soumises à des changements homographiques de paramètre. Le diagnostic géométrique porte sur les courbes de l’esquisse. Qualifier leur extrusion et leur intersection demande aussi des contrôles de régularité volumique, de raccords, de couverture et de fidélité au solide demandé. [Tests exécutables](../tests/sketch_junction.rs)
