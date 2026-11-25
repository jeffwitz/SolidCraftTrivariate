# Construire une extrusion depuis deux courbes NURBS

La famille `nurbs_sketch_strip_v1` prend deux courbes de conception, `lower` et `upper`, et une hauteur d’extrusion en millimètres. Chaque courbe contient son degré (1 à 3), son vecteur de nœuds entier clampé et ses contrôles homogènes entiers `[wx, wy, w]`. Les poids doivent être positifs. Les nœuds peuvent représenter des paramètres rationnels après application d’un même facteur d’échelle. La recette d’exemple emploie deux frontières cubiques non rectilignes, deux travées et des poids variables. [Recette](../tests/references/nurbs_strip.recipe.json)

L’insertion de nœuds se fait en arithmétique rationnelle vérifiée. Chaque travée est ramenée à une jauge entière commune. Le constructeur vérifie que les deux frontières ont les mêmes intervalles, le même degré, la même abscisse rationnelle et des poids proportionnels par travée. Il exige une abscisse strictement monotone dans la même direction sur toute la chaîne et un écart strictement positif entre les ordonnées. Les bornes de Bernstein du numérateur de la dérivée et de l’écart établissent ces propriétés sur les intervalles fermés. La monotonie permet de retrouver la travée et le paramètre longitudinal depuis x. La position entre les deux ordonnées fixe le paramètre transverse, et la hauteur fixe le dernier paramètre. Les cartes partagent leurs faces de bout. Le modèle indique pour chaque carte l’intervalle du paramètre source dans le sens local u croissant.

Depuis la racine du dépôt, les commandes suivantes utilisent le même modèle audité.

```sh
mkdir -p build/nurbs-demo
cargo run --locked --offline -- build tests/references/nurbs_strip.recipe.json build/nurbs-demo/strip.json
cargo run --locked --offline -- audit build/nurbs-demo/strip.json
cargo run --locked --offline -- export-step build/nurbs-demo/strip.json build/nurbs-demo/strip.step
cargo run --locked --offline -- export-analysis build/nurbs-demo/strip.json 16 build/nurbs-demo/space.json
```

Le paquet `solidcraft_nurbs_h1_space_v1` élève exactement les réseaux au degré maximal rencontré et identifie les degrés de liberté scalaires des faces selon leur orientation. Il fournit une base rationnelle conforme C0 entre les cartes, les contrôles homogènes, la numérotation et un volume intégré à l’ordre de Gauss demandé. Le bilan mécanique Q2 sur les quadriques et les mesures géométriques des nouvelles cubiques ont des protocoles distincts. Les profils fermés arbitraires, les courbes non monotones et la reparamétrisation générale de deux NURBS restent en dehors du certificat de cette bande. [Croisement des profils](nurbs-cross.md)

![Section, polygones de contrôle et limite des deux travées](figures/nurbs-esquisse-20261125.png)

Pour la recette fournie, le volume de référence vaut 45,96534510602271 mm³. À l’ordre 8, l’intégration du paquet H¹ donne 45,9653838425 mm³. La relecture de masse OCC utilise une intégration adaptative à tolérance 1e−12. Les trois cas admis et trois refus, les mesures par défaut et adaptatives ainsi que les réglages de quadrature sont enregistrés dans les rapports. [NURBS et STEP](results/nurbs-sketch-20261125.json) · [Équilibre du patch](results/nurbs-patch-equilibrium-20261125.json)
