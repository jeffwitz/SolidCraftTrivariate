# Perçages transverses définis par des courbes NURBS d’esquisse

## Contrat du lot

Partir de la branche de campagne après `ce372ff`. Conserver les tests de Steinmetz. L'entrée comporte deux esquisses dans les plans XY et ZY. Chaque profil fermé comprend deux flancs NURBS et deux segments de fermeture à ordonnée constante. L’enlèvement de la première esquisse traverse suivant Z, celui de la seconde suivant X. Les coordonnées de contrôle homogènes et les nœuds sont donnés par la recette ; les données décimales rationnelles peuvent être encodées sur une échelle entière commune. Les bornes du brut sont exprimées en mm.

Les flancs sont des NURBS clamped de degré 1 à 3, à poids positifs, orientées dans le sens des y croissants. Dans une bande où plusieurs flancs sont actifs, leurs bases sont mises sur un découpage et un degré communs. Leur paramétrisation homogène de y doit alors être identique à un facteur d’échelle près. Les flancs d’un même perçage ont le même intervalle de paramètres. Les deux perçages peuvent commencer et se terminer à des niveaux différents. Leurs plages doivent employer des paramètres globalement compatibles avec l’ordre physique des événements.

## Critères fixés avant essai

- Tester des courbes cubiques rationnelles à plusieurs travées, avec nœud intérieur simple, poids variables et flancs différents. Comparer à une évaluation indépendante sur les NURBS d’entrée. Tester le refus d’une désynchronisation ou d’un croisement de flancs et la conservation du fichier cible.
- Vérifier sur tout le domaine le signe des poids, la croissance stricte de y et l’ordre des quatre frontières dans chaque section. Le volume de référence est l’intégrale des longueurs matérielles en x et z multipliées entre elles. Comparer par quadrature indépendante à 10⁻⁹ en relatif. Les tests de positions et d’incidences complètent ce contrôle.
- Conserver les degrés et les coefficients homogènes des courbes de conception. Refuser un dépassement de capacité de l’arithmétique exacte. Séparer les capacités réellement vérifiées du CLI, de l’export STEP et de l’analyse.

## Construction du croisement à pôles lisses

La recette `nurbs_cross_oval_v1` fournit le brut, le centre et deux quarts de contour NURBS, l’un dans XY et l’autre dans ZY. Les réflexions des courbes ferment les deux sections. Le constructeur contrôle les tangentes aux axes de réflexion, les poids, l’ordre des contrôles et la garde au brut. Après insertion des nœuds et élévation exacte du degré, les deux méridiens partagent l’ordonnée et le poids à chaque paramètre. Si leurs points sont (a(t), e(t)) et (b(t), e(t)), l’arête commune en 3D est exactement (a(t), e(t), b(t)) en coordonnées locales. Les largeurs a et b peuvent varier indépendamment.

Le constructeur place des lignes dans les deux esquisses au même événement de nœud. Pour une courbe à une seule travée, il la subdivise à mi-paramètre. Dans chaque octant il prolonge les deux méridiens jusqu’à une garde rectangulaire intérieure, puis prolonge cette surface commune suivant X ou Z jusqu’au brut. Le bloc de coin reprend le paramétrage de la garde. Les deux directions se rejoignent sur une face entière. Chaque octant comporte deux blocs par travée et un bloc de coin. Les changements de topologie aux pôles sont traités par les jonctions entre octants.

## Mesures et reproduction

Le croisement à bouts plats donne 26 blocs et 36 interfaces. Le volume obtenu par intégration des blocs est de 50 835,02450412389 mm³. Une intégration indépendante sur les courbes sources évaluées par de Boor donne 50 835,02450413309 mm³. Sur chacune des quatre premières travées, sept points rationnels donnent un rang exact de six pour les coefficients d’une conique. Aucune conique ne contient ces travées. L’oracle vérifie également 702 points pris dans la matière.

Pour les deux profils à pôles lisses, la construction donne 40 blocs et 68 interfaces. Le paquet IGA de degré trois contient 1 588 degrés de liberté scalaires avec raccords C0. À l’ordre de Gauss 24, son volume vaut 63 390,98047047747 mm³ contre 63 390,98047048173 mm³ par inversion des courbes sources et intégration des sections.

La variante explicite `nurbs_cross_bounded_v1` accepte des ordonnées paramétrées indépendamment sur des flancs monotones à bouts plats. Elle certifie chaque corde par enveloppe des contrôles, superpose les niveaux physiques et vérifie les gardes. La tolérance demandée borne le déplacement latéral de chaque flanc. Le déplacement 3D de la frontière a pour majorant conservatif deux fois cette tolérance. Avec une tolérance de 0,25 mm, le témoin donne 82 blocs et un écart latéral maximal mesuré de 0,147 mm.

```sh
mkdir -p build/nurbs-demo
cargo run --locked -- build tests/references/nurbs_cross_oval.recipe.json build/nurbs-demo/oval.json
cargo run --locked -- audit build/nurbs-demo/oval.json
cargo run --locked -- export-step build/nurbs-demo/oval.json build/nurbs-demo/oval.step
cargo run --locked -- export-analysis build/nurbs-demo/oval.json 16 build/nurbs-demo/oval-space.json
```

Les conditions de monotonie et de synchronisation sont des conditions suffisantes de ce constructeur. Une esquisse à plusieurs boucles libres, des axes obliques ou une poche borgne qui intercepte le croisement demandent une autre décomposition. En général, deux surfaces rationnelles peuvent se couper suivant une courbe de genre positif. Une représentation exacte de leur arête par une chaîne finie de NURBS régulières ne peut alors pas être exigée. [Katz et Sederberg, *Genus of the intersection curve of two rational surface patches*](https://www.sciencedirect.com/science/article/pii/0167839688900064)
