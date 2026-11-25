# Anneau à connectivité fixe

L’essai part de la révision `d25b6b6`. L’anneau a pour rayons 2 et 4 mm et pour hauteur 8 mm. L’opération traite les quatre arêtes circulaires ensemble, avec un rayon ou un retrait de 0,25 mm. Deux expériences répondent à deux décisions que le constructeur peut prendre dès l’esquisse.

| Plan de construction | Vif | Chanfrein | Congé | Cellules de Bézier pour l’audit par état |
| --- | --- | --- | --- | --- |
| Éditer les quatre secteurs de la galerie | 4 volumes | 4 volumes | Obstruction au jacobien strict avec ces incidences | 12 |
| Préparer une section à cinq cartes et deux demi-tours | 10 volumes | 10 volumes | 10 volumes | 128 |

## Édition des quatre secteurs d’origine

Le prototype lit les quatre blocs de `full_annulus_v1` et conserve leurs identifiants, leurs orientations, leurs pôles angulaires et leurs poids à une échelle homogène globale près. Il prépare deux nœuds axiaux internes à 0,25 et 7,75 mm dans chacun des quatre blocs de degrés (2,1,1). Les deux états ont 24 pôles par volume. À hauteur z, le rayon varie entre 2+d et 4−d. Pour le chanfrein, d=max(0, 0,25−min(z,8−z)). L’état vif utilise d=0. L’édition déplace deux rangées radiales à chaque extrémité. Le volume du profil idéal chanfreiné vaut 95,25π mm³.

Les quatre raccords entre secteurs coïncident coefficient par coefficient. Le certificat de Bernstein appliqué aux douze cellules fermées donne une borne de déterminant supérieure à 0,82079 mm³ dans l’état vif et à 0,69254 mm³ dans l’état chanfreiné, en paramètres locaux de cellule. Le rayon, la hauteur et l’angle parcourent chacun leur intervalle une seule fois, avec la même approximation dyadique du cercle que le constructeur d’origine. [Constructeur et tests](../examples/annulus_four_block_chamfer.rs) · [Vif](results/annulus-fixed-topology/four-sharp.json) · [Chanfrein](results/annulus-fixed-topology/four-chamfer.json)

## Obstruction pour le congé sur ces quatre volumes

Chaque secteur d’origine est le produit d’un quadrant angulaire par une section carrée. Ses deux faces d’extrémité angulaire sont des raccords internes. Deux faces physiques extérieures sont incidentes à chacune des quatre arêtes longitudinales issues des coins de la section. Après traitement des quatre arêtes, le bord méridien entier est lisse au sens du plan tangent, y compris aux jonctions des arcs et des droites.

Dans une cellule tensorielle qui touche une telle arête, prenons les deux faces u=0 et v=0 et un point intérieur de l’arête hors des nœuds angulaires. Les dérivées tangentielles des deux faces sont respectivement (Fv,Fw) et (Fu,Fw). Si le bord physique a un seul plan tangent, les trois colonnes Fu, Fv, Fw appartiennent à ce plan. Le déterminant s’annule donc sur le bord. Le raisonnement vaut pour les dérivées unilatérales de la cellule adjacente, quel que soit le nombre de nœuds internes. Il suffit d’une seule telle arête pour faire échouer le critère du projet. Le chanfrein laisse des cassures de tangente auxquelles ces arêtes paramétriques peuvent aboutir. [Témoin Rust du défaut de rang et des changements de paramètres](../tests/annulus_four_sector_obstruction.rs)

## Gabarit préparé depuis la section et son constructeur de révolution

Le constructeur joint les milieux des quatre côtés du profil à un losange intérieur centré sur (3,4), de demi-diagonales 0,5 et 2 mm. Les quatre bandes entre ce cœur et le contour portent chacune un coin physique au milieu de leur **bord de carte**. Le cœur et les bandes existent dans l’état vif. Chaque bande utilise quatre travées tangentielles. Deux d’entre elles couvrent la zone modifiable. Dans cette zone, les pôles extérieurs représentent suivant l’état deux branches du coin vif, deux segments colinéaires ou deux restrictions d’un arc circulaire rationnel tangent. Le recul et le rayon valent 0,25 mm sur les quatre coins.

Le balayage utilise deux demi-tours injectifs, avec deux travées rationnelles par demi-tour et une vitesse angulaire raccordée. Chaque volume retrouve le même voisin sur deux faces angulaires distinctes. Les dix identifiants, incidences, degrés, nœuds, multiplicités et poids sont identiques dans les trois états. Seules les coordonnées des pôles de la paroi changent. Le format JSON du prototype enregistre les dix réseaux NURBS complets. [Constructeur et audit Rust](../examples/annulus_prepared_volumes.rs) · [Vif](results/annulus-fixed-topology/ten-sharp.json) · [Chanfrein](results/annulus-fixed-topology/ten-chamfer.json) · [Congé](results/annulus-fixed-topology/ten-fillet.json)

Pour chaque état, les 128 cellules extraites ont toutes un déterminant strictement supérieur à 0,0048828 mm³ en coordonnées locales. Le contrôle des coefficients homogènes retrouve 352 interfaces de cellules et les 64 faces de cellules attendues sur la peau. Dans la section, les quatre bandes couvrent des quadrants méridiens disjoints. Dans chaque quadrant, le bord de bande est simple (côté du cœur, deux rails, contour extérieur). Le déterminant positif dans chaque cellule et les raccords orientés, avec ce bord de degré topologique un, donnent une seule couverture. La révolution préserve cette couverture puisque r reste strictement positif et que les deux demi-tours couvrent le cercle une fois.

L’oracle lit les réseaux JSON, extrait les éléments depuis les nœuds et intègre le déterminant en quadrature tensorielle de Gauss aux ordres 12 et 16. L’état vif est comparé à 96π. Les deux états traités sont comparés aux volumes préalablement mesurés après relecture des STEP indépendants. [Oracle de volume](../tests/references/annulus_prepared_volumes.py)

| État | Volume q16 (mm³) | Référence (mm³) |
| --- | ---: | ---: |
| Vif | 301,592894744620 | 301,592894744620 |
| Chanfrein | 299,236700254428 | 299,236700254422 |
| Congé | 300,581607414644 | 300,581607414643 |

L’audit en fractions exactes retrouve les quatre profils et leurs raccords sur les données enregistrées. Pour chaque congé, l’identité polynomiale du cercle de rayon 1/4 mm vaut sur les deux travées complètes. Les tangentes aux deux contacts et à la jonction entre les travées ont la direction attendue. Sur les 1 530 pôles stockés, le passage du vif au congé modifie 80 positions, et le passage du vif au chanfrein en modifie 120. [Contrôle du profil](../tests/references/annulus_prepared_profile.py)

## Dans quelle classe dix volumes sont-ils minimaux ?

Considérons une quadrangulation conforme de la section en disque, avec bord lisse, sans double arête ni identification, et des coins de carte réguliers. Soient I le nombre de sommets internes, B le nombre de sommets de bord et q le nombre de quadrilatères incidents au sommet. Euler et le comptage des côtés donnent F=I+B/2−1 et 4=Σᵢ(4−q)+Σᵦ(2−q). L’image d’un coin régulier a un angle strictement inférieur à π. Il faut donc au moins trois cartes en un sommet interne, et deux en un sommet du bord lisse. Il en résulte I≥4. Le contour simple a un nombre pair B≥4, donc **F≥5**. Le losange et ses quatre bandes atteignent cette borne.

Dans une construction produit de cette section par des intervalles angulaires, une seule carte angulaire sur un intervalle fermé identifierait ses deux extrémités. L’injectivité sur le fermé impose au moins deux intervalles. Les deux demi-tours atteignent cette borne et donnent **dix volumes dans cette classe**. Imposer les quatre quadrants comme volumes distincts conduit à vingt. La preuve ne porte ni sur les cartes périodiques avec auto-identification, ni sur toutes les décompositions héxaédriques 3D non produites par révolution.

## Rejeu et portée

```sh
cargo test --offline --test annulus_four_sector_obstruction
cargo test --offline --example annulus_four_block_chamfer
cargo test --offline --release --example annulus_prepared_volumes -- --nocapture
cargo run --offline --example annulus_prepared_volumes -- fillet 2 > build/annulus-ten-fillet.json
python3 tests/references/annulus_prepared_profile.py
python3 tests/references/annulus_prepared_volumes.py docs/results/annulus-fixed-topology/ten-{sharp,chamfer,fillet}.json
```

La qualification porte sur la géométrie et ses raccords. L’espace IGA, les degrés de liberté et les performances mécaniques du gabarit à dix volumes restent à mesurer. Le parcours de production doit encore intégrer les cartes logiques multi-travées du prototype. Le planificateur devra compter à la fois les volumes persistants, les travées, les pôles et le coût de calcul. L’essai antérieur à 96 blocs mono-travée utilisait 96 cellules géométriques. Le gabarit à dix volumes en utilise 128 dans cette réalisation.
