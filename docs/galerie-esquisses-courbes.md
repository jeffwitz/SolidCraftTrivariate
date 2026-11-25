# Esquisses courbes de la galerie

Les six recettes reprennent la plaque à poche borgne, les deux cavités fermées, les deux sections en C et le prisme en L de la série polygonale. Le constructeur Rust conserve les identifiants des blocs et des raccords du plan d’esquisse. Il remplace les cartes trilinéaires par des volumes de Bézier de degrés (2, 2, 1), à contrôles homogènes entiers de poids commun. Ces cartes polynomiales sont des cas particuliers de NURBS. Les arcs et les surfaces figurent dans les coordonnées du modèle et dans l’export STEP.

## Règles de construction

Le rectangle de chaque opération donne le centre et les dimensions de son profil courbe. Pour un trou traversant, le constructeur rapproche les quatre coins du centre par un facteur 2/3. Il relie les coins successifs par quatre arcs quadratiques. Le contrôle central de chaque arc est décalé vers l’extérieur de la moitié de la dimension transverse du rectangle réduit. Les tangentes paramétriques coïncident aux jonctions. L’enveloppe en X et Y du trou reste celle du rectangle source.

Pour une poche borgne ou une cavité fermée, les facteurs correspondants sont 4/5 et 1/4. Les quatre arcs ont des tangentes distinctes aux jonctions. Cette construction conserve un bloc de fond régulier jusque dans ses coins. Sur chaque pièce, le côté extérieur de plus petite ordonnée est également courbe. Son contrôle central est déplacé vers l’extérieur d’un huitième de la hauteur du profil dans le plan XY. La flèche géométrique au milieu de l’arc vaut la moitié de ce déplacement.

Les contrôles des bords sont partagés entre les cartes voisines. Le centre du réseau biquadratique est obtenu par interpolation de Coons. Les plans sont ensuite extrudés suivant la direction prescrite par la recette. Pour le C, la face supérieure est translatée de (3, −2, 7) mm.

| Configuration | Blocs | Interfaces | Volume exact en mm³ |
|---|---:|---:|---:|
| Poche borgne | 9 | 16 | 48 976 |
| Cavité centrée | 14 | 28 | 247 952 / 25 |
| Cavité décalée | 14 | 28 | 247 952 / 25 |
| C à deux trous traversants | 9 | 10 | 45 878 / 27 |
| L avec trou traversant et cavité | 26 | 51 | 49 792 / 9 |
| C avec trou déplacé de 1 mm | 9 | 10 | 45 878 / 27 |

La borne de comptage des blocs s’applique aux partitions communes aux couches, avec raccords face entière et quatre jonctions prescrites par profil. Ces quatre jonctions peuvent être lisses sur un trou traversant. Leur conservation comme sommets du plan fait partie de l’hypothèse de comptage. [Démonstration combinatoire](partition-poche-comparaison.md)

## Rejouer une pièce

Depuis la racine du dépôt, cet exemple reconstruit le C à deux trous, lance l’audit natif et l’oracle indépendant, puis écrit le STEP.

```sh
cargo build --release --locked --offline
mkdir -p build/curved-gallery
cli=target/release/solidcraft-trivariate
recipe=docs/site/data/c_section_two_holes/c_section_two_holes.recipe.json
model=build/curved-gallery/c.model.json
"$cli" build "$recipe" "$model"
"$cli" audit "$model"
python3 tests/references/curved_sketch_gallery.py --recipe "$recipe" --model "$model" \
  --expected-block-count 9 --output build/curved-gallery/c.audit.json
"$cli" export-step "$model" build/curved-gallery/c.step
```

L’audit Rust rejoue la recette, compare les traces orientées et certifie la positivité du jacobien sur chaque cube paramétrique fermé. L’oracle Python vérifie les contrôles des arcs sur les parois, borne le déterminant en arithmétique rationnelle et intègre le volume de deux manières indépendantes. La première intègre le déterminant des cartes. La seconde intègre les contours d’esquisse par la formule de Green et somme les volumes des couches. La galerie permet de télécharger chaque recette, le modèle, l’audit et le STEP. [Données de la galerie](site/data)

Les six exports STEP ont été relus avec FreeCAD et OCCT 7.8.1. Le lecteur retrouve un solide valide par fichier, avec deux coquilles fermées pour les pièces à cavité interne et une pour les trois autres. Les écarts relatifs de volume avec l’intégrale exacte sont inférieurs à 10⁻¹⁴. L’option `--step build/curved-gallery/c.step` de l’oracle ajoute cette vérification lorsque les modules Python FreeCAD et Part sont disponibles.
