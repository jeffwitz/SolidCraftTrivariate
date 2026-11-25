# Congés et chanfreins sur la partition de l’anneau vif

Le prototype reprend les quatre volumes de `full_annulus_v1`, leurs identifiants, leurs quatre faces de raccord et les coefficients angulaires dyadiques de la pièce vive. Les variantes affichées utilisent la même base triquadratique, le même vecteur de nœuds par direction et les mêmes poids. Chaque volume a 63 pôles (3 × 3 × 7). Les nœuds supplémentaires sont internes aux volumes. Dans chacun des trois états on conserve les quatre secteurs, avec trois travées d’intégration axiales par secteur.

| État | Volume Gauss 16 (mm³) | Régularité |
| --- | ---: | --- |
| [Vif](site/data/models/annulus-four-sharp.json) | 301,592894744620 | Déterminant positif sur les travées fermées |
| [Quatre chanfreins](site/data/models/annulus-four-chamfer.json) | 299,236700254428 | Déterminant positif sur les travées fermées |
| [Quatre congés](site/data/models/annulus-four-fillet-singular.json) | 300,581607414644 | Positif dans la matière, nul sur quatre cercles de contact. IGA non qualifiée |

## Carte du congé

Dans la travée du bas, s est le paramètre radial (0 à l’extérieur, 1 à l’intérieur), t le paramètre axial local et R=1/4 mm. Poser a=2s−1, D=1+t², d=R(1−t)²/D, b=1−d, E=t(1−t)/D et k=R/2. Le rayon nominal et l’altitude sont

    ρ = 3 − a b
    z = 2R t²/D + k(1−a²)E.

Sur les deux parois a=±1, le terme en k s’annule. Les fonctions (1−t²)/D et 2t/D tracent un quart de cercle exact dans le profil méridien nominal. La révolution reprend l’approximation dyadique angulaire d’origine. La travée du haut est la symétrique de celle du bas autour de z=4 mm. Dans la travée centrale, le profil est rectangulaire.

Le terme en k régularise **l’intérieur des faces d’extrémité**. Le facteur méridien du jacobien orienté, à une constante positive près, vaut b z_t + 2 a² k b′ E. Ici b>0, b′≥0, E≥0 et

    D² z_t = 4R t + k(1−a²)(1−2t−t²).

Pour t>0 ce dernier terme est strictement positif, car 1−2t−t² ≥ −2t et 4R>2k. Pour t=0, il est positif dès que |a|<1. Les seuls zéros sont donc les deux cercles à z=0 et leurs deux symétriques à z=8, aux rayons nominaux 2,25 et 3,75 mm. Le facteur angulaire est strictement positif comme sur la pièce initiale. Pour vérifier la couverture, à rayon fixé on résout a=(3−ρ)/b(t). L’ensemble des t permis est un intervalle car b est croissante, et z y est strictement croissant dans l’intérieur par la formule du jacobien. Le bord décrit le profil arrondi, et les quatre secteurs angulaires couvrent chacun leur quadrant une fois.

Une arête de cube entre deux faces physiques tangentes impose cette perte de rang sur ces quatre secteurs. Un test du champ IGA et une qualification pour cette classe de cartes singulières restent à mener. Le rapport antérieur [sur la condition de jacobien strict](annulus-four-block-topology.md) explique pourquoi la variante régulière utilise une autre partition.

## Contrôles reproductibles

Le contrôle indépendant opère directement sur les modèles de la galerie. Sur chaque travée et pour chaque rangée angulaire, neuf comparaisons en fractions exactes déterminent l’identité des polynômes biquadratiques homogènes (numérateurs et dénominateur, composante par composante) avec la formule analytique. Les quatre faces de raccord sont identiques coefficient par coefficient. Les quadratures de Gauss 12 et 16 concordent avec le volume analytique vif ou le STEP traité de référence à moins de 10⁻¹⁰ mm³.

```sh
cargo test --offline --example annulus_original_edges
cargo run --offline --example annulus_original_edges -- fillet > build/annulus-original-fillet.json
python3 tests/references/annulus_original_edge_states.py
```

[Constructeur Rust](../examples/annulus_original_edges.rs) · [Oracle indépendant](../tests/references/annulus_original_edge_states.py)
