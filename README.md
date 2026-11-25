# SolidCraftTrivariate

Noyau volumique natif en Rust pour construire et auditer des cartes trivariées rationnelles, puis les exporter en STEP et préparer leur espace d’analyse IGA. Le parcours géométrique commence par une esquisse ou une recette, détermine ses lignes de partition, puis construit et audite les volumes dans la même crate. Python fournit des oracles indépendants et des références historiques.

La galerie présente les six configurations à profils courbes. Chaque fiche montre la coupe du solide et les arcs de son esquisse, avec la recette, le modèle, l’audit et le STEP à télécharger. [Galerie interactive](https://jeffwitz.github.io/SolidCraftTrivariate/) · [Construction et vérifications](docs/galerie-esquisses-courbes.md)

Le planificateur polygonal compare les séparatrices locales à la partition historique sur la même recette. La poche borgne du corpus passe de 17 à 9 volumes et le C à deux trous de 35 à 9. Ces deux résultats atteignent la borne topologique dans la classe de partitions conformes étudiée. Le conditionnement des cartes entre aussi en jeu pour le calcul. [Algorithme](docs/architecture.md#planificateur) · [Corpus, preuves et mesures](docs/partition-poche-comparaison.md) · [Construction depuis les esquisses](docs/partition-esquisses.md)

## Démarrage rapide : esquisse, plan, volumes

Prérequis : Rust/Cargo, Python 3.10+, Git, Bash et Make. Depuis la racine, le parcours Rust initial ne dépend pas des clones upstream ni des oracles Python :

```sh
cargo build --release --locked --offline
make demo-native
```

Pour installer les six références de `references.lock` et lancer l’ensemble des références et tests, faire une fois `make bootstrap`, puis `make check`. `make check` vérifie la présence de ces clones; il ne les installe pas implicitement.

`plan-sketch` expose les jonctions du contour et les niveaux de coupe induits par les autres opérations avant la création des blocs. Ici, le profil cubique périodique est soustrait par extrusion à un rectangle transversal. Le plan est enregistré dans un JSON consultable :

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" plan-sketch tests/references/g2_cross_sketches.recipe.json > build/demo/g2-plan.json
python3 -m json.tool build/demo/g2-plan.json
```

Le rapport rejoué contient huit jonctions G2 sur le profil périodique, les coins G0 du rectangle, quatre croisements exacts aux paramètres 7/2, 13/2, 5/2 et 15/2, ainsi que les niveaux transversaux 4 et 8 mm. Ces valeurs décrivent cette recette. Les nœuds introduits pour représenter une courbe et les croisements entre deux esquisses gardent des provenances distinctes.

Construire la représentation volumique à partir de la même recette, puis auditer et exporter le solide :

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" build tests/references/g2_cross_sketches.recipe.json build/demo/g2-cross.json
"$cli" audit build/demo/g2-cross.json
"$cli" export-step build/demo/g2-cross.json build/demo/g2-cross.step
```

Sur cette recette, l’audit rejoué rapporte 32 blocs, 46 interfaces partagées, 100 faces extérieures, une borne inférieure du déterminant de 5,729862499996598 mm³ et deux perçages traversants. La borne certifie la régularité locale des cartes sur leur domaine; elle ne suffit pas à elle seule à établir la fidélité globale d’une autre recette.

### Importer une timeline SolidCraft

`import-solidcraft-sketch` lit une timeline canonique exportée de SolidCraft; il ne lance ni n’exécute le script utilisateur. Le sous-ensemble qualifié ici est une esquisse polygonale orthogonale, son extrusion, un congé et deux perçages traversants d’axes distincts. La fixture versionnée permet de rejouer l’import, l’audit et l’export STEP :

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" import-solidcraft-sketch tests/references/solidcraft_sketch_bracket.timeline.json build/demo/bracket.json
"$cli" audit build/demo/bracket.json
"$cli" export-step build/demo/bracket.json build/demo/bracket.step
```

Le cas de référence rejoué donne 48 blocs, 74 interfaces partagées, 140 faces extérieures, deux perçages et une borne certifiée du déterminant de 0,00854492187498991 mm³. Les détails de l’entrée canonique et les refus atomiques sont dans [le test d’import de timeline](tests/solidcraft_sketch_features.rs). Pour les autres formats de scripts SolidCraft, voir aussi `import-solidcraft` et ses sous-ensembles documentés dans les [résultats de campagne](docs/campagne-20261125.md).

## Familles de construction et portée

Les recettes JSON versionnées expriment les opérations et dimensions; l’audit reconstruit et vérifie les cartes candidates avant publication. Les familles implémentées couvrent notamment les paires affines, couronnes, équerres, perçages axiaux et mobiles, contours NURBS multiboucles, esquisses NURBS croisées, et primitives rectangle/cercle/ellipse extrudées ou tournées. Chaque famille possède des hypothèses d’admission propres. Une réussite sur un gabarit ne constitue pas une garantie universelle pour toute géométrie.

Le préplan de `plan-sketch` traite les jonctions, subdivisions et événements d’intersection nécessaires à des familles précises. Les certificats de Bernstein portent sur les domaines et méthodes annoncés par le rapport. Les classifications, grilles et échantillonnages restent des diagnostics mesurés; ils ne sont pas présentés comme des preuves de couverture ou d’injectivité. Les refus et contre-exemples sont conservés dans `tests/counterexamples`, `tests/regressions` et les rapports versionnés.

Le catalogue des primitives décrit leurs conventions et les résultats rejoués : [gabarits NURBS](docs/gabarits-primitives-nurbs.md). Pour les courbes importées et les opérations croisées, voir [prépartition des esquisses](docs/prepartition-esquisse-nurbs.md) et [cycle NURBS multiboucles](docs/cycle-nurbs-multiboucles.md). Ces documents distinguent les mesures de chaque recette des propriétés démontrées pour son domaine. Le CdC corrigé demeure la spécification normative; les notes d’adaptation et les journaux de campagne apportent contexte, recettes et résultats, sans modifier les exigences du CdC.

La démo native inclut aussi une poche prismatique versionnée :

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" build tests/references/sketch_boolean_blind_pocket.recipe.json build/demo/blind-pocket.json
"$cli" audit build/demo/blind-pocket.json
```

Cette fixture produit 9 blocs au total et une borne certifiée du déterminant de 191,99999999999991 mm³. La recette du test de référence plus large (`sketch_boolean_oracle.py`) utilise un domaine et un historique d’opérations différents; ses nombres ne se substituent pas à ceux de cette démo.

## Autres parcours

### Calcul mécanique

Le solveur Rust prend un modèle audité et un job versionné. Cet exemple reprend le job de l’équerre de référence, résout puis exporte champ et contraintes pour ParaView :

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" import-solidcraft upstream/fieldcraft/tests/oracle/filleted-bracket-fillet-stress/design.json build/demo/bracket.json
python3 -c 'import json; print(json.dumps(json.load(open("tests/references/mechanical_jobs/bracket.json"))["job"]))' > build/demo/job.json
"$cli" solve-job build/demo/bracket.json build/demo/job.json build/demo/result.json
"$cli" export-field-vtu build/demo/result.json 4 build/demo/bracket-field.vtu
```

Le modèle FieldCraft importé fournit 36 blocs et 55 interfaces. Le résultat conserve déplacements en mm, réactions, bilan travail-énergie et contraintes par bloc aux interfaces. Le solveur traite l’élasticité linéaire isotrope en petites déformations. Les comparaisons et limites de chaque expérience sont consignées dans le [journal de campagne](docs/campagne-20261125.md).

### Éditer la profondeur d’un perçage

La recette polygonale `closed` du corpus sert à rejouer une poche axiale et l’audit de son édition :

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
python3 -c 'import json; print(json.dumps(next(r["recipe"] for r in json.load(open("tests/references/interval_bores_corpus.json")) if r["id"]=="closed")))' > build/demo/cavity.recipe.json
"$cli" build build/demo/cavity.recipe.json build/demo/cavity.json
cat > build/demo/edits.json <<'JSON'
{"schema":"interval_bore_edits_v1","steps":[{"bore":"pilot","z_interval_mm":[0,9],"require_features":["bore/pilot/upper"]}]}
JSON
"$cli" edit-bore-intervals build/demo/cavity.json build/demo/edits.json build/demo/blind.json
"$cli" audit build/demo/blind.json
"$cli" export-step build/demo/blind.json build/demo/blind.step
```

Une sélection exigée qui disparaît fait échouer l’édition avant publication de la cible. Voir [l’oracle indépendant des éditions](tests/references/interval_bore_edits.py).

### Référence Python et évaluation

Les références Python historiques sont rejouables séparément de la crate Rust :

```sh
python3 -m trivariate capabilities
python3 -m trivariate build tests/references/affine_pair.json build/affine_pair.json
python3 -m trivariate audit build/affine_pair.json
python3 -m trivariate evaluate tests/references/quarter_annulus.json 0.5 0.5 0.5
```

`make bootstrap` récupère les six dépendances d’oracle consignées dans `references.lock` sur une nouvelle copie; après cette étape, `make check` fonctionne hors ligne. Pour les seuls tests Rust : `cargo test --locked --offline`.

## Documents et résultats

- [Contrats et architecture](docs/architecture.md)
- [Évaluateur L1 et résultats](docs/L1.md)
- [Critères de développement de Luna](docs/LUNA.md)
- [Consignes de développement](AGENTS.md)
- [Cahier des charges scientifique corrigé](CdC_noyau_trivarie_sans_trim_corrige.md)
- [Construction de l’équerre de référence](LienSolidecraft.md)
- [Adaptation et organisation des revues](ArtcratAdaptation.md)
- [Journal de campagne et résultats mesurés](docs/campagne-20261125.md)
- [Galerie interactive](https://jeffwitz.github.io/SolidCraftTrivariate/)
