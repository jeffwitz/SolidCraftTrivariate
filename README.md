# SolidCraftTrivariate

Noyau volumique natif en Rust pour la conception paramétrique et l’IGA, avec des références de test en Python. Le cas de départ construit deux prismes unitaires en hauteur et profondeur, de largeurs réglables, et vérifie leurs coefficients homogènes sur toute la carte.

La galerie interactive permet de couper les modèles, d’examiner leurs blocs et de télécharger les recettes, les cartes rationnelles et les fichiers STEP. [Ouvrir la galerie](https://jeffwitz.github.io/SolidCraftTrivariate/)

## Démarrer avec Luna

Depuis la racine, avec Cargo, Rust, Python 3.10+, Git, Bash et Make :

```sh
make check
codex -m gpt-6-luna "Lis AGENTS.md, docs/L1.md et docs/LUNA.md. Rejoue make check puis poursuis le premier lot Rust ouvert dans le journal de campagne."
```

Sur une nouvelle copie du projet, lancer d’abord `make bootstrap` pour récupérer les six références de `references.lock`. Sur la machine préparée, `make check` fonctionne hors ligne. Les tests du noyau Rust seul se lancent avec `cargo test --locked --offline`. La commande `make check` ajoute les références Python et les démonstrateurs de comparaison.

## Construire, auditer et calculer en Rust

Le fichier de référence FieldCraft décrit une équerre de 60 × 50 mm, extrudée sur 40 mm avec un congé et deux perçages. Les commandes suivantes importent ses opérations SolidCraft, auditent les cartes rationnelles et exportent un solide STEP.

```sh
cargo build --release --locked --offline
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" import-solidcraft upstream/fieldcraft/tests/oracle/filleted-bracket-fillet-stress/design.json build/demo/bracket.json
"$cli" audit build/demo/bracket.json
"$cli" export-step build/demo/bracket.json build/demo/bracket.step
"$cli" export-geometry-vtu build/demo/bracket.json 6 build/demo/bracket-geometry.vtu
```

Le VTU géométrique contient les coordonnées paramétriques, les trois valeurs singulières du jacobien, son déterminant signé et son conditionnement à chaque point affiché. L’index de cellule renvoie au bloc dans le modèle JSON. Ces valeurs sont échantillonnées.

Le calcul suivant applique une traction de 1 000 N à l’équerre. Il utilise le modèle importé, un raffinement local et le fichier de conditions aux limites de l’oracle. L’extraction de l’entrée de calcul requiert Python 3.

```sh
python3 -c 'import json; print(json.dumps(json.load(open("tests/references/mechanical_jobs/bracket.json"))["job"]))' > build/demo/job.json
"$cli" solve-job build/demo/bracket.json build/demo/job.json build/demo/result.json
"$cli" export-field-vtu build/demo/result.json 4 build/demo/bracket-field.vtu
```

Le fichier VTU s’ouvre dans ParaView. Il conserve les valeurs de contrainte propres à chaque bloc aux interfaces. Le résultat JSON conserve les coefficients homogènes Q2, les contraintes cinématiques, les charges nodales, les déplacements en mm, les réactions et le bilan travail-énergie. Le solveur traite l’élasticité linéaire isotrope en petites déformations. Chaque famille possède son propre domaine géométrique et ses contrôles. Les mesures de convergence, les comparaisons aux solveurs indépendants et les cas refusés figurent dans le [journal de campagne](docs/campagne-20261125.md).

## Modifier la profondeur d’un perçage

Le document SolidCraft de boîte percée sous `tests/references/solidcraft-blind-bore/` accompagne un historique exporté par l’hôte. L’exemple suivant part de la recette polygonale du corpus pour suivre une cavité fermée qui débouche sur la face de départ. Les coordonnées axiales sont en millimètres.

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

Le fond supérieur reste sélectionnable par son nom. Le rapport d’édition indique la disparition de la face inférieure, le changement de volume en multiples de π et les couches réutilisées. Une séquence contenant une sélection exigée qui disparaît échoue avant publication du fichier cible. [Oracle des éditions](tests/references/interval_bore_edits.py)

## Vérifier le noyau Rust

```sh
cargo test --locked --offline --test bezier
```

Ce test évalue la position et le jacobien sur les 125 points du quart de couronne, avec une référence trigonométrique indépendante.

## Rejouer la construction Python de référence

```sh
python3 -m trivariate capabilities
python3 -m trivariate build tests/references/affine_pair.json build/affine_pair.json
python3 -m trivariate audit build/affine_pair.json
```

Modifier les deux largeurs entières positives dans une copie de la recette puis rejouer la construction. La sortie JSON contient les deux réseaux de contrôle et l’incidence orientée. Le rapport donne le déterminant minimal en mm³ et le volume en mm³. Une commande refusée renvoie un code de sortie non nul.

## Documents de travail

L’oracle Python de L1 accepte les degrés (1,1,1) et (2,1,1). Pour inspecter la position et les trois colonnes du jacobien au centre du quart de couronne, exécuter :

```sh
python3 -m trivariate evaluate tests/references/quarter_annulus.json 0.5 0.5 0.5
```

- [Contrat et résultats de l’évaluateur L1](docs/L1.md)
- [Contrats, conventions et architecture](docs/architecture.md)
- [Premier lot et critères d’acceptation](docs/LUNA.md)
- [Règles de développement](AGENTS.md)
- [Cahier des charges scientifique](CdC_noyau_trivarie_sans_trim_corrige.md)
- [Construction de l’équerre de référence](LienSolidecraft.md)
- [Organisation des preuves et des revues](ArtcratAdaptation.md)
