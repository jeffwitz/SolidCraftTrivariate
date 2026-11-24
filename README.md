# SolidCraftTrivariate

Noyau volumique natif en Rust pour la conception paramétrique et l’IGA, avec des références de test en Python. Le cas de départ construit deux prismes unitaires en hauteur et profondeur, de largeurs réglables, et vérifie leurs coefficients homogènes sur toute la carte.

## Démarrer avec Luna

Depuis la racine, avec Cargo, Rust, Python 3.10+, Git, Bash et Make :

```sh
make check
codex -m gpt-6-luna "Lis AGENTS.md, docs/L1.md et docs/LUNA.md. Rejoue make check puis poursuis le premier lot Rust ouvert dans le journal de campagne."
```

Sur une nouvelle copie du projet, lancer d’abord `make bootstrap` pour récupérer les six références de `references.lock`. Sur la machine préparée, `make check` fonctionne hors ligne. Les tests du noyau Rust seul se lancent avec `cargo test --locked --offline`. La commande `make check` ajoute les références Python et les démonstrateurs de comparaison.

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
