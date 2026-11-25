# Bilan de la campagne Luna du 25 novembre 2026

**Version de travail pendant la campagne.** La fenêtre prévue court du 24 novembre à 23 h 23 au 25 novembre à 09 h 23 UTC. Base L1 `ea5b940`, branche `campagne-20261125-10h`. Le journal conserve les protocoles datés, les seuils, les changements de méthode et la chronologie des échecs. [Journal de campagne](campagne-20261125.md)

## Ce que permet le démonstrateur

Un historique SolidCraft pris en charge peut produire un solide formé de cartes rationnelles sans découpe paramétrique. Le noyau Rust réaudite ces cartes, leurs raccords et la provenance des faces avant adoption. Le même modèle sert à exporter sa frontière en STEP, à poser une traction ou une pression sur une face nommée, à résoudre un problème d’élasticité et à exporter les champs en VTU. Le README donne une chaîne exécutable sur l’équerre de référence. [Démonstration](../README.md)

## Résultats par jalon

| Jalon | Essai et mesure | Portée |
|---|---|---|
| J1 — géométrie | Les sept géométries du corpus initial ont passé leur audit. L’équerre comporte 36 blocs, 55 interfaces, un volume de 32 144,7789254 mm³ et une borne continue de déterminant de 160 mm³. Les deux plans de tangence supplémentaires expliquent l’écart aux 30 blocs et 45 interfaces de la fixture amont. | L’audit de chaque famille précise la méthode de régularité, d’injectivité et de couverture. Une quadrature du volume ou une classification de points reste une observation. [Contrat de l’équerre](planar-certificate.md) |
| J2 — mécanique | Les réassemblages SciPy de deux équerres aux ordres 4, 6 et 8 retrouvent les six modes rigides du corps libre. Sur deux rayons, la distance énergétique native h3 à la référence C3D10 nominale 0,9 mm est voisine de 1,1 %. Le changement entre deux références EF fines est lui-même de l’ordre de 1,4 %. | Les comparaisons portent sur le déplacement et sa déformation en des points physiques communs. Le rayon, la résolution, le matériau et les appuis font partie du résultat. [Réassemblage](results/solver-package-20261125.json) · [Comparaison des champs](results/bracket-field-refinement-20261125.json) |
| J3 — modification | Les 1 681 couples de rayons de l’équerre ont passé l’audit. Le minimum de la borne certifiée de déterminant sur la grille est 151,4213562 mm³. Six couples hors domaine ont été refusés. Les dérivées de rayon, de carte et de charge sont contrôlées en régime de sélections fixes. | Pour les cinq déformations affines S2, dérivée tangente, adjoint et différences centrées concordent selon les critères principaux. Le test auxiliaire du reste quadratique admet 4/5 cas. [Grille exhaustive](results/bracket-radius-grid-20261125.json) · [Sensibilités de job](results/affine-job-derivative-20261125.json) |
| J4 — coût à précision fixée | Le corpus J4b comprend quatre couples et 76 essais. C3D10 h1,5 et l’immergé h1,25/q8 satisfont les critères de pointe et de congé sur 4/4. Le parcours natif uniforme obtient 0/4 au critère de stabilité de contrainte entre h2 et h3. Le parcours natif démarrant à h1 et raffiné autour du congé atteint ces critères sur 4/4. Les 24 répétitions J4r ont qualifié leurs observables. | Deux exécutions par méthode et par géométrie, affinité enregistrée, autres tâches de campagne actives. Le coût du niveau retenu, le cumul de recherche et le contrôle global du champ répondent à trois questions différentes. [Essais initiaux](results/bracket-cost-20261125.json) · [Répétitions](results/bracket-serial-cost-20261125.json) |

## Domaines géométriques et interopérabilité

- Les programmes de boîtes G16, les profils simples G17 et les profils à contours internes polygonaux G18 ont des contrats de construction distincts. G17 a obtenu les 156 décisions attendues (117 admissions, 39 refus), G18 les 152 décisions attendues (119 admissions, 33 refus). Les STEP de sept cas G17 et huit cas G18 ont été relus en CAO indépendante. [G16](results/orthogonal-program-20261125.json) · [G17](results/polygon-prism-20261125.json) · [G18](results/perforated-prism-20261125.json)
- L’image affine rationnelle G19 transporte cartes, déterminants, frontières et sélections. Six solides ont été comparés aux oracles, dont un tore cisaillé. Ce dernier exige le mode homogène entier i64. Le rapport conserve le refus du mode initial et trois témoins pour lesquels la classification OCC de la primitive portant directement la transformation diverge. Une reconstruction indépendante de sa frontière torique tranche ces trois cas. [Images affines](results/affine-images-20261125.json)
- SolidCraft 0.48.1 a enregistré puis rouvert les documents éditables utilisés dans les comparaisons de l’équerre, des programmes de boîtes et des extrusions de contours. Les contrôles rapprochent le solide hôte et le STEP natif par masse, géométrie et classification. Les familles non reconnues déclenchent un refus explicite de l’adaptateur. [Équerre](results/solidcraft-bridge-20261125.json) · [Historiques de boîtes](results/solidcraft-box-timeline-20261125.json) · [Profils](results/solidcraft-polygon-20261125.json)

## Où les critères mécaniques restent ouverts

Deux constructions du même prisme peuvent avoir des maillages différents. Sur le L et la plaque à deux tunnels du lot M7, les 24 résolutions ont abouti. À h3, les écarts de travail valent 0,2853 % et 0,02412 %, ceux de déplacement L² 0,3352 % et 0,03341 %, et les distances énergétiques 10,338 % et 5,0234 %. Les deux distances énergétiques dépassent le seuil de 5 %. Les nuages indépendants totalisent respectivement 49 152 et 344 064 sites physiques, tous situés une fois par partition. [Partitions équivalentes](results/equivalent-prism-mechanics-20261125.json)

Le raffinement local M8 a mis en évidence une identification d’arête manquante, puis un montage partiel qui lie plusieurs coefficients libres. Sur le témoin d’arête, 7 272 équations indépendamment reconstruites sur 7 641 coefficients locaux ont un rang rationnel démontré de 5 193. La dimension montée est donc 2 448. Le job v3 élimine aussi les équations d’appuis dépendants en rationnels exacts. [Témoin de rang](results/adaptive-edge-rank-20261125.json) · [Premiers parcours M8](results/adaptive-equivalent-prisms-20261125.json)

Sur le cylindre sous pression proche de l’incompressibilité, la correction creuse 3D au niveau h4 a passé les six cas imposés. Les mesures distinguent le résidu des équations discrètes, le déplacement, la pression et l’erreur énergétique contre la solution de Lamé. La qualification du solveur discret et la convergence de la solution physique ont chacune leur critère. [Lamé et corrections](results/pressure-cylinder-full-correction-20261125.json)

## Reproduction et provenance

Depuis un clone préparé avec ses six références verrouillées, la vérification de base s’exécute ainsi.

```sh
cargo test --locked --offline
cargo clippy --locked --offline --all-targets -- -D warnings
make check
```

Les études longues ont leurs recettes sous `tests/references/`, leurs résultats numériques sous `docs/results/` et leurs données régénérables sous `build/`. Le noyau et les nouveaux constructeurs sont en Rust. Les scripts Python servent aux références de portage, aux oracles, aux relances de solveurs et aux mesures indépendantes. [Architecture](architecture.md) · [Références verrouillées](../references.lock)

Le test EF amont de FieldCraft à la révision `35e8346` requérait des objets et des dépendances Cargo absents du checkout local. La revue automatique d’autorisation a rejeté le téléchargement des dépendances. Les nombres publiés par cet essai restent identifiés comme valeurs amont. Les comparaisons exécutées ici utilisent la fixture de conception verrouillée, SolidCraft, Open CASCADE, Netgen et CalculiX. Le choix de licence du code original avant diffusion relève toujours du porteur du projet (NF-8).
