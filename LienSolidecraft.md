**J’ai retrouvé l’exemple, et je pense qu’il fournit une base concrète pour démarrer notre noyau.** Il s’agit de **`bracket-hex`**, une équerre avec un congé et deux perçages, construite dans SolidCraft, partitionnée en **30 cellules**, puis maillée par FieldCraft.

Ton idée permet de fixer un premier périmètre réaliste : **construire des solides comme des assemblages de blocs paramétrables sur un cube**. Pour obtenir notre noyau trivarié, le changement central serait de **conserver ces blocs comme géométrie native, avec leurs cartes splines et leurs interfaces, dès la construction de la pièce**.

La lecture du code confirme plusieurs briques utiles. Elle montre aussi précisément où intervenir : FieldCraft produit aujourd’hui des maillages à partir de ces blocs ; SolidCraft conserve des corps B-rep. **Nous devons introduire entre les features et leurs représentations dérivées un véritable modèle volumique trivarié.**

## 1. L’exemple exact et ce qu’il démontre

### L’équerre partitionnée en 30 blocs

Voici les sources directement utiles :

- [La recette de construction et de partition](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket-blocks/recipe.json).
- [L’image du découpage en 30 blocs](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/docs/media/bracket-partition.png).
- [Le code qui prépare le calcul `bracket-hex`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/engine/src/sample.rs#L556-L580).
- [Le mailleur multibloc `hex_mesh_cells`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L1-L138).

La pièce est un profil en L de **60 × 50 mm**, avec des jambes de **8 mm**, extrudé sur **40 mm**. Elle comporte un congé intérieur de **rayon 6 mm**, puis deux perçages de **diamètre 8 mm**, un dans chaque jambe. Les trous ont des axes différents : la pièce terminée dépasse donc le cas d’une simple section extrudée dans une seule direction.

La recette prescrit **29 opérations `solid.split_body`** :

| Étape | Partition obtenue |
|---|---|
| Coupes à \(z=10\) et \(z=30\) | Trois tranches dans la largeur |
| Coupes dans les jambes et coupe diagonale dans le coude | Huit cellules dans chaque tranche |
| Découpage autour de chacun des deux trous, dans la tranche centrale | Chaque cellule contenant un trou devient quatre cellules disposées autour du trou |

On arrive à :

\[
8+(8+3+3)+8=30\ \text{cellules}.
\]

**Cette partition est explicitement écrite dans la recette.** Les plans de coupe, leur ordre et les noms successifs des corps sont fournis. C’est précisément ce qui rend cet exemple intéressant pour commencer : les décisions de construction des blocs sont déjà disponibles et examinables.

Dans l’application, l’exemple s’ouvre avec :

```sh
fieldcraft --sample bracket-hex
```

La documentation décrit ensuite un calcul sur **5 778 éléments hex20**, avec une taille demandée de 3 mm. Il faut bien distinguer les **30 blocs de partition** des milliers d’éléments qui les discrétisent. [FieldCraft](https://tangency-labs.github.io/fieldcraft/examples.html?utm_source=chatgpt.com)

### Ce que FieldCraft fait avec cette partition

Le code de `cells.rs` suit une organisation qui nous intéresse directement :

1. Reconnaître chaque cellule comme un **cube logique**.
2. Identifier les sommets, arêtes et faces partagés.
3. Affecter des nombres de subdivisions compatibles entre les blocs.
4. Construire chaque arête une seule fois.
5. Construire chaque grille de face une seule fois.
6. Remplir chaque intérieur par **interpolation transfinie trivariée**.
7. Contrôler le maillage obtenu.

Une face commune possède donc une grille commune, utilisée par les deux cellules. Une fois la topologie établie, les nœuds sont partagés par identité. La fixture du corpus enregistre **30 corps et 45 interfaces**. [Sources : algorithme multibloc](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L1-L32), [fixture](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket-blocks/fixture.json).

**C’est cette organisation que je reprendrais : une topologie commune, puis une géométrie commune sur les interfaces.**

### Il existe aussi un exemple de décomposition automatique

Le second exemple à étudier est **`crank-arm`**, une manivelle. Après une partition suivant son plan de symétrie, `mesh.hex.generate` décompose automatiquement sa section en régions adaptées au maillage structuré.

Le code correspondant est [`hex/decompose.rs`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/decompose.rs#L1-L48). Il sépare les trous, propose des découpages autour d’eux et cherche des coupes depuis les coins rentrants et les congés.

Deux limites de cet exemple orientent notre stratégie :

- ses partitions restent **internes au mailleur** et ne sont pas conservées dans le modèle géométrique ;
- la variante maillée est préparée pour rester balayable : certains congés et chanfreins sont retirés, et certaines rencontres géométriques sont modifiées pour éviter des dégénérescences. [Source : description de `crank-arm`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/docs/examples.md#L142-L165).

Je commencerais donc par **la partition maîtrisée de l’équerre**, puis utiliserais la manivelle pour étendre l’automatisation.

## 2. Ta contrainte est-elle suffisante ?

**Elle est suffisante pour définir le périmètre d’un premier noyau.** Je la formulerais cependant plus précisément :

> Le noyau construit des solides constitués de blocs trivariés réguliers, dont les interfaces sont explicitement partagées et paramétriquement compatibles.

« Décomposable en solides maillables en hexas » constitue un excellent point de départ, mais il reste quatre propriétés à établir :

| Propriété | Ce qu’elle apporte à notre noyau |
|---|---|
| Chaque bloc possède une topologie de cube | Un domaine paramétrique à trois coordonnées |
| Chaque bloc possède une carte spline régulière et injective | Une géométrie volumique exploitable pour l’analyse |
| Deux blocs voisins utilisent des traces compatibles | Un raccord géométrique et un assemblage des champs bien définis |
| Les blocs forment une partition du solide | Une couverture complète, sans recouvrement intérieur |

La distinction essentielle est entre **la topologie du bloc** et **sa paramétrisation**. Six faces logiques donnent une structure de construction ; elles ne garantissent pas à elles seules la qualité d’une carte \(X(u,v,w)\).

L’exemple contient d’ailleurs un détail instructif : ses 30 corps totalisent **188 faces B-rep**, alors que 30 cubes logiques possèdent 180 côtés. FieldCraft sait regrouper plusieurs faces géométriques en une face logique composite, par exemple au voisinage d’un raccord tangent entre un plan et un congé. [Sources : fixture](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket-blocks/fixture.json), [construction des faces composites](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L9-L15).

Pour nous, cela implique qu’une face de bloc pourra être une surface spline comportant plusieurs intervalles de nœuds. Et, si une cellule logique exige une subdivision supplémentaire pour obtenir des cartes régulières, le noyau devra pouvoir la produire explicitement.

**Je retiendrais donc le plan à 30 cellules comme première cible de construction, avec la régularité des cartes comme critère d’acceptation.** Le nombre 30 reste un objectif à vérifier, pas une contrainte qui justifierait de conserver une carte défectueuse.

## 3. La modification décisive dans FieldCraft : conserver les cartes

### L’état actuel

Les types de sortie sont très explicites :

- `HexMesh` contient des nœuds, une connectivité, des faces frontières et un rapport ;
- `CellsMesh.blocks` contient les nombres de divisions suivant les trois directions de chaque bloc.

Il n’y a pas, dans ces objets, de réseau de contrôle volumique, de poids NURBS ni de vecteurs de nœuds. [Sources : `HexMesh`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/mod.rs#L213-L223), [`CellsMesh`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L52-L72).

FieldCraft possède bien des cartes continues pour ses éléments finis et des évaluateurs analytiques pour certaines constructions. **Ce qui manque est la conservation des blocs de conception comme volumes splines indépendants du maillage.**

La fonction [`transfinite_volume`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/quad/tfi.rs#L204-L240), par exemple, remplit un tableau de points à partir de grilles de frontières. Elle constitue une référence utile pour les formules, mais son résultat reste une grille.

### L’objet à introduire

Je proposerais trois objets distincts :

| Objet proposé | Responsabilité |
|---|---|
| `BlockPlan` | Recettes de construction des blocs, incidences, orientations et provenance |
| `TrivariateComplex` | Cartes splines effectivement construites, interfaces canoniques, frontière et validité |
| `AnalysisSpace` | Bases d’analyse, raffinement, degrés de liberté, contraintes et quadrature |

Le cœur géométrique contiendrait notamment :

```rust
// Schéma proposé : ces types n'existent pas encore dans les projets.
struct Trivariate {
    degree: [usize; 3],
    knots: [Vec<f64>; 3],
    homogeneous_control_points: TensorGrid<[f64; 4]>,
}

struct NativeSolid {
    blocks: Vec<Trivariate>,
    interfaces: Vec<SharedInterface>,
    design_faces: Vec<DesignFace>,
    provenance: ConstructionHistory,
    validity: ValidityReport,
}
```

Le point crucial est que **les divisions demandées pour le calcul ne définissent plus la géométrie**. On construit d’abord les cartes. On peut ensuite les raffiner exactement, les échantillonner pour l’affichage ou en dériver un maillage EF.

### Ce qui se reprend directement dans l’organisation

Le mailleur multibloc possède déjà les bonnes idées pour :

- porter une arête commune ;
- porter une face commune ;
- convertir les orientations d’une face entre deux blocs ;
- propager des contraintes de subdivision.

Par exemple, [`FaceGrid::map_from`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L939-L971) gère les correspondances d’orientation des grilles. Nous pouvons reprendre cette logique pour les directions paramétriques.

En revanche, la compatibilité spline doit porter sur **les degrés, les nœuds, leurs multiplicités, les poids et les coefficients de trace**. Avoir le même nombre de points de chaque côté ne suffit pas.

**Pour le premier démonstrateur, je choisirais des raccords par faces entières, avec traces identiques.** C’est cohérent avec le mode multibloc de FieldCraft, qui refuse déjà les cellules se rencontrant sans partager leurs faces entières. Les raccords emboîtés prévus dans notre CdC peuvent venir ensuite. [Source : contrat de `hex_mesh_cells`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L100-L103).

## 4. Comment construire l’équerre directement en trivarié

Je transformerais la recette existante en **recette paramétrique de blocs natifs**.

Aujourd’hui, elle construit un corps, le perce, puis le découpe. Dans notre noyau, les features devront produire ou reconstruire un ensemble de blocs dont la réunion représente leur résultat.

### Une première palette d’opérateurs très limitée

Pour cette équerre, je viserais les constructions suivantes :

| Opérateur natif | Utilisation |
|---|---|
| Quadrilatère plan et extrusion | Zones droites des jambes |
| Quadrilatère à frontière circulaire et extrusion | Zones autour des perçages |
| Section avec raccord circulaire et extrusion | Coude avec congé |
| Découpage isoparamétrique | Subdivision exacte d’un bloc existant |
| Assemblage par faces partagées | Constitution du solide |
| Reconstruction d’une famille de blocs | Mise à jour d’un trou, d’une épaisseur ou d’un rayon |

L’extrusion est particulièrement favorable. Pour une section rationnelle \(\mathbf S(u,v)\), on construit :

\[
\mathbf X(u,v,w)=\mathbf S(u,v)+wh\,\mathbf d.
\]

Si les points de contrôle de la section sont \(\mathbf P_{ij}\), avec poids \(\omega_{ij}\), les deux couches de contrôle volumique sont :

\[
\mathbf P_{ij0}=\mathbf P_{ij},
\qquad
\mathbf P_{ij1}=\mathbf P_{ij}+h\mathbf d,
\]

avec les mêmes poids sur les deux couches et un degré 1 dans la direction d’extrusion.

**On obtient immédiatement un volume NURBS natif.** Sa définition ne dépend d’aucun nombre d’éléments.

### Les trous deviennent des constructions de blocs autour d’un vide

Pour chaque perçage, la partition de FieldCraft fournit quatre régions autour du trou. Nous pouvons nous en servir comme patron :

1. Définir exactement les arcs du trou.
2. Définir les frontières extérieures de la zone réservée au perçage.
3. Construire les quatre surfaces de section.
4. Les extruder à travers l’épaisseur de la jambe.
5. Partager leurs interfaces avec les blocs voisins.

L’opération « créer un trou » devient alors, dans cette famille admissible, **une reconstruction volumique locale avec une topologie connue**.

Cela répond directement à ton objectif : l’intérieur est construit et conservé dans le noyau. La surface cylindrique du trou est une restriction de cartes volumiques.

Il faudra cependant paramétrer ensemble les faces extérieures de ces blocs et celles de leurs voisins. Même quand une interface est géométriquement plane, ses deux paramétrisations peuvent différer. La face canonique partagée doit donc être une donnée de construction.

### Les arcs de FieldCraft demandent une adaptation réelle

Dans [`mesh/block.rs`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/block.rs#L177-L188), les arcs sont évalués par sinus et cosinus d’un paramètre angulaire.

Nous construirions les mêmes arcs géométriques avec des segments rationnels quadratiques. Leur paramétrisation change : les correspondances de bord doivent être recalculées en conséquence.

Pour les blocs plus généraux, les formules transfinies restent utiles, mais il faut construire des coefficients splines avec des frontières compatibles, puis contrôler les poids et le Jacobien. **Une TFI appliquée à des points projetés ne devient pas automatiquement une carte NURBS exacte.**

### Le congé doit commencer comme une famille qualifiée

Le congé intérieur de cette équerre provient d’une section plane arrondie, extrudée. C’est une bonne première opération native.

Je commencerais donc par :

> Un raccord circulaire dans une section, associé à une partition volumique définie et à des conditions de validité explicites.

L’extension aux congés spatiaux entre surfaces quelconques viendrait après. Cette progression permet de livrer des opérations de conception utiles tout en gardant des constructions volumiques maîtrisées.

## 5. Où remplacer Solidcore dans SolidCraft

### La façade existe, mais le remplacement n’est pas transparent

SolidCraft possède une crate `solidcraft-kernel` qui regroupe les opérations usuelles. Toutefois, son type `Body` repose explicitement sur un record B-rep et expose notamment :

```rust
pub fn record(&self) -> &Brep
pub fn smlib(&self) -> Result<Arc<solidcraft_smlib::Brep>>
```

Des couches supérieures utilisent aussi des types de Solidcore. Le changement porte donc sur le modèle géométrique et ses consommateurs, au-delà de la substitution d’une dépendance Cargo. [Sources : définition de `Body`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/kernel/src/body.rs#L165-L183), [accès publics](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/kernel/src/body.rs#L799-L808).

### Le bon point d’entrée est l’évaluation des features

SolidCraft possède déjà :

- des features sérialisées, dont esquisses, extrusions, révolutions et congés ;
- des paramètres et expressions ;
- des identifiants de features ;
- une timeline avec références ;
- un dispatch d’évaluation.

Le point d’intégration se trouve autour de `FeatureKind` et de `eval_feature_kind`. C’est là que l’on connaît encore **l’intention de construction**, les profils et les paramètres. [Sources : définition des features](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/doc/src/document.rs#L887-L934), [évaluation](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/doc/src/eval.rs#L14094-L14131).

Je créerais à ce niveau un contrat restreint :

> Une opération de conception résolue reçoit des paramètres et des références sémantiques, puis retourne un solide trivarié valide ou un refus expliqué.

Le premier ensemble d’opérations serait celui nécessaire à l’équerre. Cela permettrait de conserver progressivement l’éditeur et son historique tout en donnant au corps natif sa propre représentation.

### Une découverte utile : Solidcore contient des types de volumes

Solidcore contient effectivement un type `BSplineVolume`. Mais le fichier explique que le port suit un build SMLib sans la fonctionnalité `USE_VOLUMES` : l’évaluation du volume n’est pas disponible. La méthode `evaluate_simple` retourne explicitement `BSplineVolumeUnsupported`. [Source : implémentation de `BSplineVolume`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/surface/b_spline_volume.rs#L31-L38), [évaluation](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/surface/b_spline_volume.rs#L307-L324).

Nous avons donc des structures à examiner, mais pas un noyau volumique opérationnel à activer.

En revanche, des routines élémentaires existent réellement : évaluation des bases B-splines, dérivées de courbes et raffinement de réseaux de contrôle. Elles méritent d’être qualifiées et éventuellement reprises derrière notre API. [Sources : bases](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/nurbs/knots.rs#L444-L512), [raffinement](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/nurbs/crv_adv/p2.rs#L440-L479).

**Remplacer Solidcore comme modèle central reste compatible avec le réemploi de certaines de ses briques mathématiques.**

### Conserver le nommage et le rendu, en adaptant leur identité

Le nommage de SolidCraft porte déjà la provenance : feature, rôle de la face, subdivision, copie, etc. C’est une base intéressante. [Source : noms persistants](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/doc/src/persistent_name.rs#L1-L10).

Je garderais impérativement notre distinction :

- **face de conception** : « paroi du perçage », « face chargée », « face d’appui » ;
- **face de bloc** : frontière paramétrique d’une carte ;
- **interface interne** : raccord entre deux cartes.

La paroi d’un trou peut réunir quatre faces de blocs tout en conservant une seule identité de conception. Les conditions aux limites s’attachent à cette identité.

Le rendu est plus facile à raccorder : SolidCraft affiche déjà des maillages munis d’identifiants de faces. On peut produire ces triangles depuis la frontière du complexe trivarié et conserver une grande partie du moteur d’affichage. [Source : `SceneMesh`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/viewport/src/scene.rs#L24-L40).

## 6. Comment utiliser FieldCraft pour l’IGA

FieldCraft constitue un hôte intéressant parce que ses dépendances géométriques sont concentrées dans `fieldcraft-geometry`, tandis que le mailleur est séparé. [Sources : dépendances géométriques](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/geometry/Cargo.toml#L11-L29), [dépendances du mailleur](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/Cargo.toml#L11-L17).

Il faut toutefois ajouter une sortie volumique native à cet adaptateur. Sa sortie actuelle expose essentiellement des surfaces triangulées, des courbes discrétisées et certaines surfaces analytiques ; elle ne transporte pas des cartes NURBS générales.

### Le calcul IGA demande aussi un raccordement

Le solveur actuel est organisé autour de familles d’éléments finis. Ses fonctions de forme sont sélectionnées par type d’élément et certaines structures fixent le nombre maximal de fonctions locales à 27. Une base tensorielle cubique en utilise déjà 64 par élément. [Sources : `ElementType`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/fem/src/element.rs#L57-L74), [`MAX_NODES`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/fem/src/shape.rs#L24-L34).

Je passerais par une **extraction de Bézier** munie des poids, des opérateurs d’extraction et de la connectivité des degrés de liberté de contrôle.

Cela permettrait de réutiliser une partie de la mécanique et de l’algèbre linéaire, en adaptant :

- l’évaluation des bases et du Jacobien ;
- les tailles des supports locaux ;
- l’assemblage des interfaces ;
- les conditions aux limites ;
- l’évaluation des résultats aux points physiques.

L’espace d’analyse pourrait commencer avec des raccords \(C^0\) entre blocs et une continuité spline à l’intérieur des blocs. Notre premier résultat serait déjà une chaîne complète de conception et d’analyse volumiques natives.

### Les contrôles existants sont utiles, mais changent d’objet

FieldCraft contrôle les Jacobiennes de ses éléments par représentation de Bernstein et subdivision. Cette démarche est pertinente pour nous. Cependant, elle porte sur les cartes polynomiales des éléments hex8/hex20 et utilise des coefficients en `f64`. [Sources : contrôle du Jacobien](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/quality/jacobian.rs#L1-L15), [coefficients Bernstein](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/quality/bernstein.rs#L42-L50).

Nous devons étendre cette logique aux cartes rationnelles, avec l’encadrement numérique prévu dans notre CdC, et contrôler séparément la validité globale de la partition.

**La qualité des hexas échantillonnés restera un indicateur utile ; la validité géométrique sera portée par les cartes natives elles-mêmes.**

## 7. La stratégie de démarrage que je recommande

Je construirais **un noyau autonome restreint**, puis deux adaptateurs vers SolidCraft et FieldCraft. Cela permet de faire avancer le cœur géométrique sans devoir immédiatement reprendre toute l’application CAO.

### Étape 1 — Une construction native élémentaire

Premier résultat :

- une section rationnelle à quatre côtés ;
- son extrusion ;
- l’évaluation de \(X\) et \(DX\) ;
- les traces de frontière ;
- la subdivision exacte ;
- le contrôle de validité ;
- la sauvegarde et la relecture.

Puis un assemblage de quelques blocs partageant leurs faces.

**Critère de réussite :** changer la résolution d’affichage ou d’analyse conserve la géométrie.

### Étape 2 — L’équerre de FieldCraft comme démonstrateur

Reprendre ses dimensions, ses perçages, son congé et son plan de partition, puis exprimer les blocs directement avec nos constructeurs.

Je remplacerais les noms de corps issus des coupes successives par des identités explicites : bloc de jambe, bloc de coude, secteur de perçage, tranche, etc.

La recette permet déjà un premier contrôle indépendant du volume. À partir de ses dimensions :

\[
\begin{aligned}
V
&=8(60+50-8)\,40\\
&\quad+6^2\left(1-\frac{\pi}{4}\right)40
-2\pi\,4^2\,8\\
&=34080-616\pi\\
&\approx 32144{,}7789254\ \mathrm{mm}^3.
\end{aligned}
\]

J’ai recalculé cette valeur : elle correspond à celle enregistrée dans la recette. Ce contrôle ne démontre pas à lui seul l’égalité de forme, mais fournit une référence analytique utile.

**Critère de réussite :** reproduire la pièce avec des cartes natives, des interfaces explicites et une frontière conforme aux primitives prescrites.

### Étape 3 — La modification paramétrique

Faire varier successivement :

- l’épaisseur des jambes ;
- la largeur extrudée ;
- les diamètres et positions des trous ;
- le rayon du congé.

Le noyau doit déterminer si le plan de blocs reste valide, nécessite une reconstruction ou sort du domaine admissible.

L’historique prend ici une fonction précise : il porte les paramètres, les dépendances et la provenance des entités. La partition devient une construction dérivée explicite, avec ses propres conditions de validité.

### Étape 4 — Le calcul IGA sur le même modèle

Reprendre le chargement de l’exemple : acier, encastrement en bout de jambe et effort de 1 kN sur l’autre extrémité.

Comparer les déplacements, l’énergie, les réactions et des contraintes évaluées aux mêmes positions physiques. Les valeurs de FieldCraft serviront de comparaison initiale, accompagnées des études de convergence.

Le README évoque un accord à 0,1 %, mais cela concerne le déplacement de pointe documenté. Les tests inspectés utilisent d’autres seuils et ne recalculent pas tous la référence tétraédrique. Je traiterais donc ces chiffres comme des résultats publiés par le projet à reproduire, avec leur périmètre précis. [Source : test de convergence de l’équerre](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/testkit/tests/bracket_blocks.rs).

### Étape 5 — Automatiser la proposition de partition

À ce stade, `hex/decompose.rs` devient une base utile pour proposer des plans de blocs :

- séparation des trous ;
- secteurs autour des perçages ;
- coupes depuis les coins rentrants ;
- prise en compte des congés ;
- classement de partitions candidates.

Il faudra adapter son critère de décision : aujourd’hui, son évaluation dépend notamment de la taille de maillage demandée et peut conserver une solution de maillage libre. Pour notre constructeur, il devra proposer une partition admissible pour des cartes splines, puis la conserver indépendamment des raffinements d’analyse. [Source : principes et évaluation du décomposeur](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/decompose.rs#L19-L48).

## 8. Ce que cette étude change dans notre proposition

Elle nous donne une répartition concrète du travail :

| Base existante | Réemploi envisagé | Travail propre au projet |
|---|---|---|
| SolidCraft | Esquisses, paramètres, features, commandes, provenance | Évaluation native des features et nouveau type de corps |
| FieldCraft multibloc | Incidences, orientations, contraintes de subdivision, exemples | Plans persistants et interfaces de cartes splines |
| Constructions transfinies | Formules et référence d’échantillonnage | Construction rationnelle et contrôle des cartes |
| Solidcore, partie NURBS | Certaines routines élémentaires après qualification | Complexe volumique, constructeurs et opérations natives |
| FieldCraft calcul | Infrastructure mécanique, unités, solveurs et visualisation | Bases IGA, extraction, assemblage et post-traitement adaptés |

Un point pratique concerne l’interface graphique : les composants `ui-kit`, `look` et `theme` de SolidCraft ont une licence distincte du code géométrique. Le réemploi du moteur et celui de l’apparence doivent donc être examinés séparément au moment de définir le produit distribué. [Source : licence de l’interface](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/LICENSE-BRAND).

**Ma recommandation est de commencer par l’équerre à 30 cellules et de faire de son plan de blocs un objet géométrique persistant.** Nous aurions ainsi un premier noyau capable de construire, modifier, afficher et analyser une pièce percée et arrondie en restant trivarié.

L’étude a porté sur le code de FieldCraft au commit `35e8346`, son pin SolidCraft v0.42.0, et l’architecture plus récente de SolidCraft au commit `765f26e`. J’ai contrôlé les sources, la recette, les structures de données et les assertions des tests ; je n’ai pas exécuté le mailleur ou le solveur, Rust et Cargo étant absents de cet environnement.

**Le premier résultat à viser est donc précis : cette équerre construite directement par nos cartes rationnelles, avec ses faces partagées et un calcul IGA sur le même modèle.** C’est un jalon suffisamment riche pour vérifier l’architecture, tout en gardant un ensemble d’opérations que l’on peut réellement maîtriser.
