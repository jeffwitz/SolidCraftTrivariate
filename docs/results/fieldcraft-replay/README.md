# Équerres FieldCraft reconstruites depuis les esquisses

[Télécharger les cinq modèles, leurs recettes, leurs STEP et les contrôles](./fieldcraft-native.zip)

L’équerre d’origine a deux ailes de 60 et 50 mm, une épaisseur de paroi de 8 mm, une extrusion de 40 mm et un congé intérieur de rayon 6 mm. Les quatre variantes de la campagne SolidCraft modifient seulement les diamètres des perçages. Le premier trou est centré en x = 40 mm et percé suivant Y. Le second est centré en y = 32 mm et percé suivant X. Tous deux sont centrés en z = 20 mm.

| Cas | Ø suivant Y (mm) | Ø suivant X (mm) | Volume analytique (mm³) |
|---|---:|---:|---:|
| x64-y64, original | 8 | 8 | 32 144,778925 |
| x40-y80 | 5 | 10 | 32 163,628481 |
| x80-y40 | 10 | 5 | 32 163,628481 |
| x50-y74 | 6,25 | 9,25 | 32 165,984676 |
| x51-y74 | 6,375 | 9,25 | 32 156,069024 |

Le constructeur du noyau Rust part du contour en L, du congé et des centres, axes et rayons des deux perçages. Pour les coordonnées fractionnaires, la recette utilise une grille en seizièmes de millimètre, puis applique une homothétie rationnelle exacte. Les deux zones de construction autour des trous partagent les mêmes niveaux en Z. Chaque résultat compte 48 blocs et 74 interfaces conformes. Ces nombres décrivent la partition obtenue. Sa minimalité reste ouverte.

FreeCAD/OCCT a relu chaque STEP comme un solide valide à une coquille. Les cinq STEP et les cinq références SolidCraft concordent aux 1 920 points de la grille et aux 80 sondes placées à ±0,02 mm des deux parois de perçage. Avec les propriétés de masse adaptatives, le plus grand écart de volume entre deux STEP correspondants est inférieur à 5 × 10⁻⁹ mm³.

Une fois l’archive décompressée, le CLI permet de reconstruire chaque solide avec `build`, de le contrôler avec `audit` et de l’exporter avec `export-step`. Le README de l’archive précise les entrées du contrôle mécanique sur la pièce d’origine.

[Contrôle indépendant des cinq exports](../fieldcraft-kernel-replay.md) · [Origine du modèle FieldCraft](https://github.com/tangency-labs/fieldcraft/tree/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket) · [Rejeu mécanique de l’équerre](../fieldcraft-bracket-pa5-mechanics.md)
