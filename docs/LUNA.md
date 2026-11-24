# Feuille de route de Luna

La campagne développe le noyau et ses audits en Rust. Le lot L1 est porté dans `src/bezier.rs` et ses tests natifs sont dans `tests/bezier.rs`. Rejouer `cargo test --locked --offline` avant les contrôles croisés de `make check`. Le descriptif Python ci-dessous documente le point de départ historique.

## L1 — Évaluer une carte rationnelle native

**Point de départ.** Rejouer `make check`, consigner `git rev-parse HEAD` et inspecter `build/affine_pair.json`. Première modification proposée dans `trivariate/spline.py`, avec tests et recettes sous `tests/`. Conserver l’API et les résultats du cas affine.

**Livrable.** Une carte de Bézier tensorielle rationnelle à degrés explicites, supportant (1,1,1) et (2,1,1), avec évaluation de X et DX sur le cube fermé. Fixer la convention des colonnes du jacobien, l’ordre des coefficients, les unités et les paramètres de face avant l’implémentation. Distinguer une entrée invalide d’une famille encore hors domaine.

**Recette indépendante.** Extruder suivant z un quart de couronne de rayons 1 et 2 mm sur 3 mm. Pour chaque rayon r, les contrôles euclidiens de l’arc sont (r,0), (r,r), (0,r), de poids 1, √2/2, 1. Choisir l’ordre des axes pour un déterminant positif, par exemple angle croissant en u, rayon décroissant en v, altitude croissante en w. Enregistrer comment le poids irrationnel est représenté.

**Acceptation à figer pour ce lot.** Vérifier le cas affine sur sommets, faces et intérieur. Pour le cas courbe, vérifier les extrémités, le milieu angulaire, le rayon et z sur une grille incluant les bords (erreur absolue ≤ 10⁻¹² mm en double précision pour ces dimensions). Comparer les colonnes de DX à des dérivées indépendantes, préciser la méthode et son seuil avant essai. Tester les dimensions incohérentes et les poids nuls. Rapporter séparément les tests ponctuels et toute preuve sur le domaine continu.

## Suite du démonstrateur

1. **Deux blocs courbes.** Découper le quart de couronne suivant z, partager la face canonique, vérifier les permutations et inversions, sauvegarder/recharger, modifier la hauteur et tester un refus atomique. Contrôler la subdivision exacte et la reproduction d’un champ affine. Qualifier bornes de poids, déterminant, injectivité et partition pour annoncer J1 sur cette famille.
2. **Équerre.** Construire par gabarits les jambes, le congé de rayon 6 et les deux trous de diamètre 8. Comparer à la fixture FieldCraft (30 cellules, 45 interfaces) et au volume analytique 34080 − 616π mm³. Justifier les subdivisions supplémentaires. Le volume seul ne démontre pas l’égalité de forme.
3. **J2.** Tester conformité H¹, modes rigides, champ affine, équilibre et convergence sous raffinement, puis comparer à la référence EF.
4. **J3.** Régénération par dépendances et sensibilités vérifiées à régime fixé. Signaler les événements de topologie ou de sélection.
5. **J4.** Geler corpus et seuils, publier couverture avec refus, temps et mémoire à précision comparable.

Chaque livraison indique commandes, révision, paramètres, résultats, limites et prochain critère ouvert. Les jalons scientifiques restent ouverts au démarrage.
