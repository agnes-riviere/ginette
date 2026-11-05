# Calage de lomos231 avec Ginette : note de passation

Point d'observation : `lomos231` (dispositif LOMOS-mini, 5 capteurs de température de 0 à 40 cm sous le lit, une mesure de différence de charge deltaP).
Tout se règle dans `config_lomos.py`. Les scripts se lancent depuis ce dossier.
Pour travailler sur ce point, mettre `POINT_NAME = "lomos231"` dans `config_lomos.py` (ou exporter `GINETTE_POINT_NAME=lomos231`) : la valeur par défaut du dépôt est `lomos230`.
Les données du point (`OBS_point/lomos231/donnees_plot_lomos231.csv`) ne sont pas dans le dépôt : à demander séparément et à placer dans ce dossier avant de lancer.

## 1. Ce qui est fixé, et pourquoi

| Élément | Valeur | Justification |
|---|---|---|
| Période comparée | 20/04/2026 13:00, 40 jours (+ 15 jours de spin-up avant) | Période recommandée par `stallman_diffusivity.py` : charge, température et méthode de Stallman y indiquent le même sens d'écoulement (infiltration). |
| Domaine | Temp1 (CL haute) à Temp4 (CL basse), 0,30 m | Temp2 et Temp3 servent de capteurs de calage. |
| Offset de deltaP | −2,65 cm | Déduit de Stallman, pas calé sur les températures : sur les deux fenêtres où les diffusivités d'amplitude et de phase s'accordent à 1 % (11–31/03 et 20/04–10/05), le flux vaut K·(dp + offset)/0,4 avec le même K. Cela donne −2,65 cm et K ≈ 7·10⁻⁶ m/s. Cet offset est une estimation à partir de deux fenêtres : son incertitude est de l'ordre du centimètre. |
| Densité du solide | 2650 kg/m³ | Milieu minéral. |
| Capacité calorifique du solide | 800 J/kg/K | Milieu minéral. Le jeu « sédiment organique » (1180 et 2000) a été testé : il donne le même misfit mais aucun point compatible avec la diffusivité de Stallman, il est donc écarté. |

## 2. Paramètres calés (grille, 1890 combinaisons)

| Paramètre | Plage | Pas |
|---|---|---|
| log_k (perméabilité intrinsèque, m²) | −13,0 à −11,0 (42 valeurs ; la grille déborde d'un pas, jusqu'à −10,95, à cause de l'arrondi des flottants) | 0,05 |
| lam (conductivité de la fraction solide, W/m/K) | 2,0 à 4,0 | 0,25 |
| n (porosité) | 0,20 à 0,40 | 0,05 |

## 3. Lancer

```
python3 0_boundary_conditions_real_case.py   # conditions aux limites (écrit GINETTE_SENSI/)
python3 1_define_grid_search.py              # grille results/lomos231/grid_search.csv
python3 2_run_real_case.py                   # simulations, puis misfit, puis figures (flags dans la config)
```

- Environ 6 minutes pour 1890 simulations avec 6 processus. Le nombre de processus se règle avec `GINETTE_MAX_WORKERS`.
- Après tout changement de point, de période, d'offset ou de bornes, relancer les trois étapes dans l'ordre : `GINETTE_SENSI/` et `grid_search.csv` ne se régénèrent pas seuls.
- Le nom du point doit être exactement `lomos231` (un nom absent des tables prend en silence les bornes de `lomos230`).

## 4. Résultat de référence

- Meilleur point : log_k = −11,3, lam = 2,0, n = 0,25, misfit_tot = 13,47. KGE ≈ 0,98 (Temp2) et 0,97 (Temp3).
- Péclet = −0,2 : infiltration, conduction encore dominante.
- **log_k a un minimum net dans la grille** : le meilleur misfit selon log_k vaut 21,7 à −13,0, 13,73 à −11,4, 13,47 à −11,3 et 15,78 à −11,0. Le plateau (à moins de 5 % du minimum, 89 points sur 1890) est entre −11,45 et −11,15, soit K de l'ordre de 4·10⁻⁵ m/s. Attention : à Pe ≈ 0,2 (régime conductif), le diagnostic de `3_misfit.py` indique que log_k n'est identifié que **conditionnellement à l'offset de deltaP**, qui fixe le flux.
- **lam et n ne sont pas identifiables** : pour n entre 0,20 et 0,40, le meilleur misfit varie de 13,47 à 13,54. Pour lam, il est de 13,47 à 13,50 entre 2,0 et 2,75, puis monte (14,08 à 3,5 et 14,75 à 4,0) : lam est borné vers le haut, pas vers le bas. Le meilleur lam (2,0) est à la borne basse de la grille uniquement parce que la courbe est plate, pas parce que l'optimum est plus bas.
- misfit_tot et mae_tot donnent le même log_k et le même lam, et ne diffèrent que sur n (0,25 et 0,35) : sans conséquence, puisque n n'est pas identifiable.
- Les capteurs ne mesurent que la diffusivité thermique Γ = λ_bulk / Cv (avec λ_bulk = 0,598·n + lam·(1−n)), pas lam et n séparément.

## 5. Ce qu'on peut annoncer

- Perméabilité intrinsèque k ≈ 5·10⁻¹² m² (log_k ≈ −11,3 ; plateau −11,45 à −11,15, soit k de 3,5 à 7,1·10⁻¹² m²), soit une conductivité hydraulique K ≈ 4·10⁻⁵ m/s (2,6 à 5,4·10⁻⁵ m/s) pour de l'eau à environ 10 °C (K = k·ρg/μ, μ = 1,3·10⁻³ Pa·s). C'est la gamme d'un sable fin à limoneux. Sous réserve de l'offset de deltaP.
- Comparaison avec Stallman (calcul indépendant, K = uz/gradient) : fenêtre 20/04–10/05, log_k = −13,1 ; fenêtre 10/05–30/05, log_k = −11,6 ; K = 6,8·10⁻⁶ m/s (log_k ≈ −12,0) pour les deux fenêtres fiables prises ensemble. Le calage Ginette (−11,3) est du côté de la fenêtre de mai ; l'écart avec Stallman reste à expliquer et fait partie de l'incertitude.
- lam entre 2,25 et 3,5 W/m/K et n entre 0,2 et 0,4 (20 points du plateau, sur 89), compatibles avec la diffusivité mesurée par Stallman (0,74 à 1,01·10⁻⁶ m²/s). Pour ces points, λ_bulk ≈ 2,0 à 2,3 W/m/K et Cv ≈ 2,5 à 2,9 MJ/m³/K.
- Ne pas présenter la valeur de lam ou de n au minimum du misfit comme un résultat.

## 6. Limites et pistes

- log_k dépend de l'offset de deltaP : avec l'offset à 0 sur la même période, on obtient log_k ≲ −12,5. Un balayage de l'offset (−2 à −3,5 cm, via la variable `GINETTE_DELTAP_OFFSET`, en mètres) donnerait l'incertitude sur log_k. Une mesure manuelle niveau rivière/piézomètre donnerait le vrai zéro du capteur.
- Les fenêtres de février à avril indiquent une exfiltration que ni deltaP brut ni l'offset retenu ne reproduisent : ne pas élargir la période sans le traiter.
- Pour séparer lam et n, il faut une information indépendante : porosité mesurée sur un échantillon de sédiment, conductivité mesurée, ou une période à Péclet plus grand.
- Les valeurs de conductivité et de porosité de la littérature citées plus haut sont données de mémoire et à vérifier avant publication.
