# calcul des écoulements en zone non saturée

## Ce que METIS fait en 2D (cas non saturé)

En 2D, METIS résout l'écoulement sur une **coupe** du terrain (plan $x$-$z$), discrétisée en éléments finis (triangles).

### 1) Domaine 2D et inconnue principale

- Le maillage 2D représente la géométrie de la coupe (surface, fond, berges).
- L'inconnue hydraulique principale est la **charge** `h` (ou potentiel hydraulique) aux nœuds du maillage.

### 2) Passage charge → succion → saturation → perméabilité relative

Pour chaque nœud/élément, METIS enchaîne :

$$
h \rightarrow \psi = h-z \rightarrow S(\psi) \rightarrow k_r(S)
$$

- $\psi$ : succion (ou pression capillaire selon convention interne),
- $S$ : saturation en eau,
- $k_r$ : perméabilité relative.

Ensuite, la conductivité hydraulique effective est mise à jour (typiquement $K_{eff}=K_{sat}\,k_r$).

### 3) Équation résolue en 2D

METIS résout en 2D la formulation en charge de Richards avec coefficients non linéaires dépendant de l'état hydrique.

L'inconnue primaire est la charge $h$ aux noeuds.

#### 3.1 Forme continue (niveau modèle)

- Permanent:

$$
\nabla \cdot \mathbf{q}(h) = Q
$$

- Transitoire:

$$
C(h)\,\partial_t h + \nabla \cdot \mathbf{q}(h) = Q
$$

avec:

$$
\mathbf{q}(h) = -\mathbf{K}_{sat}\,k_r(h)\,\nabla h
$$

et un terme capacitif effectif construit à partir de l'emmagasinement et de la dérivée de saturation:

$$
C(h) \sim S_s\,S(h) + \phi\,\frac{dS}{d\psi}(h)
$$

#### 3.2 Forme faible et discrétisation EF P1 triangulaire

Le coeur du calcul est construit dans `cteptns.f` puis assemblé dans `assco.f`/`asscog.f`.

Pour un élément triangulaire $e$, METIS calcule une matrice locale de diffusion `cdiff`:

$$
K^e_{ij} = -\frac{D_{xx}\beta_i\beta_j + D_{xy}(\beta_i\gamma_j+\gamma_i\beta_j) + D_{yy}\gamma_i\gamma_j}{12A_e}
\left(\sum_{p=1}^{3} e_p\,k_{r,p}\right)
$$

où $\beta_i,\gamma_i$ viennent de la géométrie P1 et $k_{r,p}$ est évalué aux points d'intégration.

En transitoire, METIS calcule une matrice locale de stockage `cstor` (lumpée):

$$
M^e_{ii} = f_i\left(S_s\,S_i + \phi\,\frac{dS}{d\psi}\Big|_i\right),
\qquad M^e_{ij}=0\ (i\neq j)
$$

La correction dite de Celia est introduite via des termes de second membre (`tsC`, `tsCI`, `ttres`) dans `cteptns.f`, puis propagée lors de l'assemblage.

#### 3.2.1 Détail EF: de la forme forte a la forme faible

Sur un domaine $\Omega \subset \mathbb{R}^2$, avec frontière de Dirichlet $\Gamma_D$ et de Neumann $\Gamma_N$, la forme faible cible (transitoire) est:

$$
\int_{\Omega} w\,C(h)\,\partial_t h\,d\Omega
+ \int_{\Omega} \nabla w \cdot \mathbf{K}(h)\nabla h\,d\Omega
= \int_{\Omega} w\,Q\,d\Omega
+ \int_{\Gamma_N} w\,\bar q\,d\Gamma
$$

avec $\mathbf{K}(h)=\mathbf{K}_{sat}k_r(h)$.

METIS utilise ensuite une approximation P1 nodale:

$$
h_h(x,z,t)=\sum_{j=1}^{N} N_j(x,z)\,H_j(t),
\qquad
w_h=N_i
$$

Le système semi-discret s'écrit:

$$
\mathbf{M}(\mathbf{H})\,\dot{\mathbf{H}} + \mathbf{K}(\mathbf{H})\,\mathbf{H}=\mathbf{F}
$$

où la non-linéarité est dans $\mathbf{M}$ et $\mathbf{K}$ via $S(\psi)$, $dS/d\psi$ et $k_r$.

#### 3.2.2 Géométrie triangle P1 et coefficients locaux

Pour un triangle $e$ d'aire $A_e$:

- les gradients des fonctions de forme sont constants sur l'élément,
- les coefficients géométriques $\beta_i,\gamma_i$ sont ceux calculés dans `cteptns.f`,
- les termes diffusifs locaux `cdiff` représentent l'intégrale
  $\int_{\Omega_e} \nabla N_i\cdot \mathbf{K}(h)\nabla N_j\,d\Omega$.

La pondération
$\sum_p e_p k_{r,p}$
correspond à une évaluation élémentaire des propriétés hydrauliques aux points d'intégration utilisés par METIS.

#### 3.2.3 Terme capacitif et mass-lumping

La matrice `cstor` est diagonalisée (mass-lumping), ce qui stabilise et simplifie la résolution non linéaire.

Interprétation:

- contribution "élastique"/emmagasinement: $S_s\,S$,
- contribution capillaire: $\phi\,dS/d\psi$,
- contributions hors-diagonale nulles dans la forme lumpée utilisée.

Ce point est important pour l'analyse numérique: la structure diagonale du stockage modifie les propriétés spectrales du système (et donc la convergence du solveur itératif) par rapport à une masse consistante complète.

#### 3.2.4 Assemblage global et structure algébrique

L'assemblage est fait en stockage Morse via `assco.f` / `asscog.f`:

- insertion des blocs élémentaires dans la matrice globale,
- gestion séparée ou globale diffusion+stockage selon le mode,
- construction du second membre incluant les termes de linéarisation.

Au niveau itératif, METIS résout successivement des systèmes linéaires préconditionnés (CG via `sGC`).

#### 3.2.5 Conditions aux limites dans la formulation EF

- Dirichlet (charge imposée): imposée algébriquement sur lignes/colonnes du système via `conlie`.
- Neumann (flux imposé): contribution au second membre frontière.

Dans ce cas AvAv, la majorité des CL hydrauliques sont de type potentiel imposé sur contours (listes de noeuds de contour).

#### 3.3 Schéma en temps et système algébrique

Le centrage temporel est piloté par `cmixt` (noté $\theta$):

- en permanent: implicite (pas de terme transitoire),
- en transitoire non saturé: typiquement $\theta=0.5$ (Crank-Nicolson) sauf forçage implicite.

Après assemblage global Morse, la forme linéarisée à l'itération $k$ peut se lire comme:

$$
\left[\theta K(h^{k}) - \frac{1}{\Delta t}M(h^{k}) + J(h^{k})\right] \delta h = -R(h^{k})
$$

ou, selon l'étape de l'algorithme, directement en variable absolue $h$ à la première itération.

- $J$ provient des jacobiens `cjacp` (partie "diffusive" non linéaire) et `cjact` (partie transitoire non linéaire) en Newton.
- $R$ est le résidu non linéaire traité dans les pilotes (`rpecns.f`, `calfecns.f`).

#### 3.4 Linéarisation et stratégie numérique

- Picard / Newton-Raphson (ou alternance) selon `moditer`.
- Contrôle de convergence par résidu max, norme L2, et variation de charge.
- Réduction adaptative du step non linéaire (`c_pstep`) en cas de non-décroissance du résidu.
- Résolution linéaire par gradient conjugué préconditionné (`sGC`).

#### 3.5 Parcours de code recommandé (niveau recherche)

Pour entrer rapidement dans le coeur numérique:

1. `ecperm.f` ou `calfecns.f`/`ectran.f`: orchestration du calcul et des conditions limites.
2. `rpecns.f` (permanent non saturé) et `calfecns.f` (transitoire non saturé): boucle non linéaire, critères, pas de temps.
3. `acmec.f`: passerelle assemblage global, choix 2D/3D, non saturé/saturé.
4. `cteptns.f`: définition explicite de `cdiff`, `cstor`, `cjacp`, `cjact`, `tsC`, `ttres`.
5. `assco.f` / `asscog.f`: insertion dans la matrice globale Morse et écriture finale du système linéaire.

Chemin minimal à suivre:

`ecperm.f` -> `rpecns.f` -> `acmec.f` -> `cteptns.f` -> `assco.f`/`asscog.f` -> `conlie` -> `sGC`

### 4) Conditions limites en 2D

Les conditions limites sont appliquées sur des **lignes de contour** du maillage 2D :

- `potentiel contour` / `potentiel boucle` : charge imposée,
- flux imposé : si option de flux activée.

Dans votre cas, `trace_riv_berge` est traité comme une **liste de nœuds de contour** (pas une liste de faces, sauf option explicite `arete`).

### 5) Ce que produit le calcul

À chaque pas de temps (ou à convergence en permanent), METIS calcule notamment :

- charge/potentiel,
- pression ou succion,
- saturation,
- vitesses/flux,
- grandeurs dérivées (ex. perméabilité relative/effective selon sorties demandées).

En résumé, en 2D METIS fait un calcul éléments finis sur une coupe verticale, avec une boucle non linéaire qui met à jour en continu les propriétés hydrauliques du milieu non saturé.



# Comprendre le fichier `.COMM` du cas AvAv

Ce document résume le comportement de METIS pour le fichier [AvAv_bck.COMM](AvAv_bck.COMM), à partir de l'analyse des sources METIS utilisées localement dans :

`/home/ariviere/Programmes/metis2017_devel_diabolix_060320_modif/`

L'objectif n'est pas de décrire tout METIS, mais de comprendre **comment le solveur lit ce `.COMM`**, **dans quel ordre il exécute les blocs**, et **ce que cela implique pour le cas AvAv**.

---

## 1. Architecture générale

Le fichier `.COMM` est lu en plusieurs niveaux :

1. **Parseur principal** dans `METIS.f`
   - lit les mots-clés globaux
   - initialise le job
   - ouvre le maillage
   - fixe le type de problème

2. **Routeur des blocs de calcul** dans `ecdis.f`
   - associe les noms de blocs du `.COMM` à des routines internes
   - enchaîne les calculs permanents et transitoires

3. **Solveurs spécialisés**
   - `ecperm.f` pour l'écoulement permanent
   - `ectran.f` pour les transitoires
   - `diperm` / `ditran` pour la thermique et la dispersion selon la variable demandée

---

## 2. Mots-clés globaux lus avant les calculs

Dans [AvAv_bck.COMM](AvAv_bck.COMM), les lignes suivantes définissent le cadre général :

- `nom_job = Case_unsat_heat`
- `maillage = AvAv_mesh.mail`
- `type_probleme = 2`
- `non_sature = 1`

### Signification

- **`nom_job`**
  - définit le préfixe utilisé pour certains fichiers de sortie.

- **`maillage`**
  - demande à METIS de lire le fichier de maillage `AvAv_mesh.mail`.
  - ce chargement est géré dans le parseur principal avant le lancement des solveurs.

- **`type_probleme = 2`**
  - signifie : **écoulement + thermique**.
  - valeurs observées dans les sources :
    - `0` : écoulement seul
    - `1` : écoulement + dispersion
    - `2` : écoulement + thermique
    - `3` : écoulement + dispersion + thermique

- **`non_sature = 1`**
  - active les traitements non saturés.
  - cela a un impact sur l'allocation des variables et sur les routines hydrauliques.

---

## 3. Routage des blocs de haut niveau

Dans les sources METIS, le fichier `ecdis.f` contient les alias textuels qui relient les noms du `.COMM` aux routines internes.

### Correspondance utile pour AvAv

- `ecoulement_permanent` → routine `ecperm`
- `thermique_permanente` → routine `diperm` avec la variable `Temperature`
- `thermique_transitoire` → routine `ditran` avec la variable `Temperature`

Autrement dit, dans le cas AvAv, METIS exécute les blocs comme une **suite d'étapes**.

---

## 4. Ordre d'exécution du cas AvAv

Avec le `.COMM` actuel, la logique est la suivante :

### Étape 1 — Écoulement permanent

Bloc : `ecoulement_permanent`

METIS lance un calcul hydraulique stationnaire avec les propriétés définies plus haut dans le fichier :

- perméabilité
- emmagasinement
- porosité
- loi non saturée
- capacité aquifère
- conductivité thermique

Dans ce bloc, la section `conditions_limites` lit surtout des **conditions de potentiel** (`potentiel contour`, `potentiel boucle`).

### Étape 2 — Thermique permanente

Bloc : `thermique_permanente`

METIS enchaîne avec un calcul thermique stationnaire.

Dans ce bloc, les `temperature contour` imposent des valeurs fixes sur :

- la surface
- les côtés verticaux
- le fond

### Étape 3 — Thermique transitoire

Bloc : `thermique_transitoire`

METIS lance ensuite un transitoire thermique avec :

- un `decoupage_temporel`
- une `initialisation_limites`
- des `conditions_limites`
- des `edition_resultats`

Le calcul thermique transitoire repose donc sur un état hydraulique déjà calculé précédemment.

---

## 5. Ce que font les conditions limites dans ce `.COMM`

### 5.1 Bloc hydraulique permanent

Dans [AvAv_bck.COMM](AvAv_bck.COMM), on trouve notamment :

- `potentiel contour`
- `potentiel boucle`

Cela signifie que l'écoulement permanent est piloté par :

- des potentiels imposés sur certains contours verticaux
- un potentiel imposé sur une boucle nommée `trace_riv_berge`

Le solveur transforme ces définitions géométriques en conditions appliquées sur le maillage.

**Note importante sur `trace_riv_berge`** : ce contour est défini dans le maillage (fichier `.mail`) comme une **liste de nœuds**, pas comme une liste de faces ou d'arêtes. METIS applique la condition en imposant la valeur du potentiel (ou de la température) directement sur ces nœuds. Cela signifie que dans le fichier `AvAv_mesh.mail`, la définition du contour `trace_riv_berge` doit être une séquence ordonnée de numéros de nœuds.

### 5.2 Bloc thermique permanent

Les `temperature contour` de ce bloc imposent des températures constantes sur les bords du domaine.

### 5.3 Bloc thermique transitoire

La ligne importante est :

- `temperature contour fichier = T_R.fic`

Cela signifie que la température sur ce contour n'est **pas constante** : elle est lue dans un **fichier temporel externe**.

METIS attend donc un fichier `T_R.fic` valide, lu comme une chronique de temps/valeur.

---

## 6. Point important sur les fichiers `.fic`

Actuellement, dans le dossier `METIS`, on génère :

- `Ch_R.fic`
- `Ch_RD.fic`
- `Ch_RG.fic`

avec `cree_fic3.f90`.

Mais dans [AvAv_bck.COMM](AvAv_bck.COMM), le bloc thermique transitoire utilise :

- `T_R.fic`

Donc, à ce stade :

- **`cree_fic3.f90` fonctionne**, mais
- **les fichiers générés ne correspondent pas encore au fichier attendu par la limite thermique transitoire**.

Il faudra donc soit :

1. produire un vrai `T_R.fic`, ou
2. modifier le `.COMM` si l'intention est d'utiliser d'autres fichiers de bord.

---

## 7. Comment la thermique est calculée avec la zone non saturée

Quand `non_sature = 1`, METIS ne calcule pas la thermique indépendamment de l'hydraulique.
Le calcul thermique utilise l'état hydrique du milieu non saturé.

### Chaîne physique utilisée par METIS

Dans les sources, on retrouve la logique suivante :

1. METIS calcule la **succion** ou pression capillaire à partir de la charge et de l'altitude :
   - $\psi = h - z$

2. À partir de cette succion, METIS calcule la **saturation** $S$ avec la loi non saturée définie dans le `.COMM` (`lois_nonsat`).

3. À partir de la saturation, METIS calcule la **perméabilité relative** $k_r$.

4. Cette perméabilité relative modifie les **vitesses d'écoulement de l'eau**.

5. Ces vitesses sont ensuite utilisées dans le calcul thermique, en particulier pour le **transport advectif de chaleur**.

On peut résumer cela par :

$$
h \rightarrow \psi \rightarrow S \rightarrow k_r \rightarrow \text{vitesse d'eau} \rightarrow \text{transport thermique}
$$

### Ce que cela implique

- si le milieu est **moins saturé**, la saturation diminue ;
- la perméabilité relative diminue aussi ;
- les vitesses d'eau deviennent plus faibles ;
- le transport advectif de chaleur diminue.

Autrement dit, dans la ZNS, la chaleur est moins transportée par l'eau quand le milieu s'assèche.

### Rôle des paramètres thermiques du `.COMM`

Dans ton fichier, les mots-clés principaux sont :

- `conductivite`
- `capacite_aquifere`
- `capacite_eau`

Leur rôle est le suivant :

- **`conductivite`**
  - contrôle la **conduction thermique** dans le milieu.

- **`capacite_aquifere`**
  - correspond à la **capacité calorifique volumique du milieu** utilisée dans le stockage thermique.

- **`capacite_eau`**
  - correspond à la **capacité calorifique massique de l'eau**.
  - elle intervient dans la part advective du transport thermique.

### Point subtil important

Dans les routines thermiques, la non-saturation n'agit pas seulement comme un coefficient appliqué directement à la température.
Elle agit surtout en modifiant l'état hydraulique local :

- saturation,
- perméabilité relative,
- vitesses,
- donc transport de chaleur.

Le stockage thermique du milieu, lui, repose principalement sur la valeur de `capacite_aquifere` fournie dans le `.COMM`.

### Application au cas AvAv

Dans `AvAv_bck.COMM`, cela signifie :

- `ecoulement_permanent` calcule d'abord l'état hydrique du système ;
- `thermique_permanente` calcule ensuite un état thermique stationnaire sur ce champ hydrique ;
- `thermique_transitoire` fait évoluer la température dans le temps en utilisant ce cadre, avec une condition thermique variable via `T_R.fic`.

Donc, même si le bloc s'appelle `thermique_transitoire`, la thermique reste bien influencée par la zone non saturée parce que METIS utilise les variables hydriques issues de cette formulation.

---

## 8. Lecture simplifiée de `AvAv_bck.COMM`

On peut résumer le cas ainsi :

- chargement d'un maillage METIS `AvAv_mesh.mail`
- problème couplé **hydraulique + thermique**
- prise en compte de la **zone non saturée**
- calcul d'un **écoulement permanent**
- calcul d'un **champ thermique permanent**
- calcul d'un **transitoire thermique** piloté par une chronique en fichier
- sorties sur certaines mailles et à certaines dates

---

## 9. Pourquoi `./metis2020 < AvAv_bck.COMM` peut échouer actuellement

Même si la structure générale du `.COMM` est cohérente, l'exécution peut échouer si l'un des éléments suivants manque ou est incohérent :

- fichier de maillage absent ou mal formé
- contour inclus (`trace_riv_berge`) absent
- fichier `T_R.fic` absent
- numéros de mailles de sortie incompatibles avec le maillage actuel
- paramètres ou géométries encore hérités d'une ancienne version du cas

Autrement dit : le `.COMM` est compréhensible, mais **pas encore totalement aligné avec les fichiers disponibles**.

---

## 10. Fichiers à vérifier en priorité pour remettre le cas à jour

Pour rendre le cas exécutable de bout en bout, vérifier d'abord :

- [AvAv_bck.COMM](AvAv_bck.COMM)
- [AvAv_mesh.mail](AvAv_mesh.mail)
- [AvAv_mesh.geom](AvAv_mesh.geom)
- **`trace_riv_berge` contour exists** : dans `AvAv_mesh.mail`, ce contour doit être défini comme une **liste ordonnée de nœuds** (pas de faces). METIS lit les trois premiers nœuds pour identifier la géométrie du contour, puis applique la condition aux nœuds de ce contour.
- `T_R.fic` pour la condition thermique transitoire
- les fichiers de sortie demandés sur des numéros de mailles spécifiques

---

## 11. Clarification : nœuds vs faces dans les contours METIS

### Question : comment METIS traite-t-il un contour comme `trace_riv_berge` ?

**Réponse** : par défaut, les contours du fichier `.COMM` sont appliqués comme des **conditions aux nœuds** (sauf option `arete`).

#### Cas standard : conditions aux nœuds

Quand vous écrivez dans le `.COMM` :

```
potentiel contour = trace_riv_berge  100.0
```

METIS :
1. demande au maillage tous les **nœuds** du contour `trace_riv_berge`
2. applique `potentiel = 100.0` **sur chacun de ces nœuds**

Cela s'opère dans la routine `lecvi_n.f` du solveur METIS (ligne 853-870 environ) :

```fortran
!     condition imposee aux noeuds
if (i_arete .eq. 0) then
    noeuds4: do  i=1,lliste
    ii = itrav2(i)       ! ii = numéro du nœud i
    nvi(novi)=ii
    vi(novi)=val         ! impose val au nœud ii
enddo noeuds4
```

#### Cas optionnel : conditions aux arêtes/faces

Si la condition était à imposer sur des **arêtes** (edges en 2D, faces en 3D), il faudrait une option explicite `arete` ou `option arete` dans le bloc `conditions_limites`, et aussi `option flux` ou `option EFMH`. Mais ce n'est **pas le cas par défaut** dans AvAv.

### Application : `trace_riv_berge` doit être une liste de nœuds

Donc, quand vous vérifiez que le contour `trace_riv_berge` existe dans le maillage `AvAv_mesh.mail`, assurez-vous qu'il est défini comme une **séquence de numéros de nœuds**, dans l'ordre géométrique autour de la limite rivière-berge.

---

## 12. Résumé très court

Le `.COMM` AvAv décrit un cas **hydraulique + thermique non saturé** où METIS :

1. lit le maillage,
2. calcule un **écoulement permanent**,
3. calcule un **régime thermique permanent**,
4. lance un **transitoire thermique** avec une température imposée via un fichier externe.

Le point de vigilance principal du cas actuel est que le `.COMM` attend `T_R.fic`, alors que les conversions préparées jusqu'ici produisent surtout des fichiers `Ch_*.fic`.

---

## 13. Bloc non saturé

Dans `AvAv_hydro.COMM`, le bloc suivant pilote tout le non saturé hydraulique :

```text
sous_domaine = 1

lois_nonsat 11
1
0 1 0.35 0 0
1
3.6 1.56 0 0 0
```
Un sous-domaine dans METIS, c’est une zone matérielle du maillage à laquelle on associe un jeu de paramètres.

Concrètement :

le maillage est découpé en éléments (triangles ici),
chaque élément appartient à un numéro de sous-domaine (1, 2, 3, …),
pour chaque sous-domaine, tu peux donner des propriétés différentes :
perméabilité,
porosité,
emmagasinement,
lois non saturées (ckr, cpc), etc.
Dans ton cas :

sous_domaine = 1 signifie “un seul matériau partout”.
donc les lignes avec 1 avant les coefficients (ckr/cpc) disent : “ces coefficients s’appliquent au sous-domaine 1”.
Si tu avais 2 matériaux, tu aurais des lignes pour 1 puis pour 2, avec des coefficients potentiellement différents.

### 13.1 Ce que METIS comprend

- `sous_domaine = 1` : le modèle a un seul matériau hydraulique.
- `lois_nonsat 11` : deux lois sont activées :
  - `loikr = 1` pour `k_r(S)` (perméabilité relative)
  - `loipc = 1` pour `S(\psi)` (capillarité)

### 13.2 Interprétation des coefficients

Pour le sous-domaine `1` :

- ligne `0 1 0.35 0 0` (`ckr`) :
  - `Smin = 0`
  - `Smax = 1`
  - `m = 0.35` (paramètre de la loi VG-Mualem)
  - les deux derniers coefficients ne sont pas utilisés ici

- ligne `3.6 1.56 0 0 0` (`cpc`) :
  - `\alpha = 3.6`
  - `n = 1.56`
  - `cpc(3)=0` : pas de succion max imposée
  - les deux derniers coefficients ne sont pas utilisés ici

### 13.3 Où c'est lu dans le code METIS

- lecture de `lois_nonsat` et décodage `loikr/loipc` : `lecpar.f`
- lecture des lignes `ckr/cpc` par sous-domaine : `lecckr.f`
- calcul de `k_r` : `fkrel.f`
- calcul de `S(\psi)` : `fhtosat.f`
- calcul inverse `\psi(S)` : `fsatpsi.f`
- noms officiels des lois : `edilnsa.f`

### 13.4 Formules utilisées pour `lois_nonsat 11`

`loikr=1` (VG-Mualem) :

$$
k_r = S_e^{1/2}\left(1-\left(1-S_e^{1/m}\right)^m\right)^2,
\quad
S_e = \frac{S-S_{min}}{S_{max}-S_{min}}
$$

`loipc=1` (VG capillaire) : relation saturation–succion avec `\alpha` et `n` (implémentée dans `fhtosat.f`/`fsatpsi.f`).

### 13.5 Autres lois disponibles dans METIS

Pour `loikr` (perméabilité relative `k_r`) :

- `1` : Van Genuchten (forme VG-Mualem)
- `2` : Brooks & Corey
- `4` : linéaire
- `5` : polynomiale (Huyakorn)
- `6` : exponentielle en `psi`
- `7` : milieu saturé (`k_r = 1`)
- `8` : présente dans le code (cas spécifique), plus rare

Pour `loipc` (courbe capillaire `S(\psi)`) :

- `1` : Van Genuchten
- `2` : Brooks & Corey
- `4` : linéaire
- `5` : polynomiale (Huyakorn)
- `7` : milieu saturé
- `8` : Van Genuchten modifiée (Ippisch/Vogel/Bastian)

Les indices `3`, `9`, `10` ne sont pas des lois standard utilisées ici.

Si vous changez de loi, vérifiez toujours :

1. le nombre de coefficients attendu,
2. l'ordre des coefficients,
3. le numéro de sous-domaine avant chaque ligne de coefficients.
