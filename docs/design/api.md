# InverseProblems.jl — conception de l'API

## Principe

Toute l'API publique tient en un verbe et quatre familles de types :

```julia
prob = InverseProblem(H, y)                          # le problème
reg  = Multiplicative(Groups(r₁ => Lq(0.5),          # la formulation + les termes
                             r₂ => Lq(2.0; L = D)))  # (un opérateur de lissage par groupe)
sol  = solve(prob, reg)                              # algorithme par défaut
sol.u, sol.λ, sol.hyper                              # résultat
```

`solve(prob, formulation, algorithme)` est entièrement résolu par multiple dispatch :

| Axe de dispatch          | Types                                                          | Ce qui change                                        |
|:-------------------------|:---------------------------------------------------------------|:-----------------------------------------------------|
| Opérateur direct         | `DenseOperator`, `BlockDiagonalOperator`                       | résolution linéaire (bloc par fréquence, multithreadée) |
| Opérateur de lissage     | `IdentitySmoothing`, `MatrixSmoothing`                         | résolution diagonale rapide vs `Hᴴ H + Lᴴ D L` ; SVD vs GSVD |
| Terme                    | `Lq`, `Lpq`                                                    | poids IR et pénalité (`reweight!`)                    |
| Formulation              | `Additive{λ::Real / Vector / ParameterChoice}`, `Multiplicative{Sum/Product}`, `Bayesian` | mise à jour du paramètre de régularisation |
| Règle de choix de λ      | `LCurve`, `GCV`, `BayesEstimator`, `Heuristic`                 | critère sur le spectre (généralisé)                    |
| Corps des données        | `RealField`, `ComplexField`                                    | lois réelles / complexes du modèle bayésien            |

Le branchement à l'exécution n'est utilisé que pour la validation des arguments, le choix
`nm ≥ nr` / `nm < nr` dans la résolution diagonale (équations normales vs identité de
*push-through*), l'assemblage de `L` (tout identité ou non) et les tests d'arrêt.

## Problème

| `H`                                  | `y`             | inconnue `u`   | cas                |
|:-------------------------------------|:----------------|:---------------|:-------------------|
| matrice `nm × nr`                    | vecteur         | vecteur `nr`   | une fréquence      |
| vecteur de matrices / tableau 3D     | matrice `nm×nf` | matrice `nr×nf`| espace–fréquence   |
| matrice `nm·nt × nr·nt`              | matrice `nm×nt` | matrice `nr×nt`| espace–temps       |
| matrice **creuse** (`SparseOperator`) | vecteur ou matrice | idem      | problème assemblé (p. ex. espace–fréquence vectorisé, fréquences couplées) |

Les **lignes** de l'inconnue sont toujours les points de reconstruction : `Groups`, `L` et
`Lpq` ont donc le même sens quel que soit le domaine (`L` agit sur la dimension spatiale).

## Termes de régularisation (= lois a priori)

* `Lq(q; L = I)` — `‖L u‖_q^q`, `0 < q ≤ 2` ;
* `Lpq(p, q; L = I, Lf = nothing, samples = nothing)` — `‖L U Lfᵀ‖_{p,q}^q =
  Σᵢ (Σⱼ |(L U Lfᵀ)ᵢⱼ|^p)^(q/p)` : `q` sur l'espace, `p` sur la fréquence/le temps (MSSP 2018,
  éq. 10–15). `L` agit sur l'espace, `Lf` le long des fréquences/instants (p. ex.
  `Lf = diff_operator(nf, 1; h = Δf)`, éq. 11, pour la continuité des spectres) ; l'opérateur
  global est `kron(Lf, L)` (éq. 15). `samples` sert aux problèmes déjà vectorisés (voir plus
  bas) et ne se combine pas avec `Lf` ;
* `Groups(rows₁ => Lq(0.5), rows₂ => Lq(2.0; L = D₂))` — régularisation par régions ;
  chaque groupe a sa propre matrice `L_g` (`size(L_g, 2) == length(rows_g)`), l'identité par défaut.

Les `L_g` sont assemblés une fois en un opérateur global creux ; si tous sont l'identité,
le chemin diagonal rapide est conservé.

### Opérateurs de lissage (différences finies)

* `fd_weights(h, npts, order)` — poids de différences finies sur une grille uniforme (port
  de `ufdwt`) : une ligne par schéma, décentré en avant (ligne 1), centré (ligne du milieu)
  ou décentré en arrière (dernière ligne) ;
* `diff_operator(n, order; npts, h)` — opérateur 1D (port de `Dx_1D`) : schéma minimal
  (`npts = order + 1`, matrice `(n − order) × n` sans lignes de bord : différences avant
  pour l'ordre 1, différences secondes des points intérieurs pour l'ordre 2) ou centré
  (`npts` impair `> order + 1`, matrice `n × n`, bords décentrés). Par défaut, `npts` est le
  plus petit nombre impair de points (3 pour les ordres 1 et 2) ;
* `diff_operator(coords, dir; order, npts, atol)` — opérateur le long de la direction `dir`
  d'un nuage de points structuré (maillage 2D/3D, remplace `Dx_2D`/`Dy_2D`) : les points
  sont regroupés en lignes de maillage le long de `dir` et l'opérateur 1D est appliqué sur
  chaque ligne avec son propre pas (uniforme). Les points peuvent être dans un ordre
  quelconque ; les colonnes de `L` suivent l'ordre de `coords`.

```julia
Groups(inner => Lq(0.5),
       edges => Lq(2.0; L = diff_operator(coords[edges, :], 1)))
Lq(2.0; L = [diff_operator(coords, 1); diff_operator(coords, 2)])   # gradient
```

* `space_sample_operator(nr, ns; space = I, samples = I)` — opérateur d'un problème
  **vectorisé** espace–fréquence ou espace–temps (`u = vec(U)`, `U` : `nr × ns`, colonnes =
  fréquences ou instants) : `L = kron(samples, space)`, soit
  `L vec(U) = vec(space · U · samplesᵀ)`. `samples = Ls` seul : lissage le long des
  fréquences/instants en chaque point (le `L` de `lissage_freq`) ; `space = Lx` seul :
  lissage spatial à chaque fréquence/instant ; les deux : produit (dérivée croisée). Pour
  une somme de pénalités, empiler les opérateurs (`[L₁; L₂]`) ;
* `sample_rows(rows, nr, ns)` — indices des points `rows` à toutes les fréquences/instants
  dans `vec(U)` (port de `define_group_mixed`), pour définir les groupes d'un problème
  vectorisé ; l'opérateur du groupe est alors `space_sample_operator(length(rows), ns; …)`.

Dans un problème espace–fréquence ou espace–temps (données matricielles), l'opérateur le
long des fréquences ou du temps se donne directement par `Lpq(…; Lf)` : le problème est
alors assemblé automatiquement (voir « Stratégie de résolution »). `space_sample_operator`
et `sample_rows` servent aux problèmes que l'on vectorise soi-même.

**Problèmes vectorisés : utilisation des normes mixtes.** Dans un problème vectorisé, `L u`
est un vecteur : `Lpq` doit connaître sa disposition matricielle (points × échantillons).
Le mot-clé `samples` donne le nombre de colonnes de `L u` ; le vecteur du groupe est replié
colonne par colonne (points les plus rapides, comme `vec(U)`, `space_sample_operator` et
`sample_rows`) avant le calcul de la norme mixte :

| `L` du groupe | `samples` |
|:--|:--|
| `I` | `ns` |
| `space_sample_operator(n, ns; space = Lx)` | `ns` |
| `space_sample_operator(n, ns; samples = Ls)` (± `space`) | `size(Ls, 1)` |

```julia
Groups(sample_rows(inner, nr, nf) => Lpq(2, 0.5; samples = nf),
       sample_rows(edges, nr, nf) => Lpq(2, 1; L = space_sample_operator(length(edges), nf; samples = Ls),
                                         samples = size(Ls, 1)))
```

`Lpq(p, q; samples = nf)` sur le problème vectorisé est exactement le modèle du problème
espace–fréquence couplé avec `Lpq(p, q)`. Une incohérence (`length(L u)` non multiple de
`samples`) lève une `DimensionMismatch`, et une norme mixte sans `samples` sur une inconnue
à une seule colonne lève une `ArgumentError` (elle se réduirait sinon silencieusement à `Lq`).

Les opérateurs sont des matrices creuses et sont divisés par `h^order` (contrairement aux
opérateurs 2D du code Python, non mis à l'échelle). `selection_matrix` et
`define_group_mixed` n'ont plus d'équivalent : `Groups` et la disposition « lignes = espace »
les remplacent.

## Formulations

| Appel                                               | Méthode (article)                                   |
|:----------------------------------------------------|:----------------------------------------------------|
| `Additive(prior; λ = 1e-3)`                          | Tikhonov généralisé (JSV 2014)                      |
| `Additive(prior; λ = LCurve())`, `GCV()`, `BayesEstimator()` | λ choisi à chaque itération                  |
| `Additive(Groups(...); λ = [λ₁, λ₂])`                | additive multi-paramètres                            |
| `Multiplicative(prior)`                              | OMR (MSSP 2017, 2018 ; JSV 2020 ; MSSP 2020)         |
| `Multiplicative(Groups(...); multiparameter = true)` | MPMR (JSV 2020)                                      |
| `Bayesian(Groups(r => Lq(q₀), ...))`                 | régularisation bayésienne optimale (MSSP 2019)       |

## Choix du paramètre

Toutes les règles travaillent sur la description spectrale de
`min ‖y − Hu‖² + λ‖W^{1/2} L u‖²` : SVD de `H W^{-1/2}` si `L = I`, GSVD de `(H, W^{1/2} L)`
sinon (les composantes du noyau de `L` ne sont pas régularisées).

* `LCurve()` — courbure maximale de la courbe en L ;
* `GCV()` — validation croisée généralisée ;
* `BayesEstimator()` — minimum de `J_MAP(ℓ) = Σ log(γᵢ² + ℓ²) + (M − 2) log(Σ |βᵢ|²/(γᵢ² + ℓ²) / M)` (`reg_param_bayes`) ;
* `Heuristic()` — initialisation seulement (MSSP 2017).
* `LambdaMax()` — initialisation seulement : `λ = ‖2 Hᴴ y‖∞` (plus petite valeur pour
  laquelle la solution ℓ1 est nulle) ; initialisation de `lpq_reg` et de la RVR. Par exemple,
  `lpq_reg(A, b, p, q, R)` correspond à
  `solve(InverseProblem(A, b), Multiplicative(Lpq(p, q)); R = diag(R), init = LambdaMax(), maxiter = 50)`.

## Modèle bayésien : cas réel et cas complexe

| | réel | complexe (circulaire) |
|:--|:--|:--|
| vraisemblance | `(τn/2π)^{N/2} exp(−τn‖r‖²/2)` | `(τn/π)^N exp(−τn‖r‖²)` |
| a priori (par composante de `L_g u_g`) | `q^{1−1/q} τs^{1/q} / (2Γ(1/q)) exp(−(τs/q)|v|^q)` | `q τs^{2/q} / (2πΓ(2/q)) exp(−τs|v|^q)` |
| `τn` | `(N/2 + a − 1) / (‖r‖²/2 + b)` | `(N + a − 1) / (‖r‖² + b)` |
| `τs_g` | `(M_g/q_g + a − 1) / (Σ|v|^q/q_g + b)` | `(2M_g/q_g + a − 1) / (Σ|v|^q + b)` |
| objectif en `q_g` | `M logΓ(1/q) − (M/q) log(τs/q) − (a+M−1) log q + bq + (τs/q)Σ|v|^q` | `M logΓ(2/q) − (2M/q) log τs − (a+M−1) log q + bq + τsΣ|v|^q` |
| MAP en `u` (`λ_g = τs_g/τn`) | `½‖y − Hu‖²_R + Σ_g (λ_g/q_g)‖L_g u_g‖_q^q` | `‖y − Hu‖²_R + Σ_g λ_g‖L_g u_g‖_q^q` |

Soit, avec les facteurs `exp(−c τn‖r‖²)` et `exp(−φ_g τs_g‖L_g u_g‖_q^q)` : `c = 1/2`,
`φ_g = 1/q_g` en réel, `c = φ_g = 1` en complexe. `sol.λ` contient `λ_g = τs_g/τn` ;
l'algorithme IR résout le problème équivalent `‖y − Hu‖²_R + Σ_g λu_g‖L_g u_g‖_q^q` avec
`λu_g = φ_g τs_g / (c τn)` (`= 2λ_g/q_g` en réel, `λ_g` en complexe). Pour `q = 2`,
les deux lois a priori sont des lois normales (réelle ou complexe) de précision `τs`.

Initialisation : `τn = N / yᴴ(I + H(Lᴴλ0WL)⁻¹Hᴴ)⁻¹y` (identité de Woodbury, valable même
si `Lᴴ W L` est singulière), `τs_g = λ0 c τn / φ_g`.

### Hyperparamètres optimisés ou fixés

Chaque hyperparamètre est soit optimisé (`Estimated()`), soit fixé ; la mise à jour est
choisie par dispatch sur cette spécification :

```julia
Bayesian(prior)                                         # τs, q optimisés ; τn = 1 si R est fourni, optimisé sinon
Bayesian(prior; τn = Estimated())                       # R connue à un facteur près
Bayesian(Groups(...); τs = [τ₁, τ₂], q = Fixed())       # τs fixés par groupe, q de l'a priori
Bayesian(Lq(2.0); τn = 1.0, τs = 0.5, q = Fixed())      # tout fixé : posterior gaussien exact
```

## Covariance du bruit et incertitude a posteriori

* `solve(...; R)` : covariance du bruit **connue et diagonale** — non fournie (`I`), `σ²I`,
  variances de toutes les données (vecteur de longueur `length(y)`, matrice de la taille de
  `y`, ou matrice diagonale `length(y) × length(y)`), ou, pour des données matricielles
  (`nm × ns`), variances des `nm` points de mesure communes à tous les échantillons (vecteur
  de longueur `nm`, `Diagonal` ou matrice diagonale `nm × nm`). Une matrice de la taille de
  `y` est toujours lue comme les variances de toutes les données. L'attache aux
  données devient `‖y − Hu‖²_{R⁻¹}` ; elle est traitée une fois pour toutes par blanchiment
  (`H ← R^{-1/2}H`, `y ← R^{-1/2}y`). Dans le modèle bayésien, la covariance du bruit est
  `R/τn` ; par défaut `τn = 1` si `R` est fourni (covariance exacte), `τn` optimisé sinon.
* `solve(...; covariance = true)` : covariance a posteriori de l'approximation gaussienne au
  point fixe, `Σ = (τn (H̃ᴴH̃ + LᴴDL))⁻¹` avec `D = λW` les poids IR (exacte pour `q = 2`) ;
  `τn = 1` pour les formulations additive et multiplicative (`R` supposée exacte). Pour un
  opérateur bloc-diagonal, `Σ` est un vecteur de matrices `nr × nr` (une par fréquence).
* `posterior_std(sol)` : écarts-types marginaux (diagonale de `Σ`) ;
* `posterior_samples(sol, n; rng)` : échantillons `u = μ + S z`, `Σ = S Sᴴ` (Cholesky),
  `z ~ randn(eltype(μ))` (loi normale complexe circulaire en complexe) ; une matrice `Σ` par
  fréquence pour un opérateur bloc-diagonal ;
* `credible_interval(sol; level = 0.95, nsamples = 10_000, rng)` : quantiles empiriques
  `(1 − level)/2` et `(1 + level)/2` des échantillons, composante par composante (parties
  réelle et imaginaire séparément en complexe).

## Algorithmes

* `IR()` — condition d'optimalité du premier ordre (JSV 2020, MSSP 2020), seul algorithme
  des formulations additive et multiplicative ;
* `MAP()` — maximum a posteriori du modèle hiérarchique (MSSP 2019), défaut de `Bayesian`.
* `ElementWise()` — bayésien élément par élément (*Relevance Vector Regression*) pour
  `Bayesian` : une précision a priori gaussienne `τ_i` par composante de `L u`, estimée
  conjointement avec `u` :

      τ_i = (k/2 + a − 1) / (φ |v_i|² + b)      (= 1/(|v_i|² + 2b) en réel, a/(|v_i|² + b) en complexe)
      u   = (Hᴴ H + Lᴴ Λ L)⁻¹ Hᴴ y,   Λ = φ τ / (c τn) = τ / τn

  Initialisation : Tikhonov avec `λ = ‖2 Hᴴ y‖∞` (par fréquence pour l'opérateur
  bloc-diagonal). Seuls les opérateurs de lissage de l'a priori sont utilisés : les termes
  doivent être `Lq(2.0)` et `τs = Estimated()` ; `τn` est fixé ou estimé comme pour `MAP()`.
  `sol.λ` contient `Λ` (format de `L u`), `sol.hyper.τs` les précisions `τ_i`.

  ```julia
  sol = solve(prob, Bayesian(Lq(2.0)), ElementWise(); R = σ², covariance = true)
  ```

## Stratégie de résolution des problèmes espace–fréquence

Pour un problème espace–fréquence (opérateur bloc-diagonal), la stratégie est choisie par
dispatch sur les termes de régularisation (`sample_coupling`), comme dans les articles :

| termes de régularisation | stratégie |
|:--|:--|
| `Lq` uniquement (éventuellement en `Groups`) | **fréquence par fréquence** : chaque fréquence est un problème indépendant, avec ses propres paramètres de régularisation, hyperparamètres et critère d'arrêt, résolu en parallèle |
| au moins un `Lpq` (sans `Lf`) | **problème espace–fréquence couplé, résolu bloc par bloc** : la norme mixte couple les fréquences par les poids seulement, les systèmes linéaires restent bloc-diagonaux (MSSP 2018 avec `L = I`, MSSP 2020) |
| au moins un `Lpq` avec `Lf` | **problème assemblé** : l'opérateur `kron(Lf, L)` couple les fréquences dans les systèmes linéaires (bloc-tridiagonaux pour un `Lf` à 3 points) ; le problème est vectorisé automatiquement (`H` bloc-diagonal creux, groupes `sample_rows`, opérateurs `space_sample_operator`, `samples = size(Lf, 1)`) et résolu par Cholesky creux (MSSP 2018, éq. 10–18) |

Le problème assemblé est aussi utilisé pour un problème espace–temps (matrice de convolution
dense) dont un terme a un `Lf` (lissage temporel, MSSP 2018 §7). Le résultat a la disposition
habituelle (`u` : `nr × ns`) ; la covariance a posteriori est alors une matrice
`(nr·ns) × (nr·ns)` (les fréquences sont couplées). `Lf` sur un problème à une seule
fréquence, ou sur un problème déjà vectorisé, lève une `ArgumentError`.

```julia
Lf  = diff_operator(nf, 1; h = Δf)                                    # éq. (11)
sol = solve(InverseProblem(Hs, Y), Multiplicative(Lpq(2, 0.5; Lf)))   # éq. (18) avec L = Lf ⊗ I
```

Il n'y a pas de couplage possible avec `Lq` seul pour un problème espace–fréquence. Pour
coupler explicitement les fréquences autrement (un seul paramètre adaptatif pour `Lq`…), on
assemble soi-même le problème vectorisé avec un `H` **creux** :

```julia
prob = InverseProblem(blockdiag(sparse.(Hs)...), vec(Y))    # (ny·nf) × (nu·nf), u de taille nu·nf
sol  = solve(prob, Multiplicative(Lq(0.5)))                 # un seul α pour toutes les fréquences
```

Un `H` creux donne un `SparseOperator`, résolu en bloc par factorisations de Cholesky creuses
(CHOLMOD, analyse symbolique réutilisée d'une itération à l'autre) : les systèmes
`Hᴴ H + Lᴴ D L` restent creux (blocs + bande). Les groupes et `L` portent alors sur les
`nu·nf` composantes de l'inconnue vectorisée. Les règles fondées sur la SVD (`LCurve`,
`GCV`, `BayesEstimator`) ne sont pas disponibles (erreur explicite) ; `Heuristic()` utilise
des estimations peu coûteuses des valeurs propres extrêmes de `A = (H L⁻¹)ᴴ(H L⁻¹)` pour un
`L` carré inversible (MSSP 2018, éq. 26–27 ; `A` et `A⁻¹ = L (HᴴH)⁻¹ Lᴴ` appliqués comme
opérateurs, sans être assemblés), de `A = Hᴴ H` sinon : `s_max = ‖A‖₂` et `s_min = 1/‖A⁻¹‖₂`,
estimés par la méthode de la puissance (`normest`, `condest` en norme 2, comme
`init_alpha`). Cette heuristique sert
d'initialisation par défaut du MAP. Les problèmes à une fréquence et espace–temps sont toujours
résolus en bloc. En interne, `solve` choisit la stratégie (`Joint` / `PerSample` / `Assembled`) puis
appelle `solve_joint` (une méthode par couple formulation × algorithme) sur le problème
complet ou sur chaque fréquence.

Résultat d'une résolution fréquence par fréquence : `u` est `nr × nf` ; les grandeurs
propres à chaque fréquence sont empilées en colonnes (`λ` : vecteur de longueur `nf`, ou
matrice `G × nf` / `nr × nf` ; `hyper.τn`, `hyper.τs`, `hyper.q` de même) ; `hyper` contient
aussi `iterations`, `converged` et `history` par fréquence ; la covariance est un vecteur de
matrices `nr × nr`.

```julia
sol = solve(InverseProblem(Hs, Y), Multiplicative(Lq(0.5)); R = σ²)   # un α par fréquence
sol = solve(InverseProblem(Hs, Y), Multiplicative(Lpq(2, 0.5)))       # problème couplé
```

## Mots-clés communs de `solve`

`R` (covariance du bruit, non fournie par défaut ⇔ `I`), `covariance = false`, `tol = 1e-4`, `maxiter = 100`,
`init` (paramètre de la solution de Tikhonov initiale : `Heuristic()`, `LCurve()`, `GCV()`,
`BayesEstimator()` ou un nombre).

## Résolution linéaire : régime du problème et identité de Woodbury

Tous les algorithmes se ramènent à des résolutions de `(Hᴴ H + Lᴴ D L) u = Hᴴ y` (`D`
diagonale positive, `H` blanchie) et, pour la covariance a posteriori, à l'inverse de la même
matrice. Le régime du problème est un paramètre de type du cache (`Overdetermined` si
`nm ≥ nr`, `Underdetermined` si `nm < nr`) ; la méthode est choisie par dispatch :

| régime                              | `L`                  | méthode                                  | système factorisé |
|:------------------------------------|:---------------------|:-----------------------------------------|:------------------|
| `nm ≥ nr` (sur- ou parfaitement déterminé) | quelconque    | équations normales                       | `nr × nr`         |
| `nm < nr` (sous-déterminé)          | `I`                  | identité de Woodbury                     | `nm × nm`         |
| `nm < nr`                           | carrée inversible    | Woodbury sur la forme standard `H̃ = H L⁻¹` | `nm × nm`         |
| `nm < nr`                           | rectangulaire / singulière | équations normales (`Lᴴ D L` singulière : Woodbury inapplicable) | `nr × nr` |

avec

    (Hᴴ H + D)⁻¹ Hᴴ y = D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ y
    (Hᴴ H + D)⁻¹      = D⁻¹ − D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ H D⁻¹
    (Hᴴ H + Lᴴ D L)⁻¹ = L⁻¹ (H̃ᴴ H̃ + D)⁻¹ L⁻ᴴ                       (L carrée inversible)

Dans le cas parfaitement déterminé (`nm = nr`), les deux systèmes ont la même taille ; les
équations normales sont utilisées (`Hᴴ H` est alors calculée une seule fois). Dans le cas
sous-déterminé, `Hᴴ H` n'est ni stockée ni calculée, sauf pour un `L` rectangulaire.
L'inversibilité d'un `L` carré est testée une fois à l'assemblage (factorisation LU creuse,
`InvertibleSmoothing` / `SingularSmoothing`). Pour l'opérateur bloc-diagonal, chaque fréquence
suit ce schéma indépendamment ; pour l'espace–temps, `L` agit par blocs (`I ⊗ L`).
L'initialisation de `τn` utilise également Woodbury : `yᴴ(I + H M⁻¹ Hᴴ)⁻¹y = yᴴ(y − H û)`.

**Allocations.** Un système linéaire (`linear_system` : `DenseSystem`, ou `BlockSystem` avec un système par
fréquence pour le calcul parallèle) est construit une fois par résolution : il précalcule ce
qui ne dépend pas des poids (`Hᴴ H`, transposée de l'opérateur de lissage global, forme
standard `H L⁻¹`) et préalloue la matrice du système, `K D⁻¹`, le second membre, la solution
et le résidu. Les boucles (IR, MAP, `ElementWise`) utilisent les versions en place
`weighted_solve!`, `residual_norm2!`, `apply!`, échangent les tampons `U`/`Unew` et calculent
le critère d'arrêt sans tableau temporaire : une itération n'alloue plus de tableau de la
taille du problème. Les versions allouantes (`weighted_solve`, `residual_norm2`) restent
disponibles pour les calculs ponctuels.

## Organisation du code

```
src/
├── InverseProblems.jl      # module : include et export regroupés fichier par fichier
├── types/                  # définitions de types (API publique)
│   ├── problem.jl          # opérateurs directs (dense, bloc-diagonal), InverseProblem
│   ├── priors.jl           # termes de régularisation : Lq, Lpq, Groups
│   ├── choices.jl          # règles de choix du paramètre : LCurve, GCV, BayesEstimator, Heuristic
│   ├── formulations.jl     # Additive, Multiplicative, Bayesian, Estimated, Fixed
│   ├── algorithms.jl       # IR, MAP, default_algorithm
│   └── solution.jl         # InverseSolution
├── utils/                  # fonctions internes
│   ├── smoothing.jl        # assemblage des L_g en opérateur de lissage global
│   ├── noise.jl            # covariance du bruit R, blanchiment des données
│   ├── linsolve.jl         # caches (régime), résolutions pondérées (équations normales / Woodbury), covariance
│   ├── sparsesolve.jl      # opérateurs creux : Cholesky creux (CHOLMOD), heuristique
│   ├── minimize.jl         # minimisation 1D : Brent, recherche sur grille raffinée
│   ├── spectrum.jl         # SVD / GSVD du problème sous forme standard
│   ├── regparam.jl         # L-curve, GCV, estimateur bayésien, heuristique
│   ├── fdiff.jl            # opérateurs de différences finies : fd_weights, diff_operator
│   ├── norms.jl            # valeurs des normes, paramètre d'amortissement ε, poids IR
│   ├── bayes.jl            # modèle bayésien réel / complexe, lois conditionnelles
│   └── uncertainty.jl      # posterior_std, posterior_samples, credible_interval
└── solvers/                # algorithmes de résolution
    ├── ir.jl               # IR : formulations additive et multiplicative
    ├── map.jl              # MAP : formulation bayésienne hiérarchique
    ├── elementwise.jl      # bayésien élément par élément (RVR)
    └── strategy.jl         # `solve` : stratégie (couplée / fréquence par fréquence)
```

* `types/` ne contient que des définitions de types, leurs constructeurs et leur validation ;
  c'est là que se lit l'API.
* `utils/` regroupe les briques numériques partagées. C'est le seul endroit qui dépend de la
  structure de `H` (dense / bloc-diagonal) et de `L` (identité / matrice) : un nouvel opérateur
  se branche en ajoutant des méthodes à `weighted_solve`, `spectrum` et `posterior_covariance`.
* `solvers/` contient une méthode de `solve_joint` par couple (formulation, algorithme), `strategy.jl` le `solve` public : un nouvel
  algorithme (p. ex. `Gibbs` pour l'inférence bayésienne complète) se branche en ajoutant un
  type dans `types/algorithms.jl` et un fichier dans `solvers/`.
* Dans `InverseProblems.jl`, chaque `include` est suivi des `export` qu'il fournit ; l'ordre
  des `include` respecte les dépendances (types → utils → solvers).

## Retiré

* GMR, régularisation itérée, algorithmes IRLS ;
* du code Python : recalcul adaptatif de ε (`eps_adapt`), initialisation « 1 % de λ_max »,
  norme `ℓp` sur la fidélité aux données.

## Prochaines étapes

1. Inférence bayésienne complète (Gibbs + HMC, MSSP 2018) : `solve(prob, Bayesian(...), Gibbs(nsamples, burnin))`.
2. Construction de l'opérateur espace–temps par la méthode α-généralisée (MSSP 2019).
3. Documentation (Documenter.jl) et exemples reproduisant les cas tests des articles.
