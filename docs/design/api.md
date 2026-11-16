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

Les **lignes** de l'inconnue sont toujours les points de reconstruction : `Groups`, `L` et
`Lpq` ont donc le même sens quel que soit le domaine (`L` agit sur la dimension spatiale).

## Termes de régularisation (= lois a priori)

* `Lq(q; L = I)` — `‖L u‖_q^q`, `0 < q ≤ 2` ;
* `Lpq(p, q; L = I)` — `Σᵢ (Σⱼ |(L U)ᵢⱼ|^p)^(q/p)` : `q` sur l'espace, `p` sur la fréquence/le temps ;
* `Groups(rows₁ => Lq(0.5), rows₂ => Lq(2.0; L = D₂))` — régularisation par régions ;
  chaque groupe a sa propre matrice `L_g` (`size(L_g, 2) == length(rows_g)`), l'identité par défaut.

Les `L_g` sont assemblés une fois en un opérateur global creux ; si tous sont l'identité,
le chemin diagonal rapide est conservé.

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

## Modèle bayésien : cas réel et cas complexe

| | réel | complexe (circulaire) |
|:--|:--|:--|
| vraisemblance | `(τn/2π)^{N/2} exp(−τn‖r‖²/2)` | `(τn/π)^N exp(−τn‖r‖²)` |
| a priori (par composante de `L_g u_g`) | `q τs^{1/q} / (2Γ(1/q)) exp(−τs|v|^q)` | `q τs^{2/q} / (2πΓ(2/q)) exp(−τs|v|^q)` |
| `τn` | `(N/2 + a − 1) / (‖r‖²/2 + b)` | `(N + a − 1) / (‖r‖² + b)` |
| `τs_g` | `(M_g/q_g + a − 1) / (Σ|v|^q + b)` | `(2M_g/q_g + a − 1) / (Σ|v|^q + b)` |
| `λ_g` (étape en `u`) | `2 τs_g / τn` | `τs_g / τn` |

Initialisation : `τn = N / yᴴ(I + H(Lᴴλ₀WL)⁻¹Hᴴ)⁻¹y` (identité de Woodbury, valable même
si `Lᴴ W L` est singulière), `τs = λ₀ c τn`.

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

* `solve(...; R)` : covariance du bruit **connue et diagonale** (`I` par défaut, `σ²I`,
  `Diagonal`, vecteur des variances ou matrice des variances au format de `y`). L'attache aux
  données devient `‖y − Hu‖²_{R⁻¹}` ; elle est traitée une fois pour toutes par blanchiment
  (`H ← R^{-1/2}H`, `y ← R^{-1/2}y`). Dans le modèle bayésien, la covariance du bruit est `R/τn` ; par défaut `τn = 1` si `R` est fourni (covariance exacte), `τn` optimisé sinon.
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

## Mots-clés communs de `solve`

`R` (covariance du bruit, non fournie par défaut ⇔ `I`), `covariance = false`, `tol = 1e-4`, `maxiter = 100`,
`init` (paramètre de la solution de Tikhonov initiale : `Heuristic()`, `LCurve()`, `GCV()`,
`BayesEstimator()` ou un nombre).

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
│   ├── linsolve.jl         # caches, résolutions de Tikhonov pondérées, covariance a posteriori
│   ├── spectrum.jl         # SVD / GSVD du problème sous forme standard
│   ├── regparam.jl         # L-curve, GCV, estimateur bayésien, heuristique
│   ├── norms.jl            # valeurs des normes, paramètre d'amortissement ε, poids IR
│   ├── bayes.jl            # modèle bayésien réel / complexe, lois conditionnelles
│   └── uncertainty.jl      # posterior_std, posterior_samples, credible_interval
└── solvers/                # algorithmes de résolution
    ├── ir.jl               # IR : formulations additive et multiplicative, `solve`
    └── map.jl              # MAP : formulation bayésienne hiérarchique
```

* `types/` ne contient que des définitions de types, leurs constructeurs et leur validation ;
  c'est là que se lit l'API.
* `utils/` regroupe les briques numériques partagées. C'est le seul endroit qui dépend de la
  structure de `H` (dense / bloc-diagonal) et de `L` (identité / matrice) : un nouvel opérateur
  se branche en ajoutant des méthodes à `weighted_solve`, `spectrum` et `posterior_covariance`.
* `solvers/` contient une méthode de `solve` par couple (formulation, algorithme) : un nouvel
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
