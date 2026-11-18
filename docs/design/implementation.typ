#set document(title: "SourceEstimation.jl — points d'implémentation")
#set page(paper: "a4", margin: 2.2cm, numbering: "1")
#set text(lang: "fr", size: 10.5pt, font: "Lato")
#set heading(numbering: "1.1")
#set par(justify: true)
#show raw.where(block: true): block.with(fill: luma(245), inset: 8pt, radius: 3pt, width: 100%)
#show table: set text(size: 9.5pt)
#show math.equation: set text(font: "Lete Sans Math")

#align(center)[
  #text(size: 17pt, weight: "bold")[SourceEstimation.jl]\
  #text(size: 13pt)[Points d'implémentation importants]
]

#v(0.5em)

Ce document résume les choix d'implémentation décrits dans `docs/design/api.md`. Il insiste sur ce qui n'est pas visible depuis l'API : la disposition des données, la résolution des systèmes linéaires, la stratégie de résolution des problèmes espace–fréquence et espace–temps, et le
calcul de l'incertitude a posteriori. La dernière section discute le tirage des échantillons par perturbation–optimisation et sa relation avec un tirage gaussien « standard ».

#outline(depth: 2, indent: auto)

= Principes généraux

== Une interface, du dispatch

Toute l'API tient en `solve(problème, formulation, algorithme)`. Le choix des méthodes se fait par _multiple dispatch_ sur les types ; les branchements à l'exécution sont réservés à la validation des arguments, aux tests d'arrêt et à quelques cas particuliers (régime sur/sous-déterminé à la construction du cache, cas purement quadratique).

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([*Axe*], [*Types*], [*Ce qui change*]),
  [Opérateur direct], [`DenseOperator`, `BlockDiagonalOperator`, `SparseOperator`], [stockage et résolution des systèmes linéaires],
  [Régime], [`Overdetermined`, `Underdetermined`], [équations normales ou identité de Woodbury],
  [Opérateur de lissage], [`IdentitySmoothing`, `InvertibleSmoothing`, `SingularSmoothing`], [chemin diagonal, forme standard $H L^(-1)$, équations normales],
  [Terme], [`Lq`, `Lpq` (avec ou sans `Lc`), `Groups`], [poids IR, couplage des échantillons],
  [Formulation], [`Additive`, `Multiplicative`, `Bayesian`], [mise à jour du paramètre de régularisation],
  [Corps des données], [`RealField`, `ComplexField`], [lois réelles ou complexes du modèle bayésien],
)

== Organisation du code

- `types/` : définitions des types et constructeurs (l'API publique) ;
- `utils/` : briques numériques (systèmes linéaires, spectres, poids IR, modèle a posteriori…) ;
- `solvers/` : une méthode `solve_joint` par couple (formulation, algorithme), et `strategy.jl` qui choisit la stratégie de résolution (§4) et porte le `solve` public.

= Disposition des données et régularisation

== Lignes = espace, colonnes = échantillons

L'inconnue est toujours une matrice $U$ de taille $n_r times n_s$ : les lignes sont les points de reconstruction, les colonnes les fréquences ou les instants ($n_s = 1$ pour une seule fréquence). Les vecteurs globaux sont empilés échantillon par échantillon, $u = "vec"(U)$, soit $u_(i + (j-1) n_r) = U_(i j)$. C'est aussi l'ordre de MSSP 2018 (éq. 4).

== Termes de régularisation

$ R(U) = sum_g R_g (L_g U_g L_(c,g)^T), quad
  "Lq" : norm(V)_q^q = sum_(i,j) |V_(i j)|^q, quad
  "Lpq" : norm(V)_(p,q)^q = sum_i (sum_j |V_(i j)|^p)^(q slash p) $

`L` agit sur les lignes (espace), `Lc` sur les colonnes (fréquences ou temps). Sur $"vec"(U)$, l'opérateur global est un produit de Kronecker, grâce à l'identité $(A times.o B) "vec"(X) = "vec"(B X A^T)$ :

$
(I_(n_s) times.o L) (L_c times.o I_(n_r)) = L_c times.o L, quad
  (L_c times.o L) "vec"(U) = "vec"(L U L_c^T).
$

C'est exactement l'opérateur $overline(L) = overline(L)_s overline(L)_f$ de MSSP 2018 (éq. 15), avec $L = L_s$ et $L_c = L_f$ ; la norme de l'éq. (14) porte sur les lignes de $L_s F L_f^T$.

Les opérateurs $L_g$ des groupes sont assemblés une fois en une matrice creuse globale ; si tous valent l'identité (détecté par le type), aucune matrice n'est créée et le chemin diagonal est conservé.

= Algorithme IR

== Point fixe

À l'itération $k+1$, avec $V = L U^((k))$ : (1) poids $W$ et pénalités $P_g = R_g (V_g)$ ; (2) paramètre(s) de régularisation ; (3) résolution de

$ (H^H H + L^H D L) u = H^H y, quad D = lambda W. $

Le critère d'arrêt est la variation relative $norm(U^((k+1)) - U^((k)))_1  slash  norm(U^((k)))_1$.

== Poids issus de la condition d'optimalité

Les poids vérifient $partial R  slash  partial overline(v)_(i j) = W_(i j) v_(i j)$ :

$
"Lq" : W_i = q/2 max(epsilon, |v_i|)^(q-2), quad
"Lpq" : W_(i j) = q/2 max(epsilon^p, S_i)^(q slash p - 1) max(epsilon, |v_(i j)|)^(p-2), quad S_i = sum_j |v_(i j)|^p. $

$epsilon$ est fixé une fois pour toutes (quantile 5 % de $|L U^((0))|$, par groupe). `reweight!` remplit $W$ *en place* à travers des vues (vue sur les lignes du groupe, puis vue repliée pour un problème vectorisé) et renvoie la pénalité $R_g$.

== Paramètre de régularisation

- OMR : $alpha = norm(y - H u)^2  slash  sum_g P_g$ ; MPMR : $alpha_g = norm(y - H u)^2  slash  P_g$ ;
- additive : $lambda$ fixe, un $lambda$ par groupe, ou choisi à chaque itération (L-curve, GCV, estimateur bayésien) sur le spectre du problème pondéré.

*Cas quadratique.* Avec `Additive` et uniquement des termes `Lq(2.0)`, les poids sont constants ($W = 1$) : la solution initiale de Tikhonov est la solution finale, et `solve` fait une seule résolution linéaire, sans itération.

= Stratégie de résolution espace–fréquence et espace–temps

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([*Cas*], [*Stratégie*], [*Résolution*]),
  [Espace–fréquence, `Lq` seul], [`PerSample`], [un problème indépendant par fréquence, en parallèle (α, λ, convergence propres)],
  [Espace–fréquence, `Lpq` sans `Lc`], [`Joint` (blocs)], [poids et α communs, mais $H^H H + D$ bloc-diagonal : résolution fréquence par fréquence en parallèle],
  [Espace–temps dense, sans `Lc`], [`Joint` (dense)], [système dense de taille $n_r n_t$, opérateur $I times.o L$],
  [`Lpq` avec `Lc`], [`Assembled`], [problème vectorisé automatiquement ($H$ bloc-diagonal creux, groupes `sample_rows`, opérateur $L_c times.o L$), Cholesky creux],
)

Sans `Lc`, la structure bloc-diagonale est toujours exploitée et le couplage entre fréquences ne passe que par les poids. Avec `Lc`, $L^H D L$ relie les fréquences voisines (système tridiagonal par blocs pour un `Lc` à 3 points) : on bascule vers la résolution creuse assemblée.

= Systèmes linéaires

== Régime et identité de Woodbury

Le régime est un paramètre de type du cache ; la méthode est choisie par dispatch :

#table(
  columns: (auto, auto, 1fr, auto),
  align: left,
  table.header([*Régime*], [*`L`*], [*Méthode*], [*Taille*]),
  [$n_m >= n_r$], [quelconque], [équations normales], [$n_r$],
  [$n_m < n_r$], [$I$], [Woodbury], [$n_m$],
  [$n_m < n_r$], [carrée inversible], [Woodbury sur $tilde(H) = H L^(-1)$], [$n_m$],
  [$n_m < n_r$], [rectangulaire], [équations normales ($L^H D L$ singulière)], [$n_r$],
)

$ (H^H H + D)^(-1) H^H y = D^(-1) H^H (I + H D^(-1) H^H)^(-1) y, quad
  (H^H H + L^H D L)^(-1) = L^(-1) (tilde(H)^H tilde(H) + D)^(-1) L^(-H). $

Pour un opérateur creux (problème assemblé), on utilise CHOLMOD en réutilisant l'analyse symbolique d'une itération à l'autre.

== Allocations

Un système linéaire (`DenseSystem`, `BlockSystem`, `SparseSystem`) est construit une fois par résolution : il précalcule ce qui ne dépend pas des poids ($H^H H$, transposée de l'opérateur de lissage global, forme standard $H L^(-1)$) et préalloue ses tampons. Les boucles utilisent des versions en place (`weighted_solve!`, `residual_norm2!`, `apply!`) et échangent les
tampons `U`/`Unew`.

== Bruit

La covariance du bruit $R$ (connue, diagonale) est traitée une fois pour toutes par blanchiment : $H <- R^(-1 slash 2) H$, $y <- R^(-1 slash 2) y$. Aucune inverse n'est jamais calculée explicitement : $A^(-1)$ est toujours appliquée par une factorisation.

= Choix du paramètre et initialisation

Toutes les règles travaillent sur le spectre du problème pondéré
$min norm(y - H u)^2 + lambda norm(W^(1 slash 2) L u)^2$ : SVD de $H W^(-1 slash 2)$ si $L = I$, SVD généralisée de $(H, W^(1 slash 2) L)$ sinon. On minimise un critère en $t = log ell$ ($lambda = ell^2$) — $-kappa$ pour la L-curve, la fonctionnelle GCV, ou $J_"MAP" (ell) = sum_i log(gamma_i^2 + ell^2) + (M-2) log(1/M sum_i |beta_i|^2  slash  (gamma_i^2 + ell^2))$ — par une recherche sur grille affinée par la méthode de Brent.

L'initialisation heuristique (MSSP 2017, MSSP 2018 éq. 26–27) prend la moyenne géométrique des valeurs propres extrêmes de $A = (H L^(-1))^H (H L^(-1))$ : valeurs exactes (SVD) pour un opérateur dense, estimations par la méthode de la puissance (`normest`, `condest` en norme 2) pour un opérateur creux, $A$ et $A^(-1)$ étant appliqués comme opérateurs sans être assemblés.

= Modèle bayésien

#table(
  columns: (auto, 1fr, 1fr),
  align: left,
  table.header([], [*Réel*], [*Complexe (circulaire)*]),
  [Vraisemblance], [$(tau_n  slash  2 pi)^(N slash 2) exp(-tau_n norm(r)^2  slash  2)$], [$(tau_n  slash  pi)^N exp(-tau_n norm(r)^2)$],
  [A priori], [$prop exp(-(tau_s  slash  q) |v|^q)$], [$prop exp(-tau_s |v|^q)$],
  [MAP en $u$], [$1/2 norm(y - H u)_R^2 + sum_g lambda_g / q_g norm(L_g u_g)_q^q$], [$norm(y - H u)_R^2 + sum_g lambda_g norm(L_g u_g)_q^q$],
)

avec $lambda_g = tau_(s,g)  slash  tau_n$. L'algorithme IR résout la forme équivalente $norm(y - H u)^2 + sum_g lambda_(u,g) norm(L_g u_g)_q^q$, avec $lambda_(u,g) = 2 lambda_g  slash  q_g$ en réel et $lambda_g$ en complexe. Les hyperparamètres sont optimisés (`Estimated()`) ou fixés ; par défaut $tau_n = 1$ quand $R$ est fourni. `ElementWise()` (RVR) utilise une précision par
composante, $tau_i = (k slash 2 + a - 1)  slash  (phi |v_i|^2 + b)$.

= Incertitude a posteriori

== Précision plutôt que covariance

Au point fixe, la loi a posteriori est approchée par une gaussienne de moyenne $u$ et de précision

$ Q = tau_n A, quad A = H^H H + L^H D L, $

exacte pour $q = 2$. `solve(...; precision = true)` renvoie $Q$ (`sol.precision`) et un modèle a posteriori (`sol.posterior`) contenant une factorisation de $A$ choisie comme pour la résolution (Cholesky dense, Woodbury, Cholesky creux, ou une par bloc). $Q$ est la matrice déjà
construite par le solveur ; elle est creuse dans le cas assemblé, alors que la covariance $Q^(-1)$ est dense (1,3 Go pour 21 points × 451 fréquences). $Q^(-1)$ n'est jamais formée par défaut.

== Écarts-types et covariances marginales : solves par blocs

Le bloc diagonal de l'échantillon $j$ (fréquence ou instant) s'obtient en appliquant $A^(-1)$ aux $n_r$ colonnes correspondantes de l'identité :

$ Sigma_(j j) = 1/tau_n (A^(-1) E_j)_([j], :), $

où $E_j$ contient les colonnes de l'identité associées aux lignes $[j]$ de l'échantillon $j$. On obtient les covariances marginales exactes sans former la matrice $n times n$ (mémoire $n dot n_r$ pendant le calcul, $n_f n_r^2$ pour le résultat). `posterior_std` prend la racine de leur diagonale ; `posterior_covariance(sol; full = true)` donne la covariance complète.

== Tirage des échantillons : perturbation–optimisation

`posterior_samples` tire

$ u_s = u + 1/sqrt(tau_n) A^(-1) (H^H epsilon_1 + L^H D^(1 slash 2) epsilon_2), quad
  epsilon_1 tilde cal(N)(0, I_m), quad epsilon_2 tilde cal(N)(0, I_p) $

(lois normales complexes circulaires si les données sont complexes). Comme $epsilon_1$ et $epsilon_2$ sont indépendants,
$"Cov"(H^H epsilon_1 + L^H D^(1 slash 2) epsilon_2) = H^H H + L^H D L = A$, donc

$ "Cov"(u_s) = 1/tau_n A^(-1) A A^(-1) = (tau_n A)^(-1) = Q^(-1). $

Le nom vient de l'interprétation : $u_s$ est la solution du même problème de Tikhonov pondéré dont les données et l'a priori ont été perturbés,

$ u_s = limits(op("argmin"))_x norm(y + epsilon_1  slash  sqrt(tau_n) - H x)^2
        + norm(D^(1 slash 2) L x - epsilon_2  slash  sqrt(tau_n))^2 . $

== Est-ce un tirage standard ?

*Oui, en loi : c'est exactement un tirage de $cal(N)(u, Q^(-1))$*, comme le tirage standard $u_s = u + S z$, $z tilde cal(N)(0, I)$, avec $S S^H = Q^(-1)$. La différence porte uniquement sur la *racine carrée* utilisée et sur ce qu'il faut calculer.

Écrivons $A = B^H B$ avec la matrice rectangulaire, connue sans calcul,

$ B = mat(H; D^(1 slash 2) L) in bb(C)^((m + p) times n). $

- Le tirage standard utilise une racine *carrée* : le facteur de Cholesky $A = R^H R$, et
  $u_s = u + R^(-1) z  slash  sqrt(tau_n)$. Il faut factoriser $A$ (ou $Q$), puis faire une résolution triangulaire par tirage.
- La perturbation–optimisation utilise la racine *rectangulaire* $B$ :
  $u_s = u + A^(-1) B^H epsilon  slash  sqrt(tau_n)$, $epsilon tilde cal(N)(0, I_(m+p))$.

Les deux sont équivalents. Avec une factorisation QR $B = Q_B R$ (où $R$ est le facteur de Cholesky de $A$) :

$ A^(-1) B^H epsilon = R^(-1) R^(-H) R^H Q_B^H epsilon = R^(-1) (Q_B^H epsilon), $

et $Q_B^H epsilon tilde cal(N)(0, I_n)$ puisque $Q_B$ a des colonnes orthonormées : la perturbation–optimisation *est* le tirage standard $R^(-1) z$, sans jamais calculer $R$ ni $Q_B$. C'est le sens de la remarque : on ne change pas la loi, on change la façon de l'atteindre.

Ce qui justifie ce choix dans le paquet :

+ *Aucune factorisation supplémentaire* : il suffit d'appliquer $A^(-1)$, avec la factorisation déjà utilisée par le solveur, quelle qu'elle soit. Dans le régime sous-déterminé, le solveur n'a pas de Cholesky de $A$ ($n_r times n_r$) mais une factorisation de $I + H D^(-1) H^H$
  ($n_m times n_m$) : le tirage standard imposerait un Cholesky $n_r times n_r$ (coût $n_r^3$), la perturbation–optimisation coûte $O(n_m n_r)$ par tirage (Woodbury). Il en va de même pour la forme standard $H L^(-1)$ et pour les factorisations par bloc.
+ *Pas de covariance* : on ne forme ni $Q^(-1)$ ni sa racine ; dans le cas creux assemblé, les tirages sont joints (corrélations entre fréquences conservées) avec la seule factorisation creuse de $A$.
+ *Un seul code* pour tous les régimes et tous les opérateurs (dense, blocs, creux), puisqu'il ne repose que sur `inverse_apply`.

En contrepartie, quand un facteur de Cholesky $R$ de $A$ est déjà disponible (équations normales, CHOLMOD), le tirage standard $R^(-1) z$ est un peu moins coûteux par tirage : une résolution triangulaire au lieu de deux, plus les produits $H^H epsilon_1$ et $L^H D^(1 slash 2) epsilon_2$. La différence est faible devant le coût de la factorisation, et le
choix d'une méthode unique a été privilégié. Les deux méthodes donnent la même loi.

*Références*

+ F. Orieux, O. Féron, J.-F. Giovannelli, « Sampling high-dimensional Gaussian distributions for general linear inverse problems », IEEE Signal Processing Letters 19(5):251–254, 2012. (C'est l'algorithme que j'utilise).

+ G. Papandreou, A. Yuille, « Gaussian sampling by local perturbations », NIPS 23, 2010, p. 1858–1866. Même principe, appliqué aux champs gaussiens.

+ C. Gilavert, S. Moussaoui, J. Idier, « Efficient Gaussian sampling for solving large-scale inverse problems using MCMC », IEEE Trans. Signal Processing 63(1):70–80, 2015. Variante où l'optimisation est faite de façon approchée puis corrigée par Metropolis–Hastings.

+ M. Vono, N. Dobigeon, P. Chainais, « High-dimensional Gaussian sampling: a review and a unifying approach… », SIAM Review 64

== Intervalles de crédibilité

`credible_interval` calcule les quantiles empiriques $(1 - "level") slash 2$ et $(1 + "level") slash 2$ des tirages, composante par composante (parties réelle et imaginaire séparément en complexe).
