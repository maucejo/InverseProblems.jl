#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Bibliography])
#metadata((description: "Bibliography of the methods implemented in SourceEstimation.jl.")) <website-metadata>

#title()

The methods implemented in SourceEstimation.jl have been developed and published in the following papers. The table maps each feature of the package to its main reference.

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([*Feature*], [*API*], [*Reference*]),
  [Generalized Tikhonov regularization, $ell_q$ priors], [`Additive`, `Lq`], [@Aucejo2014],
  [Bayesian regularization with local priors], [`Bayesian`, `Groups`], [@Aucejo2016],
  [Multiplicative regularization (OMR), heuristic initialization], [`Multiplicative`, `Heuristic`], [@Aucejo2017 @Aucejo2018a],
  [Space–frequency mixed norms, operators along the frequencies], [`Lpq`, `Lc`, `space_sample_operator`], [@Aucejo2018b],
  [Full Bayesian inference], [`Bayesian`], [@Aucejo2018c],
  [Optimal Bayesian regularization (hierarchical MAP)], [`Bayesian`, `MAP`], [@Aucejo2019a],
  [Space–time regularization], [space–time `InverseProblem`], [@Aucejo2019b],
  [General iteratively reweighted algorithm], [`IR`], [@Aucejo2019c],
  [Multi-parameter multiplicative regularization (MPMR)], [`multiparameter = true`], [@Aucejo2020a],
  [Multiplicative mixed-norm algorithm], [`Multiplicative(Lpq(…))`], [@Aucejo2020b],
  [Generalized multiplicative regularization], [`Multiplicative`], [@Aucejo2021],
  [L-curve, GCV, Bayesian estimator], [`LCurve`, `GCV`, `BayesEstimator`], [@Hansen1993 @Golub1979 @Pereira2015],
  [$lambda_"max"$ initialization], [`LambdaMax`], [@Kim2007],
  [Relevance vector regression], [`ElementWise`], [@Tipping2001],
  [Perturbation–optimization sampling], [`posterior_samples`], [@Orieux2012 @Papandreou2010],
  [Trueness, precision and accuracy], [`relative_error`, `precision_metric`, `accuracy_metric`], [@ISO5725 @VIM2012],
)

The general background on discrete ill-posed problems can be found in @Hansen2010, and recent developments on high-dimensional Gaussian sampling in @Gilavert2015 @Vono2022.

#bibliography("/refs.bib", title: [All references], style: "ieee", full: true)
