# CML Algorithm

`cml()` applies the CML (Coordinated Multi-Neighborhood Learning)
algorithm to a dataset over certain target neighborhoods, which may be
provided by the user through the true DAG or may be estimated by a
Markov Blanket estimation algorithm.

## Usage

``` r
cml(
  data = NULL,
  true_dag = NULL,
  targets,
  node_names = NULL,
  lmax = 3,
  tol = 0.01,
  mb_tol = 0.05,
  method = "MMPC",
  test = "testIndFisher",
  verbose = TRUE
)
```

## Arguments

- data:

  A data matrix or data frame, or `NULL` for the population (oracle)
  version.

- true_dag:

  A 0/1 adjacency matrix of the true DAG, or `NULL` to estimate the
  neighborhoods from `data` instead.

- targets:

  A vector of 1-based target node indices.

- node_names:

  Character vector of node names; defaults to `V0`, `V1`, ... if not
  supplied.

- lmax:

  Maximum size of the conditioning set considered during the skeleton
  search.

- tol:

  Significance level for the conditional independence tests used by the
  skeleton search.

- mb_tol:

  Significance level used for Markov Blanket estimation (see
  `get_all_mbs()`); only relevant when `true_dag` is `NULL`.

- method:

  Markov Blanket estimation algorithm (see `get_mb()`); only relevant
  when `true_dag` is `NULL`.

- test:

  The conditional independence test to use.

- verbose:

  Whether to provide detailed output.

## Value

A list with the estimated adjacency matrix (`amat`), separating sets
(`S`), the number of conditional independence tests used (`NumTests`,
`MBNumTests`), the nodes considered (`Nodes`, in 1-based R numbering),
which FCI orientation rules were used (`RulesUsed`), timing information,
the reference DAG used (`referenceDAG`), the estimated Markov Blanket
list (`mbList`, if applicable), and, when `data` was supplied, the
pre-scaling column means/covariance (`data_means`, `data_cov`).

## Details

Exactly one of `data` or `true_dag` (population/oracle version) must
ultimately be usable to identify each target's neighborhood: if
`true_dag` is `NULL`, the neighborhoods are first estimated from `data`
via `get_all_mbs()`; otherwise `true_dag` is used directly (with `data`,
if also supplied, used only for the conditional independence tests
themselves – the "semi-sample" version).
