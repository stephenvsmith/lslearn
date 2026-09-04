# CML Algorithm (MAG-only output)

`cml_mag()` applies the CML algorithm identically to
[`cml()`](https://stephenvsmith.github.io/lslearn/reference/cml.md),
except it skips the within-neighborhood mixed-graph conversion rules,
leaving the result as an ancestral (MAG-style) graph rather than
converting it to CML's own neighborhood notation.

## Usage

``` r
cml_mag(
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

Same shape as
[`cml()`](https://stephenvsmith.github.io/lslearn/reference/cml.md).
