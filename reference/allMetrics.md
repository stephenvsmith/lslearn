# Score an Estimated Graph Against a Ground-Truth DAG

`allMetrics()` bundles every graph-comparison metric in this file
(skeleton recovery, v-structure recovery, parent-recovery accuracy for
`targets`, cross-neighborhood ancestral-edge recovery, and the overall
F1 score) into a single one-row data frame, for evaluating the output of
[`cml()`](https://stephenvsmith.github.io/lslearn/reference/cml.md)/[`snl()`](https://stephenvsmith.github.io/lslearn/reference/snl.md)
against a known ground-truth network.

## Usage

``` r
allMetrics(
  est,
  ref_graph,
  targets,
  true_dag,
  nbhd,
  verbose = FALSE,
  algo = "pc",
  which_nodes = ""
)
```

## Arguments

- est:

  The estimated adjacency matrix (e.g.
  [`cml()`](https://stephenvsmith.github.io/lslearn/reference/cml.md)'s
  or
  [`snl()`](https://stephenvsmith.github.io/lslearn/reference/snl.md)'s
  `amat` result).

- ref_graph:

  The ground-truth adjacency matrix restricted to the same nodes as
  `est`, used for the skeleton/v-structure/parent-recovery metrics.

- targets:

  0-based target node indices (in `est`'s numbering), used for the
  parent-recovery and F1 metrics.

- true_dag:

  The full ground-truth adjacency matrix (over all nodes in the original
  network, not just the ones in `est`), used to check ancestral
  relationships for `interNeighborhoodEdgeMetrics()`.

- nbhd:

  A mapping from `est`'s node numbering to `true_dag`'s node numbering
  (i.e. `est`'s neighborhood, in `true_dag`'s indices).

- verbose:

  Whether to print detailed output.

- algo:

  A short label for the algorithm being evaluated (e.g. `"cml"`,
  `"snl"`, `"pc"`), used as a prefix on most column names.

- which_nodes:

  An additional label distinguishing this evaluation run (e.g. which
  subset of nodes it covers), used as a prefix on the
  skeleton/v-structure column names alongside `algo`.

## Value

A one-row data frame with skeleton false positives/negatives/ true
positives, v-structure false positives/negatives/true positives,
parent-recovery false positives/negatives/true positives/potentials,
ancestral-edge recovery counts, and the overall F1 score – all column
names prefixed with `algo` (and `which_nodes` for the skeleton/
v-structure columns).
