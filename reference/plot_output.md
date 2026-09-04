# Compare the Outputs from the Estimated and True Plots

Plots the estimated graph (restricted to the nodes actually considered
by the algorithm) side-by-side with the corresponding subgraph of the
true DAG, for visual comparison.

## Usage

``` r
plot_output(local_output, true_dag)
```

## Arguments

- local_output:

  A list of the output from one of the local learning algorithms (e.g.
  [`cml()`](https://stephenvsmith.github.io/lslearn/reference/cml.md) or
  [`snl()`](https://stephenvsmith.github.io/lslearn/reference/snl.md)),
  with `amat` and `Nodes` entries.

- true_dag:

  A matrix containing the adjacency matrix for the true DAG.
