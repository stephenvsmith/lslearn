# Structure of the Andes Bayesian Network

The 0/1 adjacency matrix of the "Andes" tutoring-system Bayesian network
(Conati et al., 1997), a 223-node, 338-edge example network used to
exercise graph-metric functions on a network much larger than `asiaDAG`.

## Usage

``` r
data(andesDAG)
```

## Format

A 223x223 numeric matrix with row and column names given by the
network's node names. `andesDAG[i, j] == 1` means node `i` is a parent
of node `j`.

## Source

Ported from the CML package (<https://github.com/stephenvsmith/CML>).
