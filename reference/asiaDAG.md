# Structure of the Asia Bayesian Network

The 0/1 adjacency matrix of the well-known "Asia" (lung cancer) Bayesian
network (Lauritzen and Spiegelhalter, 1988), an 8-node example network
used as ground truth in lslearn's Markov Blanket estimation tests.

## Usage

``` r
data(asiaDAG)
```

## Format

An 8x8 numeric matrix with row and column names `asia`, `tub`, `smoke`,
`lung`, `bronc`, `either`, `xray`, `dysp`. `asiaDAG[i, j] == 1` means
node `i` is a parent of node `j`.

## Source

Ported from the CML package (<https://github.com/stephenvsmith/CML>).
