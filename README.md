# hyperDAGs
DAG skeletons on hypercubic transition networks

Given a dataset of binary string observations, can we construct an arborescence with minimal branching that is a subset of the hypercubic digraph and hits all the vertices corresponding to datapoints? Yes -- the "Gutin algorithm" will do this in polynomial time https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 , https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5 .

Given a dataset of *paired* binary string observations, can we construct a DAG with minimal branching that is a subset of the hypercubic digraph, hits all the vertices corresponding to datapoints, and has a path going from ancestor to descendant for each pair? That is what we are trying to find here. 

`paired-observations.R` currently contains Kostas' implementation of the Gutin algorithm and Iain's developments on top of this, attempting to compute this picture and also "simplify" the result in the sense that branch points are shifted down the graph as much as possible. Plot and helper functions also included. `wrapper.R` and the accompanying datafiles provide some scientific case studies, curating input data and wrapping the analysis.

Essentials
----

`paired-observations.R` contains the essential functions for the workflow:

| Function | Description |
| ----|-----|
| `simplest.DAG(ancnames, descnames)` | Estimates a simplest DAG linking ancestors `ancnames` to descendants `descnames`, where both are vectors of character strings of the same length. Runs `simplest.arborescence` then adds and prunes edges to connect required relationships. Returns a fitted model structure. |
| `simplest.arborescence(ancnames, descnames)` | Stage 1: Use Gutin's approach to get the simplest arborescence for the union of the ancestors and descendants, then rewire to attempt to maximise layer sum. Returns a fitted model structure. |
| `branching.count(g)` | Returns excess branching of graph `g` |
| `layer.sum(g)` | Returns layer sum of graph `g` |
| `fit.properties(fit)` | Outputs various statistics of a model structure `fit` |
| `transitions.spanned(g, ancnames, descnames)` | Queries whether or not graph `g` connects the ancestor-descendant relationships required |
| `plot.stage.1(fit)`, `plot.stage.2(fit)` | Plot visualisations of stage 1 and stage 2 solutions in fitted model `fit` |

Fitted model structure:

| Element | Description |
|-----|-----|
| `len` | Length of bitstrings (numeric) |
| `names` | Observations (vector of characters) |
| `raw.graph` | Output from Stage 1 (Gutin's algorithm) (igraph) |
| `rewired.graph` | Output from Stage 2 (rewiring for layer sum) (igraph) |
| `best.graph` | Output from Stage 3 (connecting anc-desc pairs) (igraph) |
| `raw.bc`, `rewired.bc`, `raw.ls`, `rewired.ls`, `best.bc` | [Excess] branching count and layer sum for the above graphs (numeric) |
| `dataset` | Dataframe containing ancestor and descendant columns (dataframe) |

![image](https://github.com/user-attachments/assets/f5a743cc-84fb-4cb8-a536-5c64e69bf199)

Stages 1 (Gutin's algorithm for union of data); 2 (rewiring to maximise layer sum); 3 (respecting ancestor-descendant relationships) for a simple test case



