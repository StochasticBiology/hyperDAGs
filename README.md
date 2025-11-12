# hyperDAGs
An R library and associated case study code for analysing DAG skeletons on hypercubic transition networks

For now, install with `remotes::install_github("StochasticBiology/hyperDAGs")`

Given a dataset of binary string observations, can we construct an arborescence with minimal branching that is a subset of the hypercubic digraph and hits all the vertices corresponding to datapoints? Yes -- the "Gutin algorithm" will do this in polynomial time https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 , https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5 .

Given a dataset of *paired* binary string observations, can we construct a DAG with minimal branching that is a subset of the hypercubic digraph, hits all the vertices corresponding to datapoints, and has a path going from ancestor to descendant for each pair? That is what we are trying to find here. 

In `R/`, `paired-observations.R` currently contains Kostas' implementation of the Gutin algorithm and Iain's developments on top of this, attempting to compute this picture and also "simplify" the result in the sense that branch points are shifted down the graph as much as possible. Plot and helper functions also included. 

In `inst/scripts/`, `dataset-analysis.R` and the accompanying datafiles provide some scientific case studies, curating input data and wrapping the analysis, as well as casting data in a consistent form and comparing an alternative approach in Python. `simulation-loop.R` contains simulation code and analysis to explore the behaviour of the approach in gold-standard cases. The repository also contains some datasets for demonstrations in `inst/extdata` (see References).

Essentials
----

`paired-observations.R` contains the essential functions for the workflow:

| Function | Description |
| ----|-----|
| `simplest_DAG(ancnames, descnames)` | Estimates a simplest DAG linking ancestors `ancnames` to descendants `descnames`, where both are vectors of character strings of the same length. Runs `simplest.arborescence` then adds and prunes edges to connect required relationships. Returns a fitted model structure. |
| `simplest_arborescence(ancnames, descnames)` | Stage 1: Use Gutin's approach to get the simplest arborescence for the union of the ancestors and descendants, then rewire to attempt to maximise layer sum. Returns a fitted model structure. |
| `branching_count(g)` | Returns excess branching of graph `g` |
| `layer_sum(g)` | Returns layer sum of graph `g` |
| `fit_properties(fit)` | Outputs various statistics of a model structure `fit` |
| `transitions_spanned(g, ancnames, descnames)` | Queries whether or not graph `g` connects the ancestor-descendant relationships required |
| `plot_stage_1(fit)`, `plot_stage_2(fit)` | Plot visualisations of stage 1 and stage 2 solutions in fitted model `fit` |

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


References
-----

Organelle genome data:

Giannakis, K., Arrowsmith, S.J., Richards, L., Gasparini, S., Chustecki, J.M., Røyrvik, E.C. and Johnston, I.G., 2022. Evolutionary inference across eukaryotes identifies universal features shaping organelle gene retention. Cell Systems, 13(11), pp.874-884.

Tuberculosis data:

Casali, N., Nikolayevskyy, V., Balabanova, Y., Harris, S.R., Ignatyeva, O., Kontsevaya, I., Corander, J., Bryant, J., Parkhill, J., Nejentsev, S. and Horstmann, R.D., 2014. Evolution and transmission of drug-resistant tuberculosis in a Russian population. Nature genetics, 46(3), pp.279-286.

Leukemia data:

Morita, K., Wang, F., Jahn, K., Hu, T., Tanaka, T., Sasaki, Y., Kuipers, J., Loghavi, S., Wang, S. A., Yan, Y., Furudate, K., Matthews, J., Little, L., Gumbs, C., Zhang, J., Song, X., Thompson, E., Patel, K. P., Bueso-Ramos, C. E., … Takahashi, K. (2020). Clonal evolution of acute myeloid leukemia revealed by high-throughput single-cell genomics. Nature Communications, 11(1), 5327. https://doi.org/10.1038/s41467-020-19119-8

and transition curation from:

Luo, X. G., Kuipers, J., & Beerenwinkel, N. (2023). Joint inference of exclusivity patterns and recurrent trajectories from tumor mutation trees. Nature Communications, 14(1), Article 1. https://doi.org/10.1038/s41467-023-39400-w

Ovarian cancer data:

Knutsen, T., Gobu, V., Knaus, R., Padilla‐Nash, H., Augustus, M., Strausberg, R.L., Kirsch, I.R., Sirotkin, K. and Ried, T., 2005. The interactive online SKY/M‐FISH & CGH database and the Entrez cancer chromosomes search database: linkage of chromosomal aberrations with the genome sequence. Genes, Chromosomes and Cancer, 44(1), pp.52-64.



