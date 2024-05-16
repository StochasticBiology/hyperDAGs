# hyperDAGs
DAG skeletons on hypercubic transition networks

![image](https://github.com/StochasticBiology/hyper-DAGs/assets/50171196/44b0b540-2c85-442c-9f19-ccbf4ae042cb)
![image](https://github.com/StochasticBiology/hyperDAGs/assets/50171196/3066dd2f-b182-4bc1-93ea-a72677bfeb68)


Given a dataset of binary string observations, can we construct an arborescence with minimal branching that is a subset of the hypercubic digraph and hits all the vertices corresponding to datapoints? Yes -- the "Gutin algorithm" will do this in polynomial time https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 , https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5 .

Given a dataset of *paired* binary string observations, can we construct a DAG with minimal branching that is a subset of the hypercubic digraph, hits all the vertices corresponding to datapoints, and has a path going from ancestor to descendant for each pair? That is what we are trying to find here. `paired-observations.R` currently contains Kostas' implementation of the Gutin algorithm and Iain's simple, bad approximation on top of this. 
