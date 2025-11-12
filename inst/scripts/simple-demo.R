#### basic example for demonstration

# simple illustrative dataset
ancnames  = c("10000","01000","00000","10010","01011","10110")
descnames = c("11110","11110","01011","11110","01111","11111")
# get STSD
s.dag = simplest_DAG(ancnames, descnames)
# plot and summarise: A is SSA, B is STSD, numbers are branching counts
plot_stage_2(s.dag)
fit_properties(s.dag)
transitions_spanned(s.dag$best.graph, ancnames, descnames)
