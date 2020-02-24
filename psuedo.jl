# A "middle of the road" approach, where instead of running an optimization problem,
#   we just set MIU in every year equal to (1 - fosslim / CCA[end]) using the CCA
#   from running the model with no abatement

include("src/utils.jl")

fosslim = 6000.     # 6,000 GtC
high_tfp = range(5, 150, length=nyears)


m5 = _get_model5(high_tfp, fosslim)
run(m5)
explore(m5, title = "m5 - 'middle of the road' approach (constant MIU)")

original = Model(m5)
update_param!(original, :MIU, zeros(100))
run(original)
explore(original, title = "original high TFP model with no abatement")

# But CCA ends up exceeding fosslim in m5:
m5[:emissions, :CCA][end] > fosslim

# Because with more abatement, there are less climate damages, and YGROSS is higher
#   in later years, causing more emissions than in the run with no abatement
m5[:grosseconomy, :YGROSS] .> original[:grosseconomy, :YGROSS]