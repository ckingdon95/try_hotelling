# A "middle of the road" approach, where instead of running an optimization problem,
#   we just reduced emissions by the ratio necessary to reach fosslim in teh final year
#   (and set that emissions path exogenously, to avoid the iterating feedback problem)

include("src/utils.jl")

fosslim = 6000.     # 6,000 GtC
high_tfp = range(5, 150, length=nyears)


m5, original = _get_model5(high_tfp, fosslim, return_original=true)
run(m5)
explore(m5, title = "m5 - 'middle of the road' approach (constant MIU)")

explore(original, title = "original high TFP model with no abatement")