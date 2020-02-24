
include("src/utils.jl")

fosslim = 6000.     # 6,000 GtC
penalty = -1e10
high_tfp = range(5, 150, length=nyears)
max_steps = 500_000

# m1, res1 = _run_optimization1(high_tfp, fosslim, penalty, max_steps)
m2, res2 = _run_optimization2(high_tfp, fosslim, penalty, max_steps)
m3, res3 = _run_optimization3(high_tfp, fosslim, max_steps)
m4, res4 = _run_optimization4(high_tfp, fosslim, max_steps, 1e10)

# explore(m1, title = "m1 - Penalty replaces total UTILITY")
explore(m2, title = "m2 - Penalty in each year that exceeds fosslim")
explore(m3, title = "m3 - Cora's way (forced abatement beyond fosslim)")
explore(m4, title = "m4 - David's way (penalty factor)")

# Gives the same solution for MIU to 3 decimal places
isapprox.(m2[:emissions, :MIU][1:end-1], m3[:emissions, :final_MIU][1:end-1], atol=1e-3)
isapprox.(m2[:emissions, :MIU][1:end-1], m4[:emissions, :MIU][1:end-1], atol=1e-3)

