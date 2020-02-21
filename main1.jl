using BlackBoxOptim
using Mimi
using MimiDICE2016

include("welfare1.jl")

years = 2015:5:2300
nyears = length(years)
fosslim = 6000.     # 6,000 GtC
high_tfp = range(5, 150, length=nyears)
penalty = -1e10

function _get_model1(; tfp, fosslim)
    m = MimiDICE2016.get_model()
    set_dimension!(m, :time, 2015:5:2300)
    delete!(m, :totalfactorproductivity)
    set_param!(m, :grosseconomy, :AL, tfp)

    update_param!(m, :MIU, zeros(100))

    replace_comp!(m, welfare1, :welfare)
    connect_param!(m, :welfare => :CCA, :emissions => :CCA)
    set_param!(m, :welfare, :fosslim, fosslim)
    set_param!(m, :welfare, :penalty, penalty)

    delete!(m, :co2cycle)
    delete!(m, :radiativeforcing)
    delete!(m, :climatedynamics)
    delete!(m, :damages)

    set_param!(m, :neteconomy, :DAMAGES, zeros(nyears))

    return m
end

m1 = _get_model1(tfp=high_tfp, fosslim=fosslim)
run(m1)

function _eval1(x)
    m1.mi.md.external_params[:MIU].values.data = [0, x...]
    run(m1)
    return -m1[:welfare, :UTILITY]
end

res = bboptimize(
    _eval1; 
    SearchRange=(0., 1.0), 
    NumDimensions=nyears-1, 
    Method=:adaptive_de_rand_1_bin_radiuslimited, 
    MaxSteps=100_000
)

# opt_miu = best_candidate(res)

run(m1)
explore(m1, title = "m1")