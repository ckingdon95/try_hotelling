using BlackBoxOptim
using Mimi
using MimiDICE2016

include("components/welfare1.jl")
include("components/welfare2.jl")
include("components/capped_emissions.jl")

years = 2015:5:2300
nyears = length(years)

# Returns a version of DICE2016 without any climate or damage components, with a welfare
#   component that sets total utility equal to a penalty when cumulative emissions in the final year
#   exceed the fosslim value. TFP is taken as exogenous from the provided values
function _get_model1(tfp, fosslim, penalty)
    m = MimiDICE2016.get_model()
    set_dimension!(m, :time, years)
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

# Returns a version of DICE2016 without any climate or damage components, with a welfare
#   component that sets utility in each period equal to a penalty when cumulative emissions in that year
#   exceed the fosslim value. TFP is taken as exogenous from the provided values
function _get_model2(tfp, fosslim, penalty)
    m = MimiDICE2016.get_model()
    set_dimension!(m, :time, years)
    delete!(m, :totalfactorproductivity)
    set_param!(m, :grosseconomy, :AL, tfp)

    update_param!(m, :MIU, zeros(100))

    replace_comp!(m, welfare2, :welfare)
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

# Returns a version of DICE2016 without any climate or damage components, with a modified emissions component,
#   where cumulativve emissions cannot exceed the fosslim. Abatement is forced to stay below, and there are abatamente costs.
function _get_model3(tfp, fosslim)
    m = MimiDICE2016.get_model()
    set_dimension!(m, :time, years)
    delete!(m, :totalfactorproductivity)
    set_param!(m, :grosseconomy, :AL, tfp)

    replace_comp!(m, capped_emissions, :emissions)
    set_param!(m, :emissions, :fosslim, fosslim)
    set_param!(m, :emissions, :initial_MIU, zeros(58))
    connect_param!(m, :neteconomy, :MIU, :emissions, :final_MIU)

    delete!(m, :co2cycle)
    delete!(m, :radiativeforcing)
    delete!(m, :climatedynamics)
    delete!(m, :damages)

    set_param!(m, :neteconomy, :DAMAGES, zeros(nyears))

    return m
end

function _get_model4(tfp)
    m = MimiDICE2016.get_model()
    set_dimension!(m, :time, years)
    delete!(m, :totalfactorproductivity)
    set_param!(m, :grosseconomy, :AL, tfp)
    update_param!(m, :MIU, zeros(100))
    set_param!(m, :neteconomy, :DAMAGES, zeros(nyears))

    return m
end

# Returns a version of DICE2016 where MIU in every year is set to 1 - fosslim/CCA[end]
function _get_model5(tfp, fosslim)
    m = MimiDICE2016.get_model()
    set_dimension!(m, :time, years)
    delete!(m, :totalfactorproductivity)
    set_param!(m, :grosseconomy, :AL, tfp)
    update_param!(m, :MIU, zeros(100))
    run(m)

    update_param!(m, :MIU, ones(100) * (1 - fosslim / m[:emissions, :CCA][end]))
    return m
end

function _run_optimization1(tfp, fosslim, penalty, max_steps)
    m = _get_model1(tfp, fosslim, penalty)
    run(m)

    function _eval(x)
        m.mi.md.external_params[:MIU].values.data = x
        run(m)
        return -m[:welfare, :UTILITY]
    end

    res = bboptimize(
        _eval; 
        SearchRange=(0., 1.0), 
        NumDimensions=nyears, 
        Method=:adaptive_de_rand_1_bin_radiuslimited, 
        MaxSteps=max_steps
    )

    return m, res
end

function _run_optimization2(tfp, fosslim, penalty, max_steps)
    m = _get_model2(tfp, fosslim, penalty)
    run(m)

    function _eval(x)
        m.mi.md.external_params[:MIU].values.data = x
        run(m)
        return -m[:welfare, :UTILITY]
    end

    res = bboptimize(
        _eval; 
        SearchRange=(0., 1.0), 
        NumDimensions=nyears, 
        Method=:adaptive_de_rand_1_bin_radiuslimited, 
        MaxSteps=max_steps
    )

    return m, res
end

function _run_optimization3(tfp, fosslim, max_steps)
    m = _get_model3(tfp, fosslim)
    run(m)

    function _eval(x)
        m.mi.md.external_params[:initial_MIU].values.data = x
        run(m)
        return -m[:welfare, :UTILITY]
    end

    res = bboptimize(
        _eval; 
        SearchRange=(0., 1.0), 
        NumDimensions=nyears, 
        Method=:adaptive_de_rand_1_bin_radiuslimited, 
        MaxSteps=max_steps
    )

    return m, res
end

function _run_optimization4(tfp, fosslim, max_steps, penalty_factor)
    m = _get_model4(tfp)
    run(m)

    function _eval(x)
        m.mi.md.external_params[:MIU].values.data = x
        run(m)

        penalty = penalty_factor * max(m[:emissions, :CCA][end] - fosslim, 0)

        return -(m[:welfare, :UTILITY] - penalty)
    end

    res = bboptimize(
        _eval; 
        SearchRange=(0., 1.0), 
        NumDimensions=nyears, 
        Method=:adaptive_de_rand_1_bin_radiuslimited, 
        MaxSteps=max_steps
    )

    return m, res
end
