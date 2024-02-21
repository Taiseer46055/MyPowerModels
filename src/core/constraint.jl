###############################################################################
# This file defines commonly used constraints for power flow models
# These constraints generally assume that the model contains p and q values
# for branches flows and bus flow conservation
###############################################################################


################################### Start Taiseer Code #########################
function constraint_system_inertia(pm::AbstractPowerModel, H_min::Float64 , f_options::Dict{String, Any})

    # Extract options from the provided dictionary
    #println(options)
    #f_options = options["f"]
    disturbance = f_options["disturbance"]
    system = f_options["system"]
    area = f_options["area"]
    bus = f_options["bus"]
    rocof = f_options["rocof"]
    weighted_area = f_options["weighted_area"]

    # Reference data from the power model
    gen_data = ref(pm, :gen)
    load_data = ref(pm, :load)
    baseMVA = ref(pm, :baseMVA)
    bus_data = ref(pm, :bus)

    # Define variable for generator status (on/off)
    z = var(pm, :z_gen)

    # Nominal frequency of the network
    f0 = 50.0
    # Calculate total load in the system by summing power demand across all loads
    P_load = sum(haskey(load, "pd") ? load["pd"] : 0.0 for (_, load) in load_data)
  
    # Ensure there is a non-zero load in the system
    @assert P_load > 0 "P_load must be greater than 0"

    # Calculate the minimum required system inertia (H_min)
    # H_min = (delta_P * f0) / (P_load * 2 * rocof)
    println("H_min: ", H_min)
    # Calculate the actual system inertia (H_sys) based on generator data and status
    H_sys = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * z[i] for i in eachindex(gen_data)) / sum(gen_data[i]["pmax"] * z[i] for i in eachindex(gen_data))
    println("H_sys: ", H_sys)
    println("Pmax: ", [gen_data[i]["pmax"] for i in eachindex(gen_data)])
    println("Pmin: ", [gen_data[i]["pmin"] for i in eachindex(gen_data)])

    if system == "true"
        # Apply the system inertia constraint to the model
        println("Adding minimum system inertia constraint to ACPModel")
        JuMP.@NLconstraint(pm.model, H_sys >= H_min)
    end
    # Apply constraints based on the type of disturbance
    if disturbance == "small"
        # hier area  einfügen. Entweder nicht da oder H_min_area>=H_min oder gewichtete Mittlewert >= H_min
        if weighted_area == "equal"
            # Area-specific inertia constraints
            println("Adding minimum equal area inertia constraint to ACPModel")
            # Initialize dictionaries to store area-specific data
            areas = unique([bus_data[j]["area"] for j in keys(bus_data)])
            for area in areas
                H_area = Dict()
                # Identify generators and loads within each area
                gens_in_area = [i for i in eachindex(gen_data) if bus_data[gen_data[i]["gen_bus"]]["area"] == area]
                H_area[area] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0) / sum(gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0)
                # Apply area-specific inertia constraints
                JuMP.@NLconstraint(pm.model, H_area[area] >= H_min)
            end
        elseif weighted_area == "load"
            # Area-specific inertia constraints
            println("Adding minimum weighted area inertia constraint to ACPModel")
            areas = unique([bus_data[j]["area"] for j in keys(bus_data)])
            sum_H_area_weighted = 0
            for area in areas
                W_v = Dict()
                H_area = Dict()
                load_sum_area = Dict()
                gens_in_area = [i for i in eachindex(gen_data) if bus_data[gen_data[i]["gen_bus"]]["area"] == area]
                loads_in_area = [i for i in eachindex(load_data) if bus_data[load_data[i]["load_bus"]]["area"] == area]
                load_sum_area[area] = sum(load_data[i]["pd"] for i in loads_in_area; init=0)
                H_area[area] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0) / sum(gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0)
                W_v[area] = load_sum_area[area] / P_load
                println("W_v von area $area: ", W_v[area])
                sum_H_area_weighted  += W_v[area] * H_area[area]
            end
            println("sum_H_area_weighted: ", sum_H_area_weighted)
            JuMP.@NLconstraint(pm.model, sum_H_area_weighted >= H_min)

        elseif weighted_area == "none"

          println("weighted_area_constraint is not added to the model")  
        end
    elseif disturbance == "large"
  
        if area == "true"
            # Area-specific inertia constraints
            println("Adding minimum area inertia constraint to ACPModel")
            
            # Initialize dictionaries to store area-specific data
            H_area = Dict()
            H_min_area = Dict()
            P_load_area = Dict()
            delta_p_area = Dict()
            pg = var(pm, :pg)
            areas = unique([bus_data[j]["area"] for j in keys(bus_data)])
            for area in areas
                
                # Identify generators and loads within each area
                gens_in_area = [i for i in eachindex(gen_data) if bus_data[gen_data[i]["gen_bus"]]["area"] == area]
                buses_in_area = [i for i in keys(bus_data) if bus_data[i]["area"] == area]
                
                # Calculate load and potential delta_p for each area
                P_load_area[area] = sum(load_data[i]["pd"] for i in eachindex(load_data) if load_data[i]["load_bus"] in buses_in_area; init=0)
                P_gen_area_expr = JuMP.@expression(pm.model, sum(pg[i] for i in gens_in_area))
                delta_p_area[area] = JuMP.@expression(pm.model, abs(P_gen_area_expr - P_load_area[area]))
                H_min_area[area] = (delta_p_area[area] * f0) / (P_load_area[area] * 2 * rocof)
                H_area[area] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0) / sum(gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0)
                # Apply area-specific inertia constraints
                JuMP.@NLconstraint(pm.model, H_area[area] >= H_min_area[area])
            end
        end

        if bus == "true"
            # Bus-specific inertia constraints
            println("Adding minimum bus inertia constraint to ACPModel")

            # Initialize dictionaries to store bus-specific data
            H_bus = Dict()
            H_min_bus = Dict()
            P_load_bus = Dict()
            delta_p_bus = Dict()
            P_gen_bus_exprs = Dict()
            pg = var(pm, :pg)

            for j in keys(bus_data)
                P_gen_bus_exprs[j] = JuMP.@expression(pm.model, sum(pg[i] for (i, gen) in ref(pm, :gen) if gen["gen_bus"] == j))
                P_load_bus[j] = sum(load_data[i]["pd"] for i in eachindex(load_data) if load_data[i]["load_bus"] == j; init=0)
                delta_p_bus[j] = JuMP.@expression(pm.model, abs(P_gen_bus_exprs[j] - P_load_bus[j]))
                H_min_bus[j] = (delta_p_bus[j] * f0) / (P_load_bus[j] * 2 * rocof)
                println("H_min_bus $j = ", H_min_bus[j])
                H_bus[j] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * z[i] for i in eachindex(gen_data) if gen_data[i]["gen_bus"] == j; init=0) / sum(gen_data[i]["pmax"] * z[i] for i in eachindex(gen_data) if gen_data[i]["gen_bus"] == j; init=0)
                println("H_bus $j = ", H_bus[j])
                # Apply bus-specific inertia constraints
                JuMP.@NLconstraint(pm.model, H_bus[j] >= H_min_bus[j])
            end
        end
    end
end




function add_reactive_power_constraints(pm::AbstractPowerModel,v_options::Dict{String, String})
    gen_data = ref(pm, :gen)
    reactive_power_limit = v_options["reactive_power_limit"]
    max_rvc = v_options["max_rvc"]
    max_pas = v_options["max_pas"]

    if reactive_power_limit == "true"
        println("Adding reactive power limit constraints to ACPModel")
        for (gen_id, gen) in gen_data
            qg = var(pm, :qg, gen_id)
            q_min = gen["qmin"]
            q_max = gen["qmax"]
            JuMP.@constraint(pm.model, q_min <= qg <= q_max)
        end
    end
    # Add maximum rapid voltage change constraints if max_rvc is number
    if max_rvc == "true"
        println("Adding maximum rapid voltage change constraints to ACPModel")
        
    end

    # Add maximum phase angle shift constraints if max_pas is number
    if max_pas == "true"
        println("Adding maximum phase angle shift constraints to ACPModel")
        
    end
   
end

################################### End Taiseer Code #########################


"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error(_LOGGER, "$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end


# Generic thermal limit constraint
"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_from(pm::AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2)
end

"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_to(pm::AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2)
end

"`[rate_a, p[f_idx], q[f_idx]] in SecondOrderCone`"
function constraint_thermal_limit_from(pm::AbstractConicModels, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    JuMP.@constraint(pm.model, [rate_a, p_fr, q_fr] in JuMP.SecondOrderCone())
end

"`[rate_a, p[t_idx], q[t_idx]] in SecondOrderCone`"
function constraint_thermal_limit_to(pm::AbstractConicModels, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    JuMP.@constraint(pm.model, [rate_a, p_to, q_to] in JuMP.SecondOrderCone())
end

# Generic on/off thermal limit constraint

"`p[f_idx]^2 + q[f_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_from_on_off(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2*z^2)
end

"`p[t_idx]^2 + q[t_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_to_on_off(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2*z^2)
end

"`p_ne[f_idx]^2 + q_ne[f_idx]^2 <= (rate_a * branch_ne[i])^2`"
function constraint_ne_thermal_limit_from(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    p_fr = var(pm, n, :p_ne, f_idx)
    q_fr = var(pm, n, :q_ne, f_idx)
    z = var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2*z^2)
end

"`p_ne[t_idx]^2 + q_ne[t_idx]^2 <= (rate_a * branch_ne[i])^2`"
function constraint_ne_thermal_limit_to(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    p_to = var(pm, n, :p_ne, t_idx)
    q_to = var(pm, n, :q_ne, t_idx)
    z = var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2*z^2)
end

"`pg[i] == pg`"
function constraint_gen_setpoint_active(pm::AbstractPowerModel, n::Int, i, pg)
    pg_var = var(pm, n, :pg, i)

    JuMP.@constraint(pm.model, pg_var == pg)
end

"`qq[i] == qq`"
function constraint_gen_setpoint_reactive(pm::AbstractPowerModel, n::Int, i, qg)
    qg_var = var(pm, n, :qg, i)

    JuMP.@constraint(pm.model, qg_var == qg)
end

"on/off constraint for generators"
function constraint_gen_power_on_off(pm::AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)
    qg = var(pm, n, :qg, i)
    z = var(pm, n, :z_gen, i)

    JuMP.@constraint(pm.model, pg <= pmax*z)
    JuMP.@constraint(pm.model, pg >= pmin*z)
    JuMP.@constraint(pm.model, qg <= qmax*z)
    JuMP.@constraint(pm.model, qg >= qmin*z)
end


"""
Creates Line Flow constraint for DC Lines (Matpower Formulation)

```
p_fr + p_to == loss0 + p_fr * loss1
```
"""
function constraint_dcline_power_losses(pm::AbstractPowerModel, n::Int, f_bus, t_bus, f_idx, t_idx, loss0, loss1)
    p_fr = var(pm, n, :p_dc, f_idx)
    p_to = var(pm, n, :p_dc, t_idx)

    JuMP.@constraint(pm.model, (1-loss1) * p_fr + (p_to - loss0) == 0)
end

"`pf[i] == pf, pt[i] == pt`"
function constraint_dcline_setpoint_active(pm::AbstractPowerModel, n::Int, f_idx, t_idx, pf, pt)
    p_fr = var(pm, n, :p_dc, f_idx)
    p_to = var(pm, n, :p_dc, t_idx)

    JuMP.@constraint(pm.model, p_fr == pf)
    JuMP.@constraint(pm.model, p_to == pt)
end


"""
do nothing, most models to not require any model-specific voltage constraints
"""
function constraint_model_voltage(pm::AbstractPowerModel, n::Int)
end

"""
do nothing, most models to not require any model-specific on/off voltage constraints
"""
function constraint_model_voltage_on_off(pm::AbstractPowerModel, n::Int)
end

"""
do nothing, most models to not require any model-specific network expansion voltage constraints
"""
function constraint_ne_model_voltage(pm::AbstractPowerModel, n::Int)
end

"""
do nothing, most models to not require any model-specific current constraints
"""
function constraint_model_current(pm::AbstractPowerModel, n::Int)
end


""
function constraint_switch_state_open(pm::AbstractPowerModel, n::Int, f_idx)
    psw = var(pm, n, :psw, f_idx)
    qsw = var(pm, n, :qsw, f_idx)

    JuMP.@constraint(pm.model, psw == 0.0)
    JuMP.@constraint(pm.model, qsw == 0.0)
end

""
function constraint_switch_thermal_limit(pm::AbstractPowerModel, n::Int, f_idx, rating)
    psw = var(pm, n, :psw, f_idx)
    qsw = var(pm, n, :qsw, f_idx)

    JuMP.@constraint(pm.model, psw^2 + qsw^2 <= rating^2)
end

""
function constraint_switch_power_on_off(pm::AbstractPowerModel, n::Int, i, f_idx)
    psw = var(pm, n, :psw, f_idx)
    qsw = var(pm, n, :qsw, f_idx)
    z = var(pm, n, :z_switch, i)

    psw_lb, psw_ub = _IM.variable_domain(psw)
    qsw_lb, qsw_ub = _IM.variable_domain(qsw)

    JuMP.@constraint(pm.model, psw <= psw_ub*z)
    JuMP.@constraint(pm.model, psw >= psw_lb*z)
    JuMP.@constraint(pm.model, qsw <= qsw_ub*z)
    JuMP.@constraint(pm.model, qsw >= qsw_lb*z)
end



""
function constraint_storage_thermal_limit(pm::AbstractPowerModel, n::Int, i, rating)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)

    JuMP.@constraint(pm.model, ps^2 + qs^2 <= rating^2)
end

""
function constraint_storage_state_initial(pm::AbstractPowerModel, n::Int, i::Int, energy, charge_eff, discharge_eff, time_elapsed)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    se = var(pm, n, :se, i)

    JuMP.@constraint(pm.model, se - energy == time_elapsed*(charge_eff*sc - sd/discharge_eff))
end

""
function constraint_storage_state(pm::AbstractPowerModel, n_1::Int, n_2::Int, i::Int, charge_eff, discharge_eff, time_elapsed)
    sc_2 = var(pm, n_2, :sc, i)
    sd_2 = var(pm, n_2, :sd, i)
    se_2 = var(pm, n_2, :se, i)
    se_1 = var(pm, n_1, :se, i)

    JuMP.@constraint(pm.model, se_2 - se_1 == time_elapsed*(charge_eff*sc_2 - sd_2/discharge_eff))
end

""
function constraint_storage_complementarity_nl(pm::AbstractPowerModel, n::Int, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model, sc*sd == 0.0)
end

""
function constraint_storage_complementarity_mi(pm::AbstractPowerModel, n::Int, i, charge_ub, discharge_ub)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    sc_on = var(pm, n, :sc_on, i)
    sd_on = var(pm, n, :sd_on, i)

    JuMP.@constraint(pm.model, sc_on + sd_on == 1)
    JuMP.@constraint(pm.model, sc_on*charge_ub >= sc)
    JuMP.@constraint(pm.model, sd_on*discharge_ub >= sd)
end


""
function constraint_storage_on_off(pm::AbstractPowerModel, n::Int, i, pmin, pmax, qmin, qmax, charge_ub, discharge_ub)
    z_storage = var(pm, n, :z_storage, i)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@constraint(pm.model, ps <= z_storage*pmax)
    JuMP.@constraint(pm.model, ps >= z_storage*pmin)
    JuMP.@constraint(pm.model, qs <= z_storage*qmax)
    JuMP.@constraint(pm.model, qs >= z_storage*qmin)
    JuMP.@constraint(pm.model, qsc <= z_storage*qmax)
    JuMP.@constraint(pm.model, qsc >= z_storage*qmin)
end
