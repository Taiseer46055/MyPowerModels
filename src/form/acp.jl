### polar form of the non-convex AC equations

################################### Start Taiseer Code #########################

# Defines a function to calculate system inertia constraints in an AC Power Model

function constraint_system_inertia(pm::AbstractACPModel, H_min::Float64 , f_options::Dict{String, Any})

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

################################### End Taiseer Code #########################


""
function variable_bus_voltage(pm::AbstractACPModel; kwargs...)
    variable_bus_voltage_angle(pm; kwargs...)
    variable_bus_voltage_magnitude(pm; kwargs...)
end

""
function sol_data_model!(pm::AbstractACPModel, solution::Dict)
    # nothing to do, this is in the data model space by default
end


""
function variable_ne_branch_voltage(pm::AbstractACPModel; kwargs...)
end

"do nothing, this model does not have complex voltage constraints"
function constraint_ne_model_voltage(pm::AbstractACPModel; kwargs...)
end


"`v[i] == vm`"
function constraint_voltage_magnitude_setpoint(pm::AbstractACPModel, n::Int, i::Int, vm)
    v = var(pm, n, :vm, i)
    JuMP.@constraint(pm.model, v == vm)
end


function constraint_current_limit_from(pm::AbstractACPModel, n::Int, f_idx, c_rating_a)
    l,i,j = f_idx

    vm_fr = var(pm, n, :vm, i)

    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= vm_fr^2*c_rating_a^2)
end

function constraint_current_limit_to(pm::AbstractACPModel, n::Int, t_idx, c_rating_a)
    l,j,i = t_idx

    vm_to = var(pm, n, :vm, j)

    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= vm_to^2*c_rating_a^2)
end


""
function constraint_power_balance(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = var(pm, n, :vm, i)
    p    = get(var(pm, n),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); _check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    # the check "typeof(p[arc]) <: JuMP.NonlinearExpression" is required for the
    # case when p/q are nonlinear expressions instead of decision variables
    # once NLExpressions are first order in JuMP it should be possible to
    # remove this.
    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    if !nl_form
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    else
        cstr_p = JuMP.@NLconstraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    end

    if !nl_form
        cstr_q = JuMP.@constraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    else
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

""
function constraint_power_balance_ls(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = var(pm, n, :vm, i)
    p    = get(var(pm, n),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); _check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    z_demand = get(var(pm, n), :z_demand, Dict()); _check_var_keys(z_demand, keys(bus_pd), "power factor", "load")
    z_shunt = get(var(pm, n), :z_shunt, Dict()); _check_var_keys(z_shunt, keys(bus_gs), "power factor", "shunt")

    # this is required for improved performance in NLP models
    if length(z_shunt) <= 0
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd*z_demand[i] for (i,pd) in bus_pd)
            - sum(gs*z_shunt[i] for (i,gs) in bus_gs)*vm^2
        )
        cstr_q = JuMP.@constraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd*z_demand[i] for (i,qd) in bus_qd)
            + sum(bs*z_shunt[i] for (i,bs) in bus_bs)*vm^2
        )
    else
        cstr_p = JuMP.@NLconstraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd*z_demand[i] for (i,pd) in bus_pd)
            - sum(gs*z_shunt[i] for (i,gs) in bus_gs)*vm^2
        )
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd*z_demand[i] for (i,qd) in bus_qd)
            + sum(bs*z_shunt[i] for (i,bs) in bus_bs)*vm^2
        )
    end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_ne_power_balance(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_arcs_ne, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = var(pm, n, :vm, i)
    p    = get(var(pm, n),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); _check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")
    p_ne = get(var(pm, n), :p_ne, Dict()); _check_var_keys(p_ne, bus_arcs_ne, "active power", "ne_branch")
    q_ne = get(var(pm, n), :q_ne, Dict()); _check_var_keys(q_ne, bus_arcs_ne, "reactive power", "ne_branch")

    JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(p_ne[a] for a in bus_arcs_ne)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*vm^2
    )
    JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(q_ne[a] for a in bus_arcs_ne)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*vm^2
    )
end



""
function expression_branch_power_ohms_yt_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    var(pm, n, :p)[f_idx] = JuMP.@NLexpression(pm.model,  (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    var(pm, n, :q)[f_idx] = JuMP.@NLexpression(pm.model, -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end

""
function expression_branch_power_ohms_yt_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    var(pm, n, :p)[t_idx] = JuMP.@NLexpression(pm.model,  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
    var(pm, n, :q)[t_idx] = JuMP.@NLexpression(pm.model, -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_yt_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n,  :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_fr ==  (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_yt_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_to  = var(pm, n,  :p, t_idx)
    q_to  = var(pm, n,  :q, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
    JuMP.@NLconstraint(pm.model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
end

"""
Creates Ohms constraints for AC models (y post fix indicates that Y values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)*(v[f_bus]/tr)^2 + -g*v[f_bus]/tr*v[t_bus]*cos(t[f_bus]-t[t_bus]-as) + -b*v[f_bus]/tr*v[t_bus]*sin(t[f_bus]-t[t_bus]-as)
q[f_idx] == -(b+b_fr)*(v[f_bus]/tr)^2 + b*v[f_bus]/tr*v[t_bus]*cos(t[f_bus]-t[t_bus]-as) + -g*v[f_bus]/tr*v[t_bus]*sin(t[f_bus]-t[t_bus]-as)
```
"""
function constraint_ohms_y_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tm, ta)
    p_fr  = var(pm, n,  :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_fr ==  (g+g_fr)*(vm_fr/tm)^2 - g*vm_fr/tm*vm_to*cos(va_fr-va_to-ta) + -b*vm_fr/tm*vm_to*sin(va_fr-va_to-ta) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b+b_fr)*(vm_fr/tm)^2 + b*vm_fr/tm*vm_to*cos(va_fr-va_to-ta) + -g*vm_fr/tm*vm_to*sin(va_fr-va_to-ta) )
end

"""
Creates Ohms constraints for AC models (y post fix indicates that Y values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + -g*v[t_bus]*v[f_bus]/tr*cos(t[t_bus]-t[f_bus]+as) + -b*v[t_bus]*v[f_bus]/tr*sin(t[t_bus]-t[f_bus]+as)
q[t_idx] == -(b+b_to)*v[t_bus]^2 +  b*v[t_bus]*v[f_bus]/tr*cos(t[f_bus]-t[t_bus]+as) + -g*v[t_bus]*v[f_bus]/tr*sin(t[t_bus]-t[f_bus]+as)
```
"""
function constraint_ohms_y_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tm, ta)
    p_to  = var(pm, n,  :p, t_idx)
    q_to  = var(pm, n,  :q, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_to ==  (g+g_to)*vm_to^2 - g*vm_to*vm_fr/tm*cos(va_to-va_fr+ta) + -b*vm_to*vm_fr/tm*sin(va_to-va_fr+ta) )
    JuMP.@NLconstraint(pm.model, q_to == -(b+b_to)*vm_to^2 + b*vm_to*vm_fr/tm*cos(va_to-va_fr+ta) + -g*vm_to*vm_fr/tm*sin(va_to-va_fr+ta) )
end


""
function constraint_switch_state_closed(pm::AbstractACPModel, n::Int, f_bus, t_bus)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@constraint(pm.model, vm_fr == vm_to)
    JuMP.@constraint(pm.model, va_fr == va_to)
end

""
function constraint_switch_voltage_on_off(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, vad_min, vad_max)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_switch, i)

    JuMP.@constraint(pm.model, z*vm_fr == z*vm_to)
    JuMP.@constraint(pm.model, z*va_fr == z*va_to)
end

""
function variable_bus_voltage_on_off(pm::AbstractACPModel; kwargs...)
    variable_bus_voltage_angle(pm; kwargs...)
    variable_bus_voltage_magnitude(pm; kwargs...)
end

"do nothing, this model does not have complex voltage constraints"
function constraint_model_voltage_on_off(pm::AbstractACPModel; kwargs...)
end

"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_yt_from_on_off(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = var(pm, n,  :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_branch, i)

    JuMP.@NLconstraint(pm.model, p_fr == z*( (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
    JuMP.@NLconstraint(pm.model, q_fr == z*(-(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_yt_to_on_off(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_to  = var(pm, n,  :p, t_idx)
    q_to  = var(pm, n,  :q, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_branch, i)

    JuMP.@NLconstraint(pm.model, p_to == z*( (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
    JuMP.@NLconstraint(pm.model, q_to == z*(-(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
end

"""
```
p_ne[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q_ne[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ne_ohms_yt_from(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = var(pm, n, :p_ne, f_idx)
    q_fr  = var(pm, n, :q_ne, f_idx)
    vm_fr = var(pm, n,   :vm, f_bus)
    vm_to = var(pm, n,   :vm, t_bus)
    va_fr = var(pm, n,   :va, f_bus)
    va_to = var(pm, n,   :va, t_bus)
    z = var(pm, n, :branch_ne, i)

    JuMP.@NLconstraint(pm.model, p_fr == z*( (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
    JuMP.@NLconstraint(pm.model, q_fr == z*(-(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
end

"""
```
p_ne[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q_ne[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ne_ohms_yt_to(pm::AbstractACPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_to = var(pm, n, :p_ne, t_idx)
    q_to = var(pm, n, :q_ne, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :branch_ne, i)

    JuMP.@NLconstraint(pm.model, p_to == z*( (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
    JuMP.@NLconstraint(pm.model, q_to == z*(-(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
end

"""
Creates Ohms constraints with variables for complex transformation ratio (y post fix indicates  Y is in rectangular form)
```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus]-ta)) + (-b)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]-ta))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus]-ta)) + (-g)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]-ta))
```
"""
function constraint_ohms_y_oltc_pst_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr)
    p_fr  = var(pm, n, :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    tm = var(pm, n, :tm, f_idx[1])
    ta = var(pm, n, :ta, f_idx[1])

    JuMP.@NLconstraint(pm.model, p_fr ==  (g+g_fr)/tm^2*vm_fr^2 + (-g)/tm*(vm_fr*vm_to*cos(va_fr-va_to-ta)) + (-b)/tm*(vm_fr*vm_to*sin(va_fr-va_to-ta)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b+b_fr)/tm^2*vm_fr^2 - (-b)/tm*(vm_fr*vm_to*cos(va_fr-va_to-ta)) + (-g)/tm*(vm_fr*vm_to*sin(va_fr-va_to-ta)) )
end

"""
Creates Ohms constraints with variables for complex transformation ratio (y post fix indicates  Y is in rectangular form)
```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus]+ta)) + (-b)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]+ta))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus]+ta)) + (-g)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]+ta))
```
"""
function constraint_ohms_y_oltc_pst_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to)
    p_to  = var(pm, n,  :p, t_idx)
    q_to  = var(pm, n,  :q, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    tm = var(pm, n, :tm, f_idx[1])
    ta = var(pm, n, :ta, f_idx[1])

    JuMP.@NLconstraint(pm.model, p_to ==  (g+g_to)*vm_to^2 + -g/tm*(vm_to*vm_fr*cos(va_to-va_fr+ta)) + -b/tm*(vm_to*vm_fr*sin(va_to-va_fr+ta)) )
    JuMP.@NLconstraint(pm.model, q_to == -(b+b_to)*vm_to^2 - -b/tm*(vm_to*vm_fr*cos(va_to-va_fr+ta)) + -g/tm*(vm_to*vm_fr*sin(va_to-va_fr+ta)) )
end


""
function constraint_ohms_y_pst_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tm)
    p_fr  = var(pm, n, :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    ta = var(pm, n, :ta, f_idx[1])

    JuMP.@NLconstraint(pm.model, p_fr ==  (g+g_fr)/tm^2*vm_fr^2 + (-g)/tm*(vm_fr*vm_to*cos(va_fr-va_to-ta)) + (-b)/tm*(vm_fr*vm_to*sin(va_fr-va_to-ta)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b+b_fr)/tm^2*vm_fr^2 - (-b)/tm*(vm_fr*vm_to*cos(va_fr-va_to-ta)) + (-g)/tm*(vm_fr*vm_to*sin(va_fr-va_to-ta)) )
end

""
function constraint_ohms_y_pst_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tm)
    p_to  = var(pm, n,  :p, t_idx)
    q_to  = var(pm, n,  :q, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    ta = var(pm, n, :ta, f_idx[1])

    JuMP.@NLconstraint(pm.model, p_to ==  (g+g_to)*vm_to^2 + -g/tm*(vm_to*vm_fr*cos(va_to-va_fr+ta)) + -b/tm*(vm_to*vm_fr*sin(va_to-va_fr+ta)) )
    JuMP.@NLconstraint(pm.model, q_to == -(b+b_to)*vm_to^2 - -b/tm*(vm_to*vm_fr*cos(va_to-va_fr+ta)) + -g/tm*(vm_to*vm_fr*sin(va_to-va_fr+ta)) )
end


"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax`"
function constraint_voltage_angle_difference_on_off(pm::AbstractACPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, z*(va_fr - va_to) <= z*angmax)
    JuMP.@constraint(pm.model, z*(va_fr - va_to) >= z*angmin)
end

"`angmin <= branch_ne[i]*(t[f_bus] - t[t_bus]) <= angmax`"
function constraint_ne_voltage_angle_difference(pm::AbstractACPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, angmin <= z*(va_fr - va_to) <= angmax)
end

"""
```
p[f_idx] + p[t_idx] >= 0
q[f_idx] + q[t_idx] >= -c/2*(v[f_bus]^2/tr^2 + v[t_bus]^2)
```
"""
function constraint_power_losses_lb(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g_fr, b_fr, g_to, b_to, tr)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    @assert(g_fr == 0 && g_to == 0)
    c = b_fr + b_to

    # TODO: Derive updated constraint from first principles
    JuMP.@constraint(m, p_fr + p_to >= 0)
    JuMP.@constraint(m, q_fr + q_to >= -c/2*(vm_fr^2/tr^2 + vm_to^2))
end


""
function constraint_storage_current_limit(pm::AbstractACPModel, n::Int, i, bus, rating)
    vm = var(pm, n, :vm, bus)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)

    JuMP.@constraint(pm.model, ps^2 + qs^2 <= rating^2*vm^2)
end


""
function constraint_storage_losses(pm::AbstractACPModel, n::Int, i, bus, r, x, p_loss, q_loss)
    vm = var(pm, n, :vm, bus)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@NLconstraint(pm.model, ps + (sd - sc) == p_loss + r*(ps^2 + qs^2)/vm^2)
    JuMP.@NLconstraint(pm.model, qs == qsc + q_loss + x*(ps^2 + qs^2)/vm^2)
end

