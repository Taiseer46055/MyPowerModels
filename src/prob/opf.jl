""

##################################################### start ##########################################

# This function is a custom version of solve_ac_opf with additional parameters for inertia constraints.
function solve_ac_opf_with_inertia(file, model_type::Type, optimizer, options::Dict{String, Dict{String}}; kwargs...)
    # Check that options has the keys "f" and "v"
    if !haskey(options, "f")
        error("Options must contain the keys 'f'")
    end

    # Check that "f" contains all required keys
    required_keys = ["inertia_constraint", "disturbance", "system", "area", "bus", "calc_delta_P", "alpha_factor", "rocof"]
    if !all(haskey(options["f"], key) for key in required_keys)
        error("Options 'f' must contain the keys: ", join(required_keys, ", "))
    end

    # Check that the values of "inertia_constraint", "system", "area", and "bus" are either "true" or "false"
    for key in ["inertia_constraint", "system", "area", "bus"]
        if options["f"][key] != "true" && options["f"][key] != "false"
            error("Option '$key' must be either 'true' or 'false'")
        end
    end

    # Check that the value of "disturbance" is either "small" or "large"
    if options["f"]["disturbance"] != "small" && options["f"]["disturbance"] != "large"
        error("Option 'disturbance' must be either 'small' or 'large'")
    end

    # Check that the value of "calc_delta_P" is either "internal" or an array of two numbers
    if !(options["f"]["calc_delta_P"] == "internal" || (isa(options["f"]["calc_delta_P"], Array) && length(options["f"]["calc_delta_P"]) == 2))
        error("Option 'calc_delta_P' must be either 'internal' or an array of two numbers")
    end

    # Check that the value of "alpha_factor" is a number between 0 and 1
    if !(isa(options["f"]["alpha_factor"], Number) && 0 <= options["f"]["alpha_factor"] <= 1)
        error("Option 'alpha_factor' must be a number between 0 and 1")
    end

    # Check that the value of "rocof" is a positive number
    if !(isa(options["f"]["rocof"], Number) && options["f"]["rocof"] > 0)
        error("Option 'rocof' must be a positive number")
    end

    # Check that the value of "weighted_area" is either "load", "equal", or "none"

    # Call a custom solve_opf function with additional parameters.

    return solve_opf_inertia(file, model_type, optimizer, options; kwargs...)

end

# Custom solve_opf function that includes inertia constraints.

function solve_opf_inertia(file, model_type, optimizer, options; kwargs...)
    model_builder = build_opf_H_min(model_type, options)
    return solve_model(file, model_type, optimizer, model_builder; kwargs...)
end

# Custom build function for OPF that includes standard constraints and additional inertia constraints.
function build_opf_H_min(model_type::Type, options::Dict{String, Dict{String}})
    # Define a function to build the OPF model.
    function build_my_opf(pm::AbstractPowerModel)
        # Extract data from the power model
        
        gen_data = ref(pm, :gen)
        load_data = ref(pm, :load)
        f_options = options["f"]
        v_options = options["v"]
        alpha = f_options["alpha_factor"]
        calc_delta_P = f_options["calc_delta_P"]
        baseMVA = ref(pm, :baseMVA)
            # Process delta_P based on the calculation method specified in options
        if isa(calc_delta_P, Array) && length(calc_delta_P) == 2
            gen_id = calc_delta_P[1]
            delta_P = calc_delta_P[2]/baseMVA
            println("delta_P: ", delta_P)

            gen_data[gen_id]["pmax"] = max(0, gen_data[gen_id]["pmax"] - delta_P)
            gen_data[gen_id]["pmin"] = max(0, gen_data[gen_id]["pmin"] - delta_P)

            println("Pmax of generator $gen_id: ", gen_data[gen_id]["pmax"])
            println("Pmin of generator $gen_id: ", gen_data[gen_id]["pmin"])

        elseif isa(calc_delta_P, String) && calc_delta_P == "internal"
            # Calculate delta_P internally based on generator data
            H_pmax_product = Dict(i => gen_data[i]["H"] * gen_data[i]["pmax"] for i in eachindex(gen_data))
            println("H_pmax_product: ", H_pmax_product)
            max_product_gen = argmax(H_pmax_product)[1]
            println("Generator with maximum product of H and Pmax: ", max_product_gen)
            println("Pmax of generator with maximum product: ", gen_data[max_product_gen]["pmax"])
            println("H of generator with maximum product: ", gen_data[max_product_gen]["H"])

            delta_P = (gen_data[max_product_gen]["pmax"])*alpha
            println("delta_P: ", delta_P)

            gen_data[max_product_gen]["pmax"] = max(0, gen_data[max_product_gen]["pmax"]*(1-alpha))
            gen_data[max_product_gen]["pmin"] = max(0, gen_data[max_product_gen]["pmin"]*(1-alpha))

        else
            error("Invalid value for calc_delta_P. It must be either an array of number or 'internal'.")
        end

        # calc H_min
        rocof = f_options["rocof"]
        f0 = 50.0
        P_load = sum(haskey(load, "pd") ? load["pd"] : 0.0 for (_, load) in load_data)
        H_min = (delta_P * f0) / (P_load * 2 * rocof)

        variable_bus_voltage(pm)

        variable_gen_indicator(pm)
        variable_gen_power_on_off(pm)

        variable_storage_indicator(pm)
        variable_storage_power_mi_on_off(pm)

        variable_branch_power(pm)
        variable_dcline_power(pm)

        objective_min_fuel_and_flow_cost(pm)

        constraint_model_voltage(pm)

        for i in ids(pm, :ref_buses)
            constraint_theta_ref(pm, i)
        end
        for i in ids(pm, :gen)
            constraint_gen_power_on_off(pm, i)
        end
        for i in ids(pm, :storage)
            constraint_storage_on_off(pm, i)
        end
        for i in ids(pm, :bus)
            constraint_power_balance(pm, i)
        end
        for i in ids(pm, :storage)
            constraint_storage_state(pm, i)
            constraint_storage_complementarity_mi(pm, i)
            constraint_storage_losses(pm, i)
            constraint_storage_thermal_limit(pm, i)
        end
        for i in ids(pm, :branch)
            constraint_ohms_yt_from(pm, i)
            constraint_ohms_yt_to(pm, i)
    
            constraint_voltage_angle_difference(pm, i)
            constraint_thermal_limit_from(pm, i)
            constraint_thermal_limit_to(pm, i)
        end

        for i in ids(pm, :dcline)
            constraint_dcline_power_losses(pm, i)
        end
        # Inertia constraints
        n = nw_id_default
        if haskey(options, "f")
            f_options = options["f"]
            inertia_constraint = f_options["inertia_constraint"]
            if inertia_constraint == "true"
                println("Add inertia constraints to the model.")
                if model_type == DCPPowerModel
                    constraint_system_inertia(pm, H_min, f_options, n)
                elseif model_type == ACPPowerModel
                    constraint_system_inertia(pm, H_min, f_options, n)
                end
            end
        else
            println("The required variables are not entered. Please check your options.")
        end

        if(haskey(options, "v"))
            v_options = options["v"]
            voltage_constraint = v_options["voltage_constraint"]
            if voltage_constraint == "true"
                println("Add voltage constraints to the model.")
                if model_type == DCPPowerModel
                    constriant_reactive_power(pm, v_options, n)
                elseif model_type == ACPPowerModel
                    constriant_reactive_power(pm, v_options, n)
                end
            end
        else
            println("The required variables are not entered. Please check your options.")
        end
        
        # add H_min and delta_P as attributes to the data model

        pm.ext[:H_min] = H_min
        pm.ext[:delta_P] = delta_P

        pm.data["H_min"] = pm.ext[:H_min]
        pm.data["delta_P"] = pm.ext[:delta_P]
    end
    return build_my_opf
end

function solve_mn_ac_opf_with_inertia(file, model_type::Type, optimizer, options::Dict{String, Dict{String}}; kwargs...)

    # Check that options has the keys "f" and "v"
    if !haskey(options, "f")
    error("Options must contain the keys 'f'")
    end

    # Check that "f" contains all required keys
    required_keys = ["inertia_constraint", "disturbance", "system", "area", "bus", "calc_delta_P", "alpha_factor", "rocof"]
    if !all(haskey(options["f"], key) for key in required_keys)
        error("Options 'f' must contain the keys: ", join(required_keys, ", "))
    end

    # Check that the values of "inertia_constraint", "system", "area", and "bus" are either "true" or "false"
    for key in ["inertia_constraint", "system", "area", "bus"]
        if options["f"][key] != "true" && options["f"][key] != "false"
            error("Option '$key' must be either 'true' or 'false'")
        end
    end

    # Check that the value of "disturbance" is either "small" or "large"
    if options["f"]["disturbance"] != "small" && options["f"]["disturbance"] != "large"
        error("Option 'disturbance' must be either 'small' or 'large'")
    end

    # Check that the value of "calc_delta_P" is either "internal" or an array of two numbers
    if !(options["f"]["calc_delta_P"] == "internal" || (isa(options["f"]["calc_delta_P"], Array) && length(options["f"]["calc_delta_P"]) == 2))
        error("Option 'calc_delta_P' must be either 'internal' or an array of two numbers")
    end

    # Check that the value of "alpha_factor" is a number between 0 and 1
    if !(isa(options["f"]["alpha_factor"], Number) && 0 <= options["f"]["alpha_factor"] <= 1)
        error("Option 'alpha_factor' must be a number between 0 and 1")
    end

    # Check that the value of "rocof" is a positive number
    if !(isa(options["f"]["rocof"], Number) && options["f"]["rocof"] > 0)
        error("Option 'rocof' must be a positive number")
    end


    return solve_mn_opf_inertia(file, model_type, optimizer, options; kwargs...)
end

function solve_mn_opf_inertia(file, model_type, optimizer, options; kwargs...)

    model_builder = build_mn_opf_inertia(model_type, options)
    return solve_model(file, model_type, optimizer, model_builder ;ref_extensions=[ref_add_connected_components!], multinetwork=true, kwargs...)
end

function build_mn_opf_inertia(model_type::Type, options::Dict{String, Dict{String}})
    # Define a function to build the OPF model for each network.
    function build_my_mn_opf_inertia(pm::AbstractPowerModel)
        # Iterate through each network
        for (n, network) in nws(pm)
            # Extract data from the current network in the power model
            gen_data = ref(pm, n, :gen)
            load_data = ref(pm, n, :load)
            baseMVA = ref(pm, n, :baseMVA)
            f_options = options["f"]
            v_options = options["v"]
            alpha = f_options["alpha_factor"]
            calc_delta_P = f_options["calc_delta_P"]            

            if isa(calc_delta_P, Array) && length(calc_delta_P) == 2
                gen_id = calc_delta_P[1]
                delta_P = calc_delta_P[2]/baseMVA
                println("delta_P: ", delta_P)
    
                gen_data[gen_id]["pmax"] = max(0, gen_data[gen_id]["pmax"] - delta_P)
                gen_data[gen_id]["pmin"] = max(0, gen_data[gen_id]["pmin"] - delta_P)
    
                println("Pmax of generator $gen_id: ", gen_data[gen_id]["pmax"])
                println("Pmin of generator $gen_id: ", gen_data[gen_id]["pmin"])
    
            elseif isa(calc_delta_P, String) && calc_delta_P == "internal"
                # Calculate delta_P internally based on generator data
                H_pmax_product = Dict(i => gen_data[i]["H"] * gen_data[i]["pmax"] for i in eachindex(gen_data))
                println("H_pmax_product: ", H_pmax_product)
                max_product_gen = argmax(H_pmax_product)[1]
                println("Generator with maximum product of H and Pmax: ", max_product_gen)
                println("Pmax of generator with maximum product: ", gen_data[max_product_gen]["pmax"])
                println("H of generator with maximum product: ", gen_data[max_product_gen]["H"])
    
                delta_P = (gen_data[max_product_gen]["pmax"])*alpha
                println("delta_P: ", delta_P)
    
                gen_data[max_product_gen]["pmax"] = max(0, gen_data[max_product_gen]["pmax"]*(1-alpha))
                gen_data[max_product_gen]["pmin"] = max(0, gen_data[max_product_gen]["pmin"]*(1-alpha))
    
            else
                error("Invalid value for calc_delta_P. It must be either an array of number or 'internal'.")
            end
    
            # calc H_min
            rocof = f_options["rocof"]
            f0 = 50.0
            P_load = sum(haskey(load, "pd") ? load["pd"] : 0.0 for (_, load) in load_data)
            H_min = (delta_P * f0) / (P_load * 2 * rocof)
    
            variable_bus_voltage(pm, nw=n)
    
            variable_gen_indicator(pm, nw=n)
            variable_gen_power_on_off(pm, nw=n)
    
            variable_storage_indicator(pm, nw=n)
            variable_storage_power_mi_on_off(pm, nw=n)
    
            variable_branch_power(pm, nw=n)
            variable_dcline_power(pm, nw=n)
        
            constraint_model_voltage(pm, nw=n)


            for i in ids(pm, :ref_buses, nw=n)
                constraint_theta_ref(pm, i, nw=n)
            end
            for i in ids(pm, :gen, nw=n)
                constraint_gen_power_on_off(pm, i, nw=n)
            end
            for i in ids(pm, :storage, nw=n)
                constraint_storage_on_off(pm, i, nw=n)
            end
            for i in ids(pm, :bus, nw=n)
                constraint_power_balance(pm, i, nw=n)
            end
            for i in ids(pm, :storage, nw=n)
                constraint_storage_state(pm, i, nw=n)
                constraint_storage_complementarity_mi(pm, i, nw=n)
                constraint_storage_losses(pm, i, nw=n)
                constraint_storage_thermal_limit(pm, i, nw=n)
            end
            for i in ids(pm, :branch, nw=n)
                constraint_ohms_yt_from(pm, i, nw=n)
                constraint_ohms_yt_to(pm, i, nw=n)
        
                constraint_voltage_angle_difference(pm, i, nw=n)
                constraint_thermal_limit_from(pm, i, nw=n)
                constraint_thermal_limit_to(pm, i, nw=n)
            end
    
            for i in ids(pm, :dcline, nw=n)
                constraint_dcline_power_losses(pm, i, nw=n)
            end
            # Inertia constraints
            if haskey(options, "f")
                f_options = options["f"]
                inertia_constraint = f_options["inertia_constraint"]
                if inertia_constraint == "true"
                    println("Add inertia constraints to the model.")
                    if model_type == DCPPowerModel
                        constraint_system_inertia(pm, H_min, f_options, n)
                    elseif model_type == ACPPowerModel
                        constraint_system_inertia(pm, H_min, f_options, n)
                    end
                end
            else
                println("The required variables are not entered. Please check your options.")
            end
    
            if(haskey(options, "v"))
                v_options = options["v"]
                voltage_constraint = v_options["voltage_constraint"]
                if voltage_constraint == "true"
                    println("Add voltage constraints to the model.")
                    if model_type == DCPPowerModel
                        constriant_reactive_power(pm, v_options, n)
                    elseif model_type == ACPPowerModel
                        constriant_reactive_power(pm, v_options, n)
                    end
                end
            else
                println("The required variables are not entered. Please check your options.")
            end
            
            # add H_min and delta_P as attributes to the data model
    
            pm.ext[:H_min] = H_min
            pm.ext[:delta_P] = delta_P
    
            pm.data["H_min"] = pm.ext[:H_min]
            pm.data["delta_P"] = pm.ext[:delta_P]
        end
        objective_min_fuel_and_flow_cost(pm)

    end
    return build_my_mn_opf_inertia
end






##################################################### End ############################################




function solve_ac_opf(file, optimizer; kwargs...)
    return solve_opf(file, ACPPowerModel, optimizer; kwargs...)
end

""
function solve_dc_opf(file, optimizer; kwargs...)
    return solve_opf(file, DCPPowerModel, optimizer; kwargs...)
end

""
function solve_opf(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_opf; kwargs...)
end

""
function build_opf(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end

end



"a toy example of how to model with multi-networks"
function solve_mn_opf(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_mn_opf; multinetwork=true, kwargs...)
end

""
function build_mn_opf(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n)
        variable_branch_power(pm, nw=n)
        variable_dcline_power(pm, nw=n)

        constraint_model_voltage(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, :branch, nw=n)
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)

            constraint_voltage_angle_difference(pm, i, nw=n)

            constraint_thermal_limit_from(pm, i, nw=n)
            constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in ids(pm, :dcline, nw=n)
            constraint_dcline_power_losses(pm, i, nw=n)
        end
    end

    objective_min_fuel_and_flow_cost(pm)
end


"a toy example of how to model with multi-networks and storage"
function solve_mn_opf_strg(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_mn_opf_strg; multinetwork=true, kwargs...)
end

""
function build_mn_opf_strg(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n)
        variable_storage_power_mi(pm, nw=n)
        variable_branch_power(pm, nw=n)
        variable_dcline_power(pm, nw=n)

        constraint_model_voltage(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, :storage, nw=n)
            constraint_storage_complementarity_mi(pm, i, nw=n)
            constraint_storage_losses(pm, i, nw=n)
            constraint_storage_thermal_limit(pm, i, nw=n)
        end

        for i in ids(pm, :branch, nw=n)
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)

            constraint_voltage_angle_difference(pm, i, nw=n)

            constraint_thermal_limit_from(pm, i, nw=n)
            constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in ids(pm, :dcline, nw=n)
            constraint_dcline_power_losses(pm, i, nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    objective_min_fuel_and_flow_cost(pm)
end




"""
Solves an opf using ptdfs with no explicit voltage or line flow variables.

This formulation is most often used when a small subset of the line flow
constraints are active in the data model.
"""
function solve_opf_ptdf(file, model_type::Type, optimizer; full_inverse=false, kwargs...)
    if !full_inverse
        return solve_model(file, model_type, optimizer, build_opf_ptdf; ref_extensions=[ref_add_connected_components!,ref_add_sm!], kwargs...)
    else
        return solve_model(file, model_type, optimizer, build_opf_ptdf; ref_extensions=[ref_add_connected_components!,ref_add_sm_inv!], kwargs...)
    end
end

""
function build_opf_ptdf(pm::AbstractPowerModel)
    Memento.error(_LOGGER, "build_opf_ptdf is only valid for DCPPowerModels")
end

""
function build_opf_ptdf(pm::DCPPowerModel)
    variable_gen_power(pm)

    for i in ids(pm, :bus)
        expression_bus_power_injection(pm, i)
    end

    objective_min_fuel_cost(pm)

    constraint_model_voltage(pm)

    # this constraint is implicit in this model
    #for i in ids(pm, :ref_buses)
    #    constraint_theta_ref(pm, i)
    #end

    for i in ids(pm, :components)
        constraint_network_power_balance(pm, i)
    end

    for (i, branch) in ref(pm, :branch)
        # requires optional vad parameters
        #constraint_voltage_angle_difference(pm, i)

        # only create these expressions if a line flow is specified
        if haskey(branch, "rate_a")
            expression_branch_power_ohms_yt_from_ptdf(pm, i)
            expression_branch_power_ohms_yt_to_ptdf(pm, i)
        end

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end
end


""
function ref_add_sm!(ref::Dict{Symbol, <:Any}, data::Dict{String, <:Any})
    apply_pm!(_ref_add_sm!, ref, data)
end


""
function _ref_add_sm!(ref::Dict{Symbol, <:Any}, data::Dict{String, <:Any})
    reference_bus(data) # throws an error if an incorrect number of reference buses are defined
    ref[:sm] = calc_susceptance_matrix(data)
end


""
function ref_add_sm_inv!(ref::Dict{Symbol, <:Any}, data::Dict{String, <:Any})
    apply_pm!(_ref_add_sm_inv!, ref, data)
end


""
function _ref_add_sm_inv!(ref::Dict{Symbol, <:Any}, data::Dict{String, <:Any})
    ref[:sm] = calc_susceptance_matrix_inv(data)
end
