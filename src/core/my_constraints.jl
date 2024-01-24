
function constraint_min_system_inertia(pm::AbstractPowerModel, gen_id::Int, delta_P::Float64, max_rocof::Float64)
    println("Adding minimum system inertia constraint")
    
    # Retrieve generator and bus data
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)

    # Check if the generator id exists
    if !haskey(gen_data, gen_id)
        error("Generator id $gen_id does not exist in the network")
    end

    # Get the specific generator
    gen_at_bus = gen_data[gen_id]
    if haskey(gen_at_bus, "pg")
        gen_at_bus["pg"] -= delta_P
    else
        error("Key 'pg' not found in generator with ID $gen_id")
    end

    # Set the base frequency f0
    f0 = 50.0
    
    # Calculate P_LOAD as the sum of Pd for all buses with a default value of 0
    P_LOAD = 0.0
    for (_, bus) in bus_data
        if haskey(bus, "pd")
            P_LOAD += bus["pd"]
        end
    end

    # Calculate the minimum system inertia H_min
    H_min = (delta_P * f0) / (P_LOAD * 2 * max_rocof)
    
    H_sys = 0.0
    total_Pg = 0.0
    for (_, gen) in gen_data
        if haskey(gen, "pg") && gen["pg"] > 0
            H_sys += gen["H"] * gen["pg"]
            total_Pg += gen["pg"]
        end
    end
    H_sys = total_Pg > 0 ? H_sys / total_Pg : 0.0

    # Add the inertia constraint to the model
    JuMP.@constraint(pm.model, H_sys >= H_min)
end












#=

function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id::Int, gen_tech::Int, delta_P::Float64, max_rocof::Float64)
    println("Adding minimum system inertia constraint")
    
    # Retrieve generator and bus data
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)
#=
    # Debugging: Output the content of gen_data and bus_data
    println("Generator data: ", gen_data)
    println("Bus data: ", bus_data)
=#
    # Check if the bus number exists
    if !haskey(bus_data, bus_id)
        error("Bus number $bus_id does not exist in the network")
    end

    # Search for the specific generator
    gen_at_bus = nothing
    for (gen_id, gen) in gen_data
        println("Check generator $gen_id at the bus $(gen["gen_bus"]) with GenTech $(gen["GenTech"])")
        if string(gen["gen_bus"]) == string(bus_id) && string(gen["GenTech"]) == string(gen_tech)
            if haskey(gen, "pg")
                gen_at_bus = gen
                gen_at_bus["pg"] -= delta_P
                break
            else
                error("Key 'pg' not found in generator with ID $gen_id")
            end
        end
    end

    if gen_at_bus === nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end
    # Set the base frequency f0
    f0 = 50.0
    
    # Calculate P_LOAD as the sum of Pd for all buses with a default value of 0
    P_LOAD = 0.0
    for (_, bus) in bus_data
        if haskey(bus, "pd")
            P_LOAD += bus["pd"]
        end
    end


    # Calculate the minimum system inertia H_min
    H_min = (delta_P * f0) / (P_LOAD * 2 * max_rocof)

#=
    # Calculate the weighted sum of inertia for all generators with Pg > 0
    sum_H_gen = sum(gen["H"] * gen["pg"] for (_, gen) in gen_data if haskey(gen, "pg") && gen["pg"] > 0)
    sum_P_gen = sum(gen["pg"] for (_, gen) in gen_data if haskey(gen, "pg") && gen["pg"] > 0)

    # Calculate H_sys as the weighted average of inertia constants
    H_sys = sum_P_gen > 0 ? sum_H_gen / sum_P_gen : 0.0

=#
    # Calculate the minimum system inertia (H_min) and system inertia (H_sys)
    H_min = delta_P / max_rocof
    H_sys = 0.0
    total_Pg = 0.0
    for (_, gen) in gen_data
        if haskey(gen, "pg") && gen["pg"] > 0
            H_sys += gen["H"] * gen["pg"]
            total_Pg += gen["pg"]
        end
    end
    H_sys = total_Pg > 0 ? H_sys / total_Pg : 0.0

    
    # Add the inertia constraint to the model
    JuMP.@constraint(pm.model, H_sys >= H_min)
end
=#
