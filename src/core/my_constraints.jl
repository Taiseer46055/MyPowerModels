
function constraint_min_system_inertia(pm::AbstractACPModel, gen_id::Int, delta_P::Float64, max_rocof::Float64)
    # Aufruf der Funktion constraint_min_system_inertia aus acp.jl
    constraint_min_system_inertia(pm, gen_id, delta_P, max_rocof)
end


#=

function constraint_min_system_inertia(pm::AbstractPowerModel, gen_id::Int, delta_P::Float64, max_rocof::Float64)
    println("Adding minimum system inertia constraint")
    
    H_sys = variable_system_inertia(pm)
    
    # Retrieve generator and bus data
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)
    load_data = ref(pm, :load)

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
    for (_, load) in load_data
        if haskey(load, "pd")
            P_LOAD += load["pd"]
        end
    end
    
    if P_LOAD <= 0
        error("Total load P_LOAD is non-positive, which is invalid")
    end
    
    # Calculate the minimum system inertia H_min
    H_min = (delta_P * f0) / (P_LOAD * 2 * max_rocof)
    println(H_min)
    println(H_sys)
    
    # Add the inertia constraint to the model
    JuMP.@constraint(pm.model, H_sys >= H_min)
end


=#



#=    
    # Initialize a dictionary to store the inertia at each bus
    H_bus = Dict{Int, Float64}()
    
    # Iterate over all buses
    for (bus_id, bus) in bus_data
        # Initialize the inertia at this bus to 0
        H_bus[parse(Int, bus_id)] = 0.0
    end
    
    # Iterate over all generators
    for (_, gen) in gen_data
        # If the generator is at this bus, add its inertia to the bus inertia
        if haskey(gen, "pg") && gen["pg"] > 0
            bus_id = gen["gen_bus"]
            H_bus[bus_id] += gen["H"] * gen["pg"]
        end
    end


    H_sys = 0.0
    total_Pg = 0.0
    for (_, gen) in gen_data
        if  gen["pg"] > 0 #haskey(gen, "pg") &&
            H_sys += gen["H"] * gen["pmax"]
            total_Pg += gen["pmax"]
        end
    end
    
    
    #H_sys = total_Pg > 0 ? H_sys / total_Pg : 0.0
    H_sys = H_sys / total_Pg
=#





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
