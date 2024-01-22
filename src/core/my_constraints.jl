function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id::Int, gen_tech::Int, delta_P::Float64, max_rocof::Float64)
    println("Adding minimum system inertia constraint")
    
    # Retrieve generator and bus data
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)

    # Debugging: Output the content of gen_data and bus_data
    println("Generator data: ", gen_data)
    println("Bus data: ", bus_data)

    # Check if the bus number exists
    if !haskey(bus_data, bus_id)
        error("Bus number $bus_id does not exist in the network")
    end

    # Search for the specific generator
    gen_at_bus = nothing
    for (gen_id, gen) in gen_data
        println("Überprüfe Generator $gen_id am Bus $(gen["gen_bus"]) mit GenTech $(gen["GenTech"])")
        if string(gen["gen_bus"]) == string(bus_id) && string(gen["GenTech"]) == string(gen_tech)
            if haskey(gen, "Pg")
                gen_at_bus = gen
                gen_at_bus["Pg"] -= delta_P
                break
            else
                error("Key 'Pg' not found in generator with ID $gen_id")
            end
        end
    end

    if gen_at_bus === nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end

    # Calculate the minimum system inertia (H_min) and system inertia (H_sys)
    H_min = delta_P / max_rocof
    H_sys = 0.0
    total_Pg = 0.0
    for (_, gen) in gen_data
        if haskey(gen, "Pg") && gen["Pg"] > 0
            H_sys += gen["H"] * gen["Pg"]
            total_Pg += gen["Pg"]
        end
    end
    H_sys = total_Pg > 0 ? H_sys / total_Pg : 0.0

    # Add the inertia constraint to the model
    JuMP.@constraint(pm.model, H_sys >= H_min)
end

