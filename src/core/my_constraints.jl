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

    # Find the specified generator at the given bus
    gen_at_bus = findfirst(gen -> string(gen["gen_bus"]) == string(bus_id) && string(gen["GenTech"]) == string(gen_tech), gen_data)
    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end
    
    # Find the specified generator at the given bus
    gen_at_bus = findfirst(gen -> gen[:gen_bus] == bus_id && gen[:GenTech] == gen_tech, gen_data)
    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end

    # Calculate the minimum system inertia (H_min) based on delta_P and max_rocof
    H_min = delta_P / max_rocof

    # Calculate the system inertia (H_sys)
    H_sys = 0.0
    total_Pg = 0.0
    for gen in values(gen_data)
        # Only consider active generators
        δ = gen[:Pg] > 0 ? 1 : 0
        H_sys += gen[:H] * gen[:Pg] * δ
        total_Pg += gen[:Pg] * δ
    end
    H_sys /= total_Pg

    # Add the inertia constraint to the model
    JuMP.@constraint(pm.model, H_sys >= H_min)
end
