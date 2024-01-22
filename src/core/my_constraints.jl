function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id::Int, gen_tech::Int, delta_P::Float64, max_rocof::Float64)
    println("Add minimum system inertia constraint")
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)

    # Check if the bus number exists
    if !haskey(bus_data, bus_id)
        error("Bus number $bus_id does not exist in the network")
    end

    # Check if required keys exist in the dictionary
    for gen in values(gen_data)
        println("1")
        for key in required_keys
            println("2")
            if !haskey(gen, key)
                println("3")
                error("Required key $key not found in Generator with ID $(gen[:index])")
            end
        end
    end


    # Find the specified generator at the given bus
    gen_at_bus = findfirst(gen -> gen[:gen_bus] == bus_id && gen[:GenTech] == gen_tech, values(gen_data))

    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end

    # Calculate H_min based on delta_P and max_rocof
    H_min = delta_P / max_rocof

    # Calculate H_sys
    H_sys = 0
    total_Pg = 0
    for gen in values(gen_data)
        δ = gen[:Pg] > 0 ? 1 : 0
        H_sys += gen[:H] * gen[:Pg] * δ
        total_Pg += gen[:Pg] * δ
    end
    H_sys /= total_Pg

    # Add the constraint
    JuMP.@constraint(pm.model, H_sys >= H_min)
end
