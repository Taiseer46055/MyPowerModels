# The constraints for inertia and reactive power are defined here.
function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id::Int, gen_tech::Int, delta_P::Float64, max_rocof::Float64)
    println("Add minimum system inertia constraint")
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)
    
    # Check if the bus number exists
    if !haskey(bus_data, bus_id)
        error("Bus number $bus_id does not exist in the network")
    end
    
    # Check if the required keys exist in gen_data
    required_keys = [:bus, :GenTech, :Pg, :H]
    for key in required_keys
        if !haskey(gen_data, key)
            error("Required key $key not found in Generator Data")
        end
    end
  
    # Find the specified generator at the specified bus
    gen_at_bus = nothing
    for (id, gen) in gen_data
        if gen[:bus] == bus_id && gen[:GenTech] == gen_tech
            gen_at_bus = gen
            break
        end
    end

    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end
    
    # Calculate H_min based on delta_P and max_rocof
    H_min = delta_P / max_rocof

    # Calculate H_sys
    H_sys = 0
    total_Pg = 0
    for (_, gen) in gen_data
        δ = gen[:Pg] > 0 ? 1 : 0
        H_sys += gen[:H] * gen[:Pg] * δ
        total_Pg += gen[:Pg] * δ
    end
    H_sys /= total_Pg

    # Add the constraint to JuMP
    JuMP.@constraint(pm.model, H_sys >= H_min)
end

