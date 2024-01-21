# The constraints for inertia and reactive power are defined here.

function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id::Int, gen_tech::Int, delta_P::Float64, max_rocof::Float64)
    println("Add minimum system inertia constraint")
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)
    
    # Chek if the busnumber exist
    if !haskey(bus_data, bus_id)
        error("Bus number $bus_id does not exist in the network")
    end
    
    # Check if the required columns exist in the DataFrame
    required_columns = [:bus, :GenTech, :Pg, :H]
    for col in required_columns
        if !(col in names(gen_data))
            error("Required column $col not found in Generator DataFrame")
        end
    end
  
    # Find the specified generator at the specified bus
    gen_at_bus = findfirst(row -> row[:bus] == bus_id && row[:GenTech] == gen_tech, eachrow(gen_data))

    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end
    
    # Calculate H_min based on delta_P and max_rocof
    H_min = delta_P / max_rocof

    # Calculate H_sys
    H_sys = 0
    total_Pg = 0
    for gen in eachrow(gen_data)
        δ = gen[:Pg] > 0 ? 1 : 0
        H_sys += gen[:H] * gen[:Pg] * δ
        total_Pg += gen[:Pg] * δ
    end
    H_sys /= total_Pg

    # Add the constraint
    JuMP.@constraint(pm.model, H_sys >= H_min)
end
