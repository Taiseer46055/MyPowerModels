# The constraints for inertia and reactive power are defined here.

using DataFrames

function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id, gen_tech, delta_P, max_rocof)
  
    # Check if the required columns exist in the DataFrame
    required_columns = [:bus, :GenTech, :Pg, :H]
    for col in required_columns
        if !(col in names(pm.gen))
            error("Required column $col not found in pm.gen DataFrame")
        end
    end
  
    # Find the specified generator at the specified bus
    gen_at_bus = findfirst(row -> row[:bus] == bus_id && row[:GenTech] == gen_tech, eachrow(pm.gen))

    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end

    # Calculate H_min based on delta_P and max_rocof
    H_min = delta_P / max_rocof

    # Calculate H_sys
    H_sys = 0
    total_Pg = 0
    for gen in eachrow(pm.gen)
        δ = gen[:Pg] > 0 ? 1 : 0
        H_sys += gen[:H] * gen[:Pg] * δ
        total_Pg += gen[:Pg] * δ
    end
    H_sys /= total_Pg

    # Add the constraint
    JuMP.@constraint(pm.model, H_sys >= H_min)
end
