function constraint_min_system_inertia(pm::AbstractPowerModel, bus_id::Int, gen_tech::Int, delta_P::Float64, max_rocof::Float64)
    println("Adding minimum system inertia constraint")
    gen_data = ref(pm, :gen)
    bus_data = ref(pm, :bus)

    # Prüfen, ob die Busnummer existiert
    if !haskey(bus_data, bus_id)
        error("Bus number $bus_id does not exist in the network")
    end

    # Überprüfen, ob die erforderlichen Spalten in der Datenstruktur existieren
    required_columns = [:gen_bus, :GenTech, :Pg, :H]
    for col in required_columns
        if !(col in names(gen_data))
            error("Required column $col not found in Generator DataFrame")
        end
    end

    # Finden des spezifizierten Generators am angegebenen Bus
    gen_at_bus = findfirst(row -> row[:gen_bus] == bus_id && row[:GenTech] == gen_tech, eachrow(gen_data))

    if gen_at_bus == nothing
        error("No generator with GenTech $gen_tech found at bus $bus_id")
    end

    # Berechnen von H_min basierend auf delta_P und max_rocof
    H_min = delta_P / max_rocof

    # Berechnen von H_sys
    H_sys = 0
    total_Pg = 0
    for gen in eachrow(gen_data)
        δ = gen[:Pg] > 0 ? 1 : 0
        H_sys += gen[:H] * gen[:Pg] * δ
        total_Pg += gen[:Pg] * δ
    end
    H_sys /= total_Pg

    # Hinzufügen der Bedingung
    JuMP.@constraint(pm.model, H_sys >= H_min)
end
