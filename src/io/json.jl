# Parse MyPowerModels data from JSON exports of MyPowerModels data structures.

function _jsonver2juliaver!(pm_data)
    if haskey(pm_data, "source_version") && isa(pm_data["source_version"], Dict)
        pm_data["source_version"] = "$(pm_data["source_version"]["major"]).$(pm_data["source_version"]["minor"]).$(pm_data["source_version"]["patch"])"
    end
end

"Parses json from iostream or string"
function parse_json(io::Union{IO,String}; kwargs...)::Dict{String,Any}
    pm_data = JSON.parse(io)

    _jsonver2juliaver!(pm_data)

    if haskey(pm_data, "conductors")
        Memento.warn(_LOGGER, "The JSON data contains the conductor parameter, but only single conductors are supported.  Consider using MyPowerModelsDistribution.")
    end

    if get(kwargs, :validate, true)
        MyPowerModels.correct_network_data!(pm_data)
    end

    return pm_data
end
