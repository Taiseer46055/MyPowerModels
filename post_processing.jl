# Script to post-process the results of the optimization problem

using MyPowerModels
using JuMP
using Gurobi
using Plots
using Ipopt
using Juniper
using JLD2
using CSV
using DataFrames
using PowerPlots
using Graphs
using GraphPlot
using LightGraphs
using Glob


bus_system = "10"
relax_integrality = false

function load_all_results(bus_system, relax_integrality)

    path = "results\\results_bus_$(bus_system)\\relax_$(relax_integrality)"
    files = filter(f -> endswith(f, ".jld2"), readdir(path))
    println("Files: ", files)
    all_data = Dict()

    for file in files
        basename = split(file, '\\')[end]
        case_label = join(split(basename, '_')[1:2], '_')
        re_inj_int = split(split(basename, '.')[1], '_')[end]
        data = JLD2.load(joinpath(path, file))
    
        if haskey(data, "results")
            case_data = data["results"]
        else
            println("Key 'results' not found in file $file")
            continue
        end
        
        if haskey(all_data, case_label)
            all_data[case_label][re_inj_int] = case_data
        else
            all_data[case_label] = Dict(re_inj_int => case_data)
        end
    end

    mn_data_main = JLD2.load(joinpath(path, "mn_data_main.jld2"))["mn_data_main"]
    
    return all_data, mn_data_main
end

all_data, mn_data_main = load_all_results(bus_system, relax_integrality);

function add_missing_values_to_sol_data(all_data::Dict{Any, Any}, mn_data_main)

    for (case_label, results_v) in all_data
        for (re_inj_int, result) in results_v
            println("Adding missing values for case_label: $case_label, re_inj: $re_inj_int")
            if  result["results"]["termination_status"] == JuMP.OPTIMAL
                for (nw, data) in result["results"]["solution"]["nw"]
                    # Get gen_data and bus_data from mn_data_main
                    gen_data_mn = mn_data_main["nw"][nw]["gen"]
                    bus_data_mn = mn_data_main["nw"][nw]["bus"]
                    branch_data_mn = mn_data_main["nw"][nw]["branch"]

                    # Get gen_data and bus_data from solution
                    gen_data_solution = data["gen"]
                    bus_data_solution = data["bus"]
                    branch_data_solution = data["branch"]

                    # Merge gen_data
                    for (gen_id, gen) in gen_data_solution
                        if haskey(gen_data_mn, gen_id)
                            gen_data_solution[gen_id] = merge((x, y) -> ismissing(x) ? y : x, gen, gen_data_mn[gen_id])
                        else
                            println("No matching gen_id found in mn_data_main for gen_id: $gen_id")
                        end
                    end

                    # Merge bus_data
                    for (bus_id, bus) in bus_data_solution
                        if haskey(bus_data_mn, bus_id)
                            bus_data_solution[bus_id] = merge((x, y) -> ismissing(x) ? y : x, bus, bus_data_mn[bus_id])
                        else
                            println("No matching bus_id found in mn_data_main for bus_id: $bus_id")
                        end
                    end

                    # Merge branch_data
                    for (branch_id, branch) in branch_data_solution
                        if haskey(branch_data_mn, branch_id)
                            branch_data_solution[branch_id] = merge((x, y) -> ismissing(x) ? y : x, branch, branch_data_mn[branch_id])
                        else
                            println("No matching branch_id found in mn_data_main for branch_id: $branch_id")
                        end
                    end

                    # Update solution data with the merged data
                    data["gen"] = gen_data_solution
                    data["bus"] = bus_data_solution
                    data["branch"] = branch_data_solution

                end
            else
                println("No optimal solution found for case_label: $case_label, re_inj: $re_inj_int")
            end
        end
    end
    return all_data, mn_data_main
end

all_sol_data = add_missing_values_to_sol_data(all_data, mn_data_main)[1];
mn_data_main = add_missing_values_to_sol_data(all_data, mn_data_main)[2];


function create_dataframes(all_sol_data::Dict{Any, Any})
    dfs = Dict()

    zeitabhaengige_variablen = Set(["pg", "pf", "pt", "gen_status", "gen_startup", "gen_shutdown"])

    for (case_label, results_v) in all_sol_data
        for (re_inj_int, result) in results_v
            case_re_label = "$(case_label)_$(re_inj_int)"
            dfs[case_re_label] = Dict()
            if result["results"]["termination_status"] == JuMP.OPTIMAL
                for comp in ["gen", "bus", "branch"]
                    alle_vars = Set()
                    for (_, comp_data) in result["results"]["solution"]["nw"]
                        if haskey(comp_data, comp)
                            for (_, details) in comp_data[comp]
                                union!(alle_vars, keys(details))
                            end
                        end
                    end
                    nicht_zeitabhaengige_vars = setdiff(alle_vars, zeitabhaengige_variablen)

                    for var in zeitabhaengige_variablen
                        data = Dict()
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    if !haskey(data, nw)
                                        data[nw] = Dict()
                                    end
                                    data[nw][id] = get(details, var, NaN)
                                end
                            end
                        end
                        if isempty(data)
                            continue
                        end
                        nw_keys = collect(keys(data))
                        id_keys = collect(keys(data[first(nw_keys)]))
                        matrix = [get(get(data, nw, Dict()), id, NaN) for nw in nw_keys, id in id_keys]
                        df = DataFrame(matrix, Symbol.(id_keys))
                        insertcols!(df, 1, :nw => nw_keys)
                        dfs[case_re_label]["$(comp)_$(var)"] = df
                    end

                    for var in nicht_zeitabhaengige_vars
                        simple_data = []
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    push!(simple_data, (nw=nw, id=id, value=get(details, var, NaN)))
                                end
                            end
                        end
                        if isempty(simple_data)
                            continue
                        end
                        df = DataFrame(simple_data, [:nw, :id, Symbol(var)])
                        dfs[case_re_label]["$(comp)_$(var)"] = df
                    end
                end
            end
        end
    end

    return dfs
end

dfs = create_dataframes(all_sol_data)














#=
function create_dataframes(all_sol_data::Dict{Any, Any})
    # Dictionary, das alle DataFrames enthält
    dfs = Dict()

    # Definieren der Menge von zeitabhängigen Variablen
    zeitabhaengige_variablen = Set(["pg", "pf", "pt", "gen_status", "gen_startup", "gen_shutdown"])

    for (case_label, results_v) in all_sol_data
        case_dict = Dict()
        
        for (re_inj_int, result) in results_v
            re_inj_dict = Dict()

            if result["results"]["termination_status"] == JuMP.OPTIMAL
                for comp in ["gen", "bus", "branch"]
                    # Alle Variablen für die Komponente sammeln
                    alle_vars = Set()
                    for (_, comp_data) in result["results"]["solution"]["nw"]
                        if haskey(comp_data, comp)
                            for (_, details) in comp_data[comp]
                                union!(alle_vars, keys(details))
                            end
                        end
                    end
                    nicht_zeitabhaengige_vars = setdiff(alle_vars, zeitabhaengige_variablen)

                    # Verarbeitung der zeitabhängigen Variablen
                    for var in zeitabhaengige_variablen
                        data = Dict()
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    if !haskey(data, nw)
                                        data[nw] = Dict()
                                    end
                                    data[nw][id] = get(details, var, NaN)
                                end
                            end
                        end
                        if isempty(data)
                            continue
                        end
                        nw_keys = collect(keys(data))
                        id_keys = collect(keys(data[first(nw_keys)]))
                        matrix = [get(get(data, nw, Dict()), id, NaN) for nw in nw_keys, id in id_keys]
                        df = DataFrame(matrix, Symbol.(id_keys))
                        insertcols!(df, 1, :nw => nw_keys)
                        re_inj_dict["$(comp)_$(var)"] = df
                    end

                    # Verarbeitung der nicht-zeitabhängigen Variablen
                    for var in nicht_zeitabhaengige_vars
                        simple_data = []
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    push!(simple_data, (nw=nw, id=id, value=get(details, var, NaN)))
                                end
                            end
                        end
                        if isempty(simple_data)
                            continue
                        end
                        df = DataFrame(simple_data, [:nw, :id, Symbol(var)])
                        re_inj_dict["$(comp)_$(var)"] = df
                    end
                end
            end
            case_dict["$re_inj_int"] = re_inj_dict
        end
        dfs["$case_label"] = case_dict
    end

    return dfs
end

dfs_t = create_dataframes(all_sol_data)

=#


#=
function create_dataframes(all_sol_data::Dict{Any, Any})
    df_dict = Dict()
    zeitabhaengige_variablen = Set(["pg", "pf", "pt", "gen_status", "gen_startup", "gen_shutdown"])

    for (case_label, results_v) in all_sol_data
        for (re_inj_int, result) in results_v
            if result["results"]["termination_status"] == JuMP.OPTIMAL
                for comp in ["gen", "bus", "branch"]
                    # Alle Variablen sammeln
                    alle_vars = Set()
                    for (_, comp_data) in result["results"]["solution"]["nw"]
                        if haskey(comp_data, comp)
                            for (_, details) in comp_data[comp]
                                union!(alle_vars, keys(details))
                            end
                        end
                    end
                    nicht_zeitabhaengige_vars = setdiff(alle_vars, zeitabhaengige_variablen)

                    # Verarbeitung der zeitabhängigen Variablen
                    for var in zeitabhaengige_variablen
                        data = Dict()
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    if !haskey(data, nw)
                                        data[nw] = Dict()
                                    end
                                    data[nw][id] = get(details, var, NaN)
                                end
                            end
                        end
                        if isempty(data)
                            continue
                        end
                        nw_keys = collect(keys(data))
                        id_keys = collect(keys(data[first(nw_keys)]))
                        matrix = [get(get(data, nw, Dict()), id, NaN) for nw in nw_keys, id in id_keys]
                        df = DataFrame(matrix, Symbol.(id_keys))
                        insertcols!(df, 1, :nw => nw_keys)
                        df_dict["$(case_label)_$(re_inj_int)_$(comp)_$(var)"] = df
                    end

                    # Verarbeitung der nicht-zeitabhängigen Variablen
                    for var in nicht_zeitabhaengige_vars
                        simple_data = []
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    push!(simple_data, (nw=nw, id=id, value=get(details, var, NaN)))
                                end
                            end
                        end
                        if isempty(simple_data)
                            continue
                        end
                        df = DataFrame(simple_data, [:nw, :id, Symbol(var)])
                        df_dict["$(case_label)_$(re_inj_int)_$(comp)_$(var)_simple"] = df
                    end
                end
            end
        end
    end

    return df_dict
end

dfs_t = create_dataframes(all_sol_data)

=#






#=
function load_all_results(bus_system, relax_integrality)

    path = "results\\results_bus_$(bus_system)\\relax_$(relax_integrality)"
    files = glob("*.jld2", path)
    
    all_cases = Dict()

    for file in files
        case_label = split(split(file, '\\')[end], '.')[1]
        case_data = JLD2.load(file)["results"]
        all_cases[case_label] = case_data
    end
    
    return all_cases
end

all_results = load_all_results(bus_system, relax_integrality)



function add_missing_values_to_mn_data(all_results::Dict{Any, Any})
    for (case_label, results_v) in all_results
        for (i_string, result) in results_v
            println("Adding missing values for case_label: $case_label, re_inj: $i_string")
            println("Results: ", keys(result["results"]))
            if result["results"]["termination_status"] == JuMP.OPTIMAL
                for (nw, data) in result["results"]["solution"]["nw"]
                    # Get gen_data and bus_data from mn_data
                    gen_data_mn = result["mn_data"]["nw"][nw]["gen"]
                    bus_data_mn = result["mn_data"]["nw"][nw]["bus"]

                    # Get gen_data and bus_data from solution
                    gen_data_solution = data["gen"]
                    bus_data_solution = data["bus"]

                    # Merge gen_data
                    for (gen_id, gen) in gen_data_mn
                        if haskey(gen_data_solution, gen_id)
                            gen_data_mn[gen_id] = merge(gen, gen_data_solution[gen_id])
                        else
                            println("No matching gen_id found in sol_gen_data for gen_id: $gen_id")
                        end
                    end

                    # Merge bus_data
                    for (bus_id, bus) in bus_data_mn
                        if haskey(bus_data_solution, bus_id)
                            bus_data_mn[bus_id] = merge(bus, bus_data_solution[bus_id])
                        else
                            println("No matching bus_id found in sol_bus_data for bus_id: $bus_id")
                        end
                    end

                    # Update mn_data with the merged data
                    result["mn_data"]["nw"][nw]["gen"] = gen_data_mn
                    result["mn_data"]["nw"][nw]["bus"] = bus_data_mn
                end
            else
                println("No optimal solution found for case_label: $case_label, re_inj: $i_string")
            end
        end
    end
    return all_results
end

all_data = add_missing_values_to_mn_data(all_results)



function convert_all_data_to_powermodels_dataframe(all_sol_data, re_inj_int)
    dfs_gen = Dict{String, DataFrame}()
    dfs_bus = Dict{String, DataFrame}()
    dfs_branch = Dict{String, DataFrame}()
    
    for (case_label, results_v) in all_sol_data
        for (re_inj, result) in results_v
            df_gen = DataFrame()
            df_bus = DataFrame()
            df_branch = DataFrame()

            for (nw, data) in result["results"]["solution"]["nw"]
                gen_data = data["gen"]
                bus_data = data["bus"]
                branch_data = data["branch"]

                for (gen_key, gen_values) in gen_data
                    df = DataFrame(value = gen_values)
                    df[!, :nw] = repeat([nw], nrow(df))
                    df[!, :gen_key] = repeat([gen_key], nrow(df))
                    df_gen = vcat(df_gen, df)
                end

                for (bus_key, bus_values) in bus_data
                    df = DataFrame(value = bus_values)
                    df[!, :nw] = repeat([nw], nrow(df))
                    df[!, :bus_key] = repeat([bus_key], nrow(df))
                    df_bus = vcat(df_bus, df)
                end

                for (branch_key, branch_values) in branch_data
                    df = DataFrame(value = branch_values)
                    df[!, :nw] = repeat([nw], nrow(df))
                    df[!, :branch_key] = repeat([branch_key], nrow(df))
                    df_branch = vcat(df_branch, df)
                end
            end

            dfs_gen["$(case_label)_$(re_inj)"] = df_gen
            dfs_bus["$(case_label)_$(re_inj)"] = df_bus
            dfs_branch["$(case_label)_$(re_inj)"] = df_branch
        end
    end
    return dfs_gen, dfs_bus, dfs_branch
end

dfs_gen, dfs_bus, dfs_branch = convert_all_data_to_powermodels_dataframe(all_sol_data, "re_inj_int")

=#







#=
function add_calculated_values_to_mn_data(all_results)
 

    if case_results["results"]["termination_status"] == JuMP.OPTIMAL

        for nw in sort(collect(keys(all_data["data"]["mn_data"]["nw"])))

            results[test_case] = case_results["results"]
            options[test_case] =  case_results["options"]
            gen_data = mn_data["nw"][nw]["gen"]
            bus_data = mn_data["nw"][nw]["bus"]
            load_data = mn_data["nw"][nw]["load"]
            f0 = 50
            rocof = options[test_case]["f"]["rocof"]

            sol_gen_data = case_results["results"]["solution"]["nw"][nw]["gen"]
            sol_bus_data = case_results["results"]["solution"]["nw"][nw]["bus"]


            for (gen_id, gen) in gen_data
                if haskey(sol_gen_data, gen_id)
                    gen_data[gen_id] = merge(gen_data[gen_id], sol_gen_data[gen_id])
                else
                    println("No matching gen_id found in sol_gen_data for gen_id: $gen_id")
                end
            end

            for (bus_id, bus) in bus_data
                if haskey(sol_bus_data, bus_id)
                    bus_data[bus_id] = merge(bus_data[bus_id], sol_bus_data[bus_id])
                else
                    println("No matching bus_id found in sol_bus_data for bus_id: $bus_id")
                end
            end
            
            # Update mn_data with the updated gen_data and bus_data
            mn_data["nw"][nw]["gen"] = gen_data
            mn_data["nw"][nw]["bus"] = bus_data
    
            for i in keys(gen_data)
                if gen_data[i]["gen_status"] < 1e-2
                    gen_data[i]["gen_status"] = 0
                end
                if gen_data[i]["pg"] < 1e-2
                    gen_data[i]["gen_status"] = 0
                end
            end

        end
    end
    return all_data
end

add_calculated_values_to_mn_data(results_m)

for test_case in test_cases

    if case_results["Results"]["termination_status"] == JuMP.OPTIMAL

        total_pg = 0.0
        renewable_injection = 0.0
        renewable_carriers = [3, 4, 5]
        for nw in sort(collect(keys(case_results["Results"]["solution"]["nw"])))
            gen_data = mn_data["nw"][nw]["gen"]
            total_pg += sum(gen["pg"] for (id, gen) in gen_data)
            renewable_injection += sum(gen["pg"] for (id, gen) in gen_data if gen["carrier"] in renewable_carriers)
        end
        print("Total_PG_$test_case: $total_pg ")
        if total_pg > 0
            re_inj = renewable_injection / total_pg
            case_results["Results"]["solution"]["re_inj"] = re_inj
        else
            println("Total generation (total_pg) for case $(test_case) is 0, therefore no calculation of the renewable share is possible.")
        end
    end
end
# create a PowerModelsDataFrame from the results for all test cases

function create_power_models_data_frames(all_data, test_cases)
    pmdfs = Dict()

    for test_case in test_cases
        for nw in sort(collect(keys(all_data["data"]["mn_data"]["nw"])))
            # Create a PowerModelsDataFrame for each network and case
            pm_data_frame = PowerModelsDataFrame(all_data["data"]["mn_data"]["nw"][nw])
            pmdfs["$test_case-nw_$nw"] = pm_data_frame
        end
    end

    return pmdfs
end


pmdfs = create_power_models_data_frames(all_data, test_cases)



=#

























#     pmdfs = Dict{String, Dict{String, PowerModelsDataFrame}}()
#     for (test_case, case_data) in test_cases
#         pmdfs[test_case] = Dict{String, PowerModelsDataFrame}()

#         for (nw_id, nw_data) in all_data["data"]["mn_data"]["nw"]

#             data_for_pmdf = Dict{String, Any}()
#             data_for_pmdf["gen"] =  all_data["cases"][test_cases]["Results"]["solution"]["nw"][nw_id]["gen"]
#             data_for_pmdf["bus"] = all_data["cases"][test_cases]["Results"]["solution"]["nw"][nw_id]["bus"]
#             data_for_pmdf["branch"] = all_data["cases"][test_cases]["Results"]["solution"]["nw"][nw_id]["branch"]

#             pmdfs[test_case][nw_id] = PowerModelsDataFrame(data_for_pmdf)
#         end
#     end
#     return pmdfs
# end

#=
function create_pmdf_for_solution(test_cases, all_data)
    pmdf = PowerModelsDataFrame(Dict{String, Any}())
    for test_case in test_cases
        println("Creating PowerModelsDataFrame for test case: $test_case")
        solution_data = case_results["Results"]["solution"]
        solution_data_string = Dict{String, Any}(solution_data)
        pmdf = PowerModelsDataFrame(solution_data_string)
        case_results["Results"]["solution"]["PowerModelsDataFrame"] = pmdf
    end
    return pmdf
end
pmdf = create_pmdf_for_solution(test_cases, all_data)












function load_results(case_name, is_multi_network)
    if is_multi_network
        results_filename = "results\\multi_network_results\\results_$case_name.jld2"
        data_filename = "results\\multi_network_results\\data_$case_name.jld2"
        options_filename = "results\\multi_network_results\\options_$case_name.jld2"
        mn_data = JLD2.load(data_filename)["data"]
        mn_options = JLD2.load(options_filename)["options"]
        mn_results = JLD2.load(results_filename)["results"]
        return mn_data, mn_options, mn_results
    else
        results_filename = "results\\single_network_results\\results_$case_name.jld2"
        data_filename = "results\\single_network_results\\data_$case_name.jld2"
        options_filename = "results\\single_network_results\\options_$case_name.jld2"
        sn_data = JLD2.load(data_filename)["data"]
        sn_options = JLD2.load(options_filename)["options"]
        sn_results = JLD2.load(results_filename)["results"]
        return sn_data, sn_options, sn_results
    end
end


if is_multi_network
    mn_data, mn_options, mn_results = load_results(case_name)
    sorted_keys = sort([parse(Int, n) for n in keys(mn_results)])
    num_networks = length(sorted_keys)
    
    for n in sorted_keys
        nw_results = mn_results[string(n)]
        mn_gen_data = mn_data["nw"][string(n)]["gen"]
        mn_bus_data = mn_data["nw"][string(n)]["bus"]
        mn_branch_data = mn_data["nw"][string(n)]["branch"]
        mn_load_data = mn_data["nw"][string(n)]["load"]
        mn_baseMVA = mn_data["nw"][string(n)]["baseMVA"]
        mn_sol_gen_data = nw_results["gen"]
        mn_sol_bus_data = nw_results["bus"]
        local options = mn_options["f"]
        local alpha = options["alpha_factor"]
        local calc_delta_P = options["calc_delta_P"]
        local rocof = options["rocof"]
        local bus = options["bus"]
        local area = options["area"]
        local weighted_area = options["weighted_area"]

        
        for (gen_id, gen) in mn_gen_data
            H = mn_gen_data[gen_id]["H"]
            pmin = mn_gen_data[gen_id]["pmin"]
            pmax = mn_gen_data[gen_id]["pmax"]
            pmin = abs(pmin)
            pmax = abs(pmax)
            pv = [pmin, pmax]
            pmin = minimum(pv)
            pmax = maximum(pv)
            gen_bus = mn_gen_data[gen_id]["gen_bus"]

            if !haskey(mn_sol_gen_data, gen_id)
                mn_sol_gen_data[gen_id] = Dict()
            end

            mn_sol_gen_data[gen_id]["H"] = H
            mn_sol_gen_data[gen_id]["pmin"] = pmin
            mn_sol_gen_data[gen_id]["pmax"] = pmax
            mn_sol_gen_data[gen_id]["gen_bus"] = gen_bus
        end

        for (bus_id, bus) in mn_bus_data
            local area = mn_bus_data[bus_id]["area"]
            mn_sol_bus_data[bus_id]["area"] = area
        end

        for i in keys(mn_sol_gen_data)
            if mn_sol_gen_data[i]["gen_status"] < 1e-2
                mn_sol_gen_data[i]["gen_status"] = 0
            end
            if mn_sol_gen_data[i]["pg"] < 1e-2
                mn_sol_gen_data[i]["gen_status"] = 0
            end
        end

        E_I_sys = sum(mn_sol_gen_data[i]["H"] * mn_sol_gen_data[i]["pmax"] * mn_sol_gen_data[i]["gen_status"] for i in keys(mn_sol_gen_data))
        local delta_P = mn_data["delta_P"]    
        E_I_min = mn_data["E_I_min"]

        E_I_weighted_area = Dict()

        if options["weighted_area"] == "load" && options["disturbance"] == "small"

            areas = unique([mn_sol_bus_data[j]["area"] for j in keys(mn_sol_bus_data)])
            local sum_E_I_area_weighted = 0
            for area in areas
                gens_in_area = [i for i in keys(mn_sol_gen_data) if mn_sol_bus_data[string(mn_sol_gen_data[i]["gen_bus"])]["area"] == area]
                E_I_weighted_area[area] = sum(mn_sol_gen_data[i]["H"] * mn_sol_gen_data[i]["pmax"] * mn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
            end
        end

        if options["weighted_area"] == "equal" && options["disturbance"] == "small"

            areas = unique([mn_sol_bus_data[j]["area"] for j in keys(mn_sol_bus_data)])
            for area in areas
                gens_in_area = [i for i in keys(mn_sol_gen_data) if mn_sol_bus_data[string(mn_sol_gen_data[i]["gen_bus"])]["area"] == area]
                E_I_weighted_area[area] = sum(mn_sol_gen_data[i]["H"] * mn_sol_gen_data[i]["pmax"] * mn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
            end
        end

        if options["weighted_area"] == "none" && options["disturbance"] == "small"

            areas = unique([mn_sol_bus_data[j]["area"] for j in keys(mn_sol_bus_data)])
            for area in areas
                E_I_weighted_area[area] = missing
            end
        end

        # Berechnung für jeden Bus
        local f0 = 50
        local H_min_bus = Dict()
        local E_I_bus = Dict()
        local P_gen_bus = Dict()
        local P_load_bus = Dict()
        local delta_p_bus = Dict()

        if options["bus"] == "true" && options["disturbance"] == "large"
            E_I_bus = Dict()
            E_I_min_bus = Dict()
            P_load_bus = Dict()
            delta_p_bus = Dict()
            P_gen_bus = Dict()

            for j in keys(mn_bus_data)
                P_gen_bus[j] = sum(mn_sol_gen_data[i]["pg"] for (i, gen) in ref(pm, n, :gen) if gen["gen_bus"] == j)
                P_load_bus[j] = sum(mn_load_data[i]["pd"] for i in eachindex(mn_load_data) if mn_load_data[i]["load_bus"] == j; init=0)
                delta_p_bus[j] = abs(P_gen_bus[j] - P_load_bus[j])
                E_I_min_bus[j] = (delta_p_bus[j] * f0 / rocof * 2)
                E_I_bus[j] = sum(mn_sol_gen_data[i]["H"] * mn_sol_gen_data[i]["pmax"] * mn_sol_gen_data[i]["gen_status"] for i in eachindex(mn_sol_gen_data) if mn_sol_gen_data[i]["gen_bus"] == j; init=0)

            end
        end

        # Berechnung für jede Area
        local E_I_min_area = Dict()
        local E_I_area = Dict()
        local P_gen_area = Dict()
        local P_load_area = Dict()
        local delta_p_area = Dict()

        global areas = unique([mn_bus_data[j]["area"] for j in keys(mn_bus_data)])

        if options["area"] == "true" && options["disturbance"] == "large"
            E_I_area = Dict()
            E_I_min_area = Dict()
            P_load_area = Dict()
            delta_p_area = Dict()
            P_gen_area = Dict()

            areas = unique([mn_bus_data[j]["area"] for j in keys(mn_bus_data)])
            for area in areas

                gens_in_area = [i for i in keys(mn_sol_gen_data) if mn_sol_bus_data[string(mn_sol_gen_data[i]["gen_bus"])]["area"] == area]
                buses_in_area = [i for i in keys(mn_sol_bus_data) if mn_sol_bus_data[i]["area"] == area]
                
                P_load_area[area] = sum(mn_load_data[i]["pd"] for i in eachindex(mn_load_data) if mn_load_data[i]["load_bus"] in buses_in_area; init=0)
                P_gen_area = sum(mn_sol_gen_data[i]["pg"] for i in gens_in_area)
                delta_p_area[area] = abs(P_gen_area - P_load_area[area])
                E_I_min_area[area] = (delta_p_area[area] * f0/ rocof * 2)
                E_I_area[area] = sum(mn_sol_gen_data[i]["H"] * mn_sol_gen_data[i]["pmax"] * mn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
            end
        end
        local df = DataFrame(
            Network = String[],  # Neue Spalte für das Netzwerk
            Options = String[],
            rocof = Float64[],
            H_Sys = Float64[],
            Delta_P = Float64[],
            H_Min = Float64[],
            Total_cost = Float64[]
        )

        function generate_options_name(options)
            option_parts = []

            if options["inertia_constraint"] == "true"
                
                if options["system"] == "true"
                    push!(option_parts, "system")
                end

                if options["disturbance"] == "small"
                    push!(option_parts, "small")
                    if options["weighted_area"] in ["load", "equal", "none"]
                        push!(option_parts, options["weighted_area"])
                    end
                elseif options["disturbance"] == "large"
                    push!(option_parts, "large")
                    if options["bus"] == "true"
                        push!(option_parts, "bus")
                    end
                    if options["area"] == "true"
                        push!(option_parts, "area")
                    end
                end

            else 
                push!(option_parts, "unit_commetment")
            end

            return join(option_parts, "_")

        end
        
        function save_results_as_csv(df, options, results, network, num_networks)
            
            options_name = generate_options_name(options["f"])
            row = Dict{Symbol, Any}(:Network => network, :Options => options_name)
            total_cost = 0
            for (gen_id, gen) in nw_results["gen"]
                pg = gen["pg"]
                pg_cost = gen["pg_cost"]
                su_cost = gen["startup_cost"]
                sd_cost = gen["shutdown_cost"]
                if network == "0"
                    investment_cost = gen["investment_cost"]
                else
                    investment_cost = 0
                end
                # investment_cost = gen["investment_cost"]
                row[Symbol("PG_$gen_id")] = pg
                row[Symbol("PG_Cost_$gen_id")] = pg_cost
                row[Symbol("SU_Cost_$gen_id")] = su_cost
                row[Symbol("SD_Cost_$gen_id")] = sd_cost
                row[Symbol("Investment_Cost_$gen_id")] = investment_cost

                total_cost += pg_cost + su_cost + sd_cost + investment_cost
            end
            row[Symbol("Total_cost")] = total_cost
            row[Symbol("E_I_sys")] = E_I_sys
            row[Symbol("Delta_P")] = delta_P
            row[Symbol("E_I_min")] = E_I_min
            row[Symbol("rocof")] = rocof

            if options["weighted_area"] == "load" && options["disturbance"] == "small"
                for area in areas
                    column_name = Symbol("H_Weighted_Area_$area")
                    row[column_name] = E_I_weighted_area[area]
                end
            end

            if options["weighted_area"] == "equal" && options["disturbance"] == "small"
                for area in areas
                    column_name = Symbol("H_Weighted_Area_$area")
                    row[column_name] = E_I_weighted_area[area]
                end
            end

            if options["weighted_area"] == "none" && options["disturbance"] == "small"
                for area in areas
                    column_name = Symbol("H_Weighted_Area_$area")
                    row[column_name] = E_I_weighted_area[area]
                end
            end

            if options["bus"] == "true" && options["disturbance"] == "large"
                for j in keys(mn_sol_bus_data)
                    row[Symbol("E_I_Bus_$j")] = E_I_bus[j]
                    row[Symbol("Delta_P_Bus_$j")] = delta_p_bus[j]
                    row[Symbol("E_I_min_bus$j")] = E_I_min_bus[j]
                end
            end

            if options["area"] == "true" && options["disturbance"] == "large"
                for r in areas
                    row[Symbol("E_I_Area_$r")] = E_I_area[r]
                    row[Symbol("Delta_P_Area_$r")] = delta_p_area[r]
                    row[Symbol("E_I_min_area$r")] = E_I_min_area[r]
                end
            end
            
            df_row = DataFrame(row)
            df_row = select(df_row, :Options, :Network, Not([:Options, :Network]))

            results_filename = "results\\multi_network_results\\results_$(case_name)_$(num_networks)_networks.csv"
            if isfile(results_filename)
                existing_df = CSV.read(results_filename, DataFrame)
                for col_name in names(df_row)
                    if !hasproperty(existing_df, col_name)
                        existing_df[!, col_name] = fill(missing, nrow(existing_df))
                    end
                end
                combined_df = vcat(existing_df, df_row, cols=:union)
                combined_df = select(combined_df, :Options, :Network, Not([:Options, :Network]))    
                CSV.write(results_filename, combined_df)
            else
                df_row = select(df_row, :Options, :Network, Not([:Options, :Network]))      
                CSV.write(results_filename, df_row)
            end
        end
        
        save_results_as_csv(df, mn_options, nw_results, string(n), num_networks)
        
    end

###########################################################################
###########################################################################

else

    sn_data, sn_options, sn_results = load_results(case_name, is_multi_network)
    sn_gen_data = sn_data["gen"]
    sn_bus_data = sn_data["bus"]
    sn_branch_data = sn_data["branch"]
    sn_load_data = sn_data["load"]
    sn_baseMVA = sn_data["baseMVA"]
    sn_sol_gen_data = sn_results["gen"]
    sn_sol_bus_data = sn_results["bus"]
    options = sn_options["f"]
    alpha = options["alpha_factor"]
    calc_delta_P = options["calc_delta_P"]
    rocof = options["rocof"]
    bus = options["bus"]
    area = options["area"]
    weighted_area = options["weighted_area"]


    for (gen_id, gen) in sn_gen_data
        H = sn_data["gen"][gen_id]["H"]

        pmin = sn_data["gen"][gen_id]["pmin"]
        pmax = sn_data["gen"][gen_id]["pmax"]
        gen_bus = sn_data["gen"][gen_id]["gen_bus"]

        if !haskey(sn_sol_gen_data, gen_id)
            sn_sol_gen_data[gen_id] = Dict()
        end

        sn_sol_gen_data[gen_id]["H"] = H
        sn_sol_gen_data[gen_id]["pmin"] = pmin
        sn_sol_gen_data[gen_id]["pmax"] = pmax
        sn_sol_gen_data[gen_id]["gen_bus"] = gen_bus
    end

    for (bus_id, bus) in sn_bus_data
    
        local area = sn_data["bus"][bus_id]["area"]
        sn_sol_bus_data[bus_id]["area"] = area

    end

    for i in keys(sn_sol_gen_data)
        if sn_sol_gen_data[i]["gen_status"] < 1e-2
            sn_sol_gen_data[i]["gen_status"] = 0
        end
        if sn_sol_gen_data[i]["pg"] < 1e-2
            sn_sol_gen_data[i]["gen_status"] = 0
        end
    end
    H_sys = sum(sn_sol_gen_data[i]["H"] * sn_sol_gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in eachindex(sn_sol_gen_data)) / sum(sn_sol_gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in eachindex(sn_sol_gen_data))
    println("H_sys = ", H_sys)

    delta_P = sn_data["delta_P"]
    println("delta_P = ", delta_P)

    H_min = sn_data["H_min"]
    println("H_min = ", H_min)

    H_min_weighted_area = Dict()
    P_gen_weighted_area = Dict()
    P_load_weighted_area = Dict()
    delta_p_weighted_area = Dict()
    H_weighted_area = Dict()

    if options["weighted_area"] == "load" && options["disturbance"] == "small"

        areas = unique([sn_bus_data[j]["area"] for j in keys(sn_bus_data)])
        local sum_E_I_area_weighted = 0
        for area in areas

            W_v = Dict()
            load_sum_area = Dict()
            gens_in_area = [i for i in keys(sn_sol_gen_data) if haskey(sn_sol_bus_data, string(sn_sol_gen_data[i]["gen_bus"])) && sn_sol_bus_data[string(sn_sol_gen_data[i]["gen_bus"])]["area"] == area]
            loads_in_area = [i for i in keys(sn_load_data) if sn_sol_bus_data[string(sn_load_data[i]["load_bus"])]["area"] == area]
            load_sum_area[area] = sum(sn_load_data[i]["pd"] for i in loads_in_area; init=0)
            E_I_weighted_area[area] = sum(sn_sol_gen_data[i]["H"] * sn_sol_gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
            P_load = sum(sn_load_data[i]["pd"] for i in keys(sn_load_data); init=0)
            W_v[area] = load_sum_area[area] / P_load
            sum_E_I_area_weighted += W_v[area] * H_weighted_area[area]
    
        end
    end

    if options["weighted_area"] == "equal" && options["disturbance"] == "small"

        areas = unique([sn_bus_data[j]["area"] for j in keys(sn_bus_data)])
        for area in areas
            gens_in_area = [i for i in keys(sn_sol_gen_data) if haskey(sn_sol_bus_data, string(sn_sol_gen_data[i]["gen_bus"])) && sn_sol_bus_data[string(sn_sol_gen_data[i]["gen_bus"])]["area"] == area]
            H_weighted_area[area] = sum(sn_sol_gen_data[i]["H"] * sn_sol_gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0) / sum(sn_sol_gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
        end
    end


    if options["weighted_area"] == "none" && options["disturbance"] == "small"

        for area in areas
            H_weighted_area[area] = missing
        end
    end

    # Berechnung für jeden Bus
    f0 = 50
    H_min_bus = Dict()
    E_I_bus = Dict()
    P_gen_bus = Dict()
    P_load_bus = Dict()
    delta_p_bus = Dict()

    if options["bus"] == "true" && options["disturbance"] == "large"
        for j in keys(sn_sol_bus_data)

            gens_at_bus = [i for i in keys(sn_sol_gen_data) if string(sn_sol_gen_data[i]["gen_bus"]) == j]
            P_gen_bus[j] = sum(sn_sol_gen_data[i]["pg"] for i in gens_at_bus; init=0)
            P_load_bus[j] = sum(sn_load_data[i]["pd"] for i in keys(sn_load_data) if string(sn_load_data[i]["load_bus"]) == j; init=0)
            delta_p_bus[j] = abs(P_gen_bus[j] - P_load_bus[j])
            H_min_bus[j] = (delta_p_bus[j] * f0) / (P_load_bus[j] * 2 * rocof)
            E_I_bus[j] = sum(sn_sol_gen_data[i]["H"] * sn_sol_gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in gens_at_bus; init=0)


        end
    end

    # Berechnung für jede Area
    E_I_min_area = Dict()
    E_I_area = Dict()
    P_gen_area = Dict()
    P_load_area = Dict()
    delta_p_area = Dict()

    global areas = unique([sn_bus_data[j]["area"] for j in keys(sn_bus_data)])

    if options["area"] == "true" && options["disturbance"] == "large"
        for area in areas

            gens_in_area = [i for i in keys(sn_sol_gen_data) if haskey(sn_sol_bus_data, string(sn_sol_gen_data[i]["gen_bus"])) && sn_sol_bus_data[string(sn_sol_gen_data[i]["gen_bus"])]["area"] == area]
            buses_in_area = [i for i in keys(sn_bus_data) if sn_bus_data[i]["area"] == area]

            P_load_area[area] = sum(sn_load_data[i]["pd"] for i in keys(sn_load_data) if string(sn_load_data[i]["load_bus"]) in buses_in_area; init=0)
            P_gen_area[area] = sum(sn_sol_gen_data[i]["pg"] for i in gens_in_area; init=0)
            delta_p_area[area] = abs(P_gen_area[area] - P_load_area[area])
            E_I_min_area[area] = (delta_p_area[area] * f0) / (P_load_area[area] * 2 * rocof)
            E_I_area[area] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * sn_sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0) 
            end
        end
        df = DataFrame(
        Options = String[],
        rocof = Float64[],
        H_Sys = Float64[],
        Delta_P = Float64[],
        H_Min = Float64[],
        Total_cost = Float64[]

        )

    function generate_options_name(options)
        option_parts = []

        if options["inertia_constraint"] == "true"
            
            if options["system"] == "true"
                push!(option_parts, "system")
            end

            if options["disturbance"] == "small"
                push!(option_parts, "small")
                if options["weighted_area"] in ["load", "equal", "none"]
                    push!(option_parts, options["weighted_area"])
                end
            elseif options["disturbance"] == "large"
                push!(option_parts, "large")
                if options["bus"] == "true"
                    push!(option_parts, "bus")
                end
                if options["area"] == "true"
                    push!(option_parts, "area")
                end
            end

        else 
            push!(option_parts, "unit_commetment")
        end

        return join(option_parts, "_")

    end

    function save_results_as_csv(df, options, results)
        options_name = generate_options_name(options["f"])
        row = Dict{Symbol, Any}(:Options => options_name)

        total_cost = 0
        for (gen_id, gen) in results["gen"]
            pg = gen["pg"]
            pg_cost = gen["pg_cost"]
            row[Symbol("PG_$gen_id")] = pg
            row[Symbol("PG_Cost_$gen_id")] = pg_cost
            total_cost += pg_cost
        end
        row[Symbol("Total_cost")] = total_cost
        row[Symbol("H_Sys")] = H_sys
        row[Symbol("Delta_P")] = delta_P
        row[Symbol("H_Min")] = H_min
        row[Symbol("rocof")] = rocof

        if options["weighted_area"] == "load" && options["disturbance"] == "small"
            for area in areas
                column_name = Symbol("H_Weighted_Area_$area")
                row[column_name] = H_weighted_area[area]
            end
        end

        if options["weighted_area"] == "equal" && options["disturbance"] == "small"
            for area in areas
                column_name = Symbol("H_Weighted_Area_$area")
                row[column_name] = H_weighted_area[area]
            end
        end

        if options["weighted_area"] == "none" && options["disturbance"] == "small"
            for area in areas
                column_name = Symbol("H_Weighted_Area_$area")
                row[column_name] = H_weighted_area[area]
            end
        end

        if options["bus"] == "true" && options["disturbance"] == "large"
            for j in keys(sn_sol_bus_data)
                row[Symbol("E_I_Bus_$j")] = E_I_bus[j]
                row[Symbol("Delta_P_Bus_$j")] = delta_p_bus[j]
                row[Symbol("E_I_Min_Bus_$j")] = E_I_min_bus[j]
            end
        end

        if options["area"] == "true" && options["disturbance"] == "large"
            for r in areas
                row[Symbol("E_I_Area_$r")] = E_I_Area[r]
                row[Symbol("Delta_P_Area_$r")] = delta_p_area[r]
                row[Symbol("E_I_Min_Area_$r")] = E_I_min_area[r]
            end
        end

        df_row = DataFrame(row)

        results_filename = "results\\single_network_results\\results_$case_name.csv"
        if isfile(results_filename)
            existing_df = CSV.read(results_filename, DataFrame)
            for col_name in names(df_row)
                if !hasproperty(existing_df, col_name)
                    existing_df[!, col_name] = fill(missing, nrow(existing_df))
                end
            end
            combined_df = vcat(existing_df, df_row, cols=:union)
            combined_df = select(combined_df, :Options, Not(:Options))     
            CSV.write(results_filename, combined_df)
        else
            df_row = select(df_row, :Options, Not(:Options))    
            CSV.write(results_filename, df_row)
        end
    end

    save_results_as_csv(df, options, results)
end





















results_filename = "results\\results_$case_name.csv"
if isfile(results_filename)
    existing_df = CSV.read(results_filename, DataFrame)
    matching_rows = filter(row -> row.Options == options_name && row.rocof == rocof && row.Delta_P == delta_P, existing_df)

    if isempty(matching_rows)
        combined_df = vcat(existing_df, df)
    else
        for row in eachrow(matching_rows)
            for (col, value) in pairs(df[1, :])
                row[col] = value
            end
        end
        combined_df = existing_df
    end

    CSV.write(results_filename, combined_df)
else
    CSV.write(results_filename, df)
end



plot_filename = "results\\results_case2.csv"
function plot_results(plot_filename; plot_size=(800, 600), plot_type=:bar, multiple_plots=false)
    df = CSV.read(plot_filename, DataFrame)
    cases = 1:size(df, 1)  # Anpassen nach Bedarf für spezifische Fälle
    gen_ids = filter(x -> startswith(string(x), "PG_"), names(df))
    colors = [:skyblue, :lightgreen, :violet]

    if multiple_plots
        p = plot(layout=(length(gen_ids), 1), size=plot_size)
    else
        p = plot(size=plot_size)
    end

    for (i, gen_id) in enumerate(gen_ids)
        if multiple_plots
            pi = p[i]
        else
            pi = p
        end
        for (j, fall) in enumerate(cases)
            fall_df = df[fall, :] # Dies erhält eine Zeile als DataFrameRow
            wert = fall_df[gen_id] # Zugriff auf den Wert der Spalte `gen_id` für die aktuelle Zeile
            if plot_type == :bar
                bar!(pi, [fall], [wert], label="Fall $fall", color=colors[mod1(j, length(colors))])
            elseif plot_type == :line
                plot!(pi, [fall], [wert], label="Fall $fall", color=colors[mod1(j, length(colors))])
            end
        end
        display(pi)
        if multiple_plots
            title!(p[i], string(gen_id))
        end
    end

    if !multiple_plots
        legend=:outerright
    end

    return p
end

plot_results(plot_filename, plot_size=(800, 600), plot_type=:bar, multiple_plots=true)





Pg_values = [
    [3.0093, 3.0093, 0.0000], # Gen 1 Werte für jeden Fall
    [0.0000, 0.0000, 3.0093], # Gen 2 Werte für jeden Fall
    [7.0000, 7.0000, 7.0000], # Gen 3 Werte für jeden Fall
    [0.0000, 0.0000, 0.0000]  # Gen 4 Werte für jeden Fall
]

Cost_values = [
    [2407.43, 2407.43, 0.0000], # Kosten für Gen 1 in jedem Fall
    [0.0000, 0.0000, 3009.28],  # Kosten für Gen 2 in jedem Fall
    [3500.00, 3500.00, 3500.00],# Kosten für Gen 3 in jedem Fall
    [0.0000, 0.0000, 0.0000]    # Kosten für Gen 4 in jedem Fall
]

H_sys_values = [0.000, 5.19, 6.75]
H_min_values = [0.000, 5.0, 5.0]
colors = [:skyblue, :lightgreen, :violet]


cases = ["unit commetment", "uc with H_sys", "uc with H_area"]
gen_ids = ["Gen 1", "Gen 2", "Gen 3", "Gen 4"]

p1 = plot(size=(800, 600), title="active power for each generator and case")
p2 = plot(size=(800, 600), title="cost for each generator and case")
p3 = plot(size=(800, 600), title="H_sys and H_min values for each case")


num_generators = length(Pg_values) 
bar_width = 0.2
spacing = 0.1 
group_spacing = bar_width + spacing 

positions = [(i-1)*(3*group_spacing) + j*group_spacing for i=1:num_generators, j=1:3]

xlabel!(p1, "Generator")
ylabel!(p1, "active power [MW]")
xlabel!(p2, "Generator")
ylabel!(p2, "cost in [€]")

tick_positions = [(i-1)*(3*group_spacing) + 2*group_spacing for i=1:num_generators]

xticks!(p1, tick_positions, string.(1:num_generators))
xticks!(p2, tick_positions, string.(1:num_generators))

H_sys_positions = [(i-1)*group_spacing + 1.5bar_width for i=1:length(cases)]
H_min_positions = [(i-1)*group_spacing + 2.5*bar_width for i=1:length(cases)]

H_sys_values = [H_sys_values[case_index] for case_index in 1:length(cases)]
H_min_values = [H_min_values[case_index] for case_index in 1:length(cases)]

for (gen_index, pg_values_gen) in enumerate(Pg_values)
    for (case_index, pg_value) in enumerate(pg_values_gen)
        x_position = positions[gen_index, case_index]
        label = gen_index == 1 ? cases[case_index] : nothing
        bar!(p1, [x_position], [pg_value], label=label, bar_width=bar_width, color=colors[case_index])
    end
end

for (gen_index, cost_values_gen) in enumerate(Cost_values)
    for (case_index, cost_value) in enumerate(cost_values_gen)
        x_position = positions[gen_index, case_index]
        label = gen_index == 1 ? cases[case_index] : nothing
        bar!(p2, [x_position], [cost_value], label=false, bar_width=bar_width, color=colors[case_index])
    end
end
for case_index in 1:length(cases)
    x_position = case_index
    bar!(p3, [x_position - 0.2], [H_sys_values[case_index]], label=false, bar_width=0.4, color=colors[case_index])
    bar!(p3, [x_position + 0.2], [H_min_values[case_index]], label=false, bar_width=0.4, color=colors[case_index])
    annotate!(p3, [(x_position - 0.2, -0.5, text("H_sys", 8, :center, :top)), (x_position + 0.2, -0.5, text("H_min", 8, :center, :top))])
end

xticks!(p3, 1:length(cases), [" " for _ in 1:length(cases)])

#bar!(p3, H_sys_positions, H_sys_values, label="H_sys", bar_width=0.5*bar_width)
#bar!(p3, H_area_positions, H_area_values, label="H_min", bar_width=0.5*bar_width)

xlabel!(p3, "case")
ylabel!(p3, "H value in [s]")
plot(p1, p2, p3, layout=(3, 1), size=(800, 600))

=#