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


is_multi_network = true
case_name = "mpc_multinetwork_10"
bus_system = case_name[end-1:end]
results_file = "results\\multi_network_results\\results_bus_$bus_system\\results_v_$case_name.jld2"

all_data = JLD2.load(results_file)["results_v"]
mn_data = all_data["data"]["mn_data"]
test_cases = sort(collect(keys(all_data["cases"])))

function add_calculated_values_to_df(test_cases, mn_data)

    for test_case in test_cases
        println("Adding calculated values to DataFrame for test case: $test_case")

        for nw in sort(collect(keys(all_data["cases"][test_case]["Results"]["solution"]["nw"])))

            case_data = all_data["cases"][test_case]
            results = Dict()
            options = Dict()
            E_I_weighted_area = Dict()
            P_load_area = Dict()
            P_gen_area = Dict()
            delta_p_area = Dict()
            E_I_min_area = Dict()
            P_gen_bus = Dict()
            P_load_bus = Dict()
            delta_p_bus = Dict()
            E_I_min_bus = Dict()
            E_I_bus = Dict()
            E_I_weighted_area = Dict()
            E_I_weighted_area = Dict()
            E_I_area = Dict()
            E_I_min_area = Dict()
    

            if haskey(case_data, "Error") && case_data["Error"] == "None"
                results[test_case] = all_data["cases"][test_case]["Results"]
                options[test_case] =  all_data["cases"][test_case]["Options"]
                gen_data = mn_data["nw"][nw]["gen"]
                bus_data = mn_data["nw"][nw]["bus"]
                load_data = mn_data["nw"][nw]["load"]
                f0 = 50
                alpha = options[test_case]["f"]["alpha_factor"]
                rocof = options[test_case]["f"]["rocof"]

                sol_gen_data = all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["gen"]
                sol_bus_data = all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["bus"]

                for (gen_id, gen) in gen_data
                    H = gen_data[gen_id]["H"]
                    pmin = gen_data[gen_id]["pmin"]
                    pmax = gen_data[gen_id]["pmax"]
                    pmin = abs(pmin)
                    pmax = abs(pmax)
                    pv = [pmin, pmax]
                    pmin = minimum(pv)
                    pmax = maximum(pv)
                    gen_bus = gen_data[gen_id]["gen_bus"]

                    if !haskey(sol_gen_data, gen_id)
                        sol_gen_data[gen_id] = Dict()
                    end

                    sol_gen_data[gen_id]["H"] = H
                    sol_gen_data[gen_id]["pmin"] = pmin
                    sol_gen_data[gen_id]["pmax"] = pmax
                    
                    sol_gen_data[gen_id]["gen_bus"] = gen_bus
                end
        
                for (bus_id, bus) in bus_data
                    area = bus_data[bus_id]["area"]
                    sol_bus_data[bus_id]["area"] = area
                end
        
                for i in keys(sol_gen_data)
                    if sol_gen_data[i]["gen_status"] < 1e-2
                        sol_gen_data[i]["gen_status"] = 0
                    end
                    if sol_gen_data[i]["pg"] < 1e-2
                        sol_gen_data[i]["gen_status"] = 0
                    end
                end
            end
            
            delta_P = mn_data["delta_P"]    
            E_I_min = mn_data["E_I_min"]
            all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["Delta_P"] = delta_P
            all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_min"] = E_I_min

            if options[test_case]["f"]["system"] == "true"
    
                E_I_sys = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in keys(sol_gen_data))
                all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_sys"] = E_I_sys
            end
    
            if options[test_case]["f"]["weighted_area"] == "load" && options[test_case]["f"]["disturbance"] == "small"
                    
                areas = unique([sol_bus_data[j]["area"] for j in keys(sol_bus_data)])
                for area in areas
                    gens_in_area = [i for i in keys(sol_gen_data) if sol_bus_data[string(sol_gen_data[i]["gen_bus"])]["area"] == area]
                    E_I_weighted_area[area] = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
                    
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Weighted_Area_$area"] = E_I_weighted_area[area]
                end
    
            end
    
            if options[test_case]["f"]["weighted_area"] == "equal" && options[test_case]["f"]["disturbance"] == "small"
                    
                areas = unique([sol_bus_data[j]["area"] for j in keys(sol_bus_data)])
                for area in areas
                    gens_in_area = [i for i in keys(sol_gen_data) if sol_bus_data[string(sol_gen_data[i]["gen_bus"])]["area"] == area]
                    E_I_weighted_area[area] = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)

                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Weighted_Area_$area"] = E_I_weighted_area[area]
                end
    
            end
    
            if options[test_case]["f"]["weighted_area"] == "none" && options[test_case]["f"]["disturbance"] == "small"
                    
                areas = unique([sol_bus_data[j]["area"] for j in keys(sol_bus_data)])
                for area in areas
                    E_I_weighted_area[area] = nothing
                    
                    # add the calculated values to the DataFrame
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Weighted_Area_$area"] = E_I_weighted_area[area]
                end
    
            end
    
            if options[test_case]["f"]["area"] == "true" && options[test_case]["f"]["disturbance"] == "large"
    
                areas = unique([sol_bus_data[j]["area"] for j in keys(sol_bus_data)])
                for area in areas
    
                    gens_in_area = [i for i in keys(sol_gen_data) if sol_bus_data[string(sol_gen_data[i]["gen_bus"])]["area"] == area]
                    buses_in_area = [i for i in keys(sol_bus_data) if sol_bus_data[i]["area"] == area]
                    P_load_area[area] = sum(load_data[i]["pd"] for i in keys(load_data) if load_data[i]["load_bus"] in buses_in_area; init=0)
                    P_gen_area[area] = sum(sol_gen_data[i]["pg"] for i in gens_in_area; init=0)
                    delta_p_area[area] = abs(P_gen_area[area] - P_load_area[area])
                    E_I_min_area[area] = (delta_p_area[area] * f0 / rocof * 2)
                    E_I_area[area] = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)
                    
                    # add the calculated values to the DataFrame
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Area_$area"] = E_I_area[area]
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["Delta_P_Area_$area"] = delta_p_area[area]
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Min_Area_$area"] = E_I_min_area[area]
                end
    
            end
    
    
            if options[test_case]["f"]["bus"] == "true" && options[test_case]["f"]["disturbance"] == "large"
    
                for j in keys(sol_bus_data)
    
                    gens_at_bus = [i for i in keys(sol_gen_data) if string(sol_gen_data[i]["gen_bus"]) == j]
                    P_gen_bus[j] = sum(sol_gen_data[i]["pg"] for i in gens_at_bus; init=0)
                    P_load_bus[j] = sum(load_data[i]["pd"] for i in keys(load_data) if string(load_data[i]["load_bus"]) == j; init=0)
                    delta_p_bus[j] = abs(P_gen_bus[j] - P_load_bus[j])
                    H_min_bus[j] = (delta_p_bus[j] * f0) / (P_load_bus[j] * 2 * rocof)
                    E_I_bus[j] = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_at_bus; init=0)
                    
                    # add the calculated values to the DataFrame
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Bus_$j"] = E_I_bus[j]
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["Delta_P_Bus_$j"] = delta_p_bus[j]
                    all_data["cases"][test_case]["Results"]["solution"]["nw"][nw]["E_I_Min_Bus_$j"] = E_I_min_bus[j]
                end
            end
        end
    end
    return all_data
end

add_calculated_values_to_df(test_cases, mn_data)

# create a PowerModelsDataFrame from the results for all test cases

function create_pmdfs(all_data, test_cases)

    pmdfs = Dict{String, Dict{String, PowerModelsDataFrame}}()
    for (test_case, case_data) in test_cases
        pmdfs[test_case] = Dict{String, PowerModelsDataFrame}()

        for (nw_id, nw_data) in all_data["data"]["mn_data"]["nw"]

            data_for_pmdf = Dict{String, Any}()
            data_for_pmdf["gen"] =  all_data["cases"][test_cases]["Results"]["solution"]["nw"][nw_id]["gen"]
            data_for_pmdf["bus"] = all_data["cases"][test_cases]["Results"]["solution"]["nw"][nw_id]["bus"]
            data_for_pmdf["branch"] = all_data["cases"][test_cases]["Results"]["solution"]["nw"][nw_id]["branch"]

            pmdfs[test_case][nw_id] = PowerModelsDataFrame(data_for_pmdf)
        end
    end
    return pmdfs
end

#=
function create_pmdf_for_solution(test_cases, all_data)
    pmdf = PowerModelsDataFrame(Dict{String, Any}())
    for test_case in test_cases
        println("Creating PowerModelsDataFrame for test case: $test_case")
        solution_data = all_data["cases"][test_case]["Results"]["solution"]
        solution_data_string = Dict{String, Any}(solution_data)
        pmdf = PowerModelsDataFrame(solution_data_string)
        all_data["cases"][test_case]["Results"]["solution"]["PowerModelsDataFrame"] = pmdf
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