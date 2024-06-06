# Plot Skript

using JLD2
using CSV
using DataFrames
using Colors
using Measures
using Statistics
using PlotlyJS




# Load Data

const BUS_SYSTEM = "37"
const RELAX_INTEGRALITY = false
const PATH = "results\\results_bus_$(BUS_SYSTEM)\\relax_$(RELAX_INTEGRALITY)"
const CARRIERS_FILE = ".\\test\\data\\matpower\\multi_nw\\bus_$(BUS_SYSTEM)\\pypsa-eur_multinetwork_$(BUS_SYSTEM)\\carrier_indices.csv"


function load_all_dataframes(path)
    files = filter(f -> endswith(f, ".jld2"), readdir(path))
    all_data = Dict()
    for file in files
        try
            basename = split(file, '\\')[end]
            case_label = join(split(split(basename, '.')[1], '_')[1:2], "_")
            re_inj_int = split(split(basename, '.')[1], '_')[end]
            data = JLD2.load(joinpath(path, file))
            all_data["$(case_label)_$(re_inj_int)"] = data
        catch e
            println("Fehler beim Laden von $file: $e")
        end
    end

    return all_data
end

all_data = load_all_dataframes(PATH)

function load_cases(all_data)
    cases_constraint = Dict{String, Dict{Any, Any}}(
            "case_1" => Dict("constraint"=> "Unit Commitment"),
            "case_2" => Dict("constraint"=> "System Inertia"),
            "case_3" => Dict("constraint"=> "SD Weighted Load"),
            "case_4" => Dict("constraint"=> "SD Weighted Equal"),
            "case_5" => Dict("constraint"=> "Region Islanding 25%"),
            "case_6" => Dict("constraint"=> "Region Islanding 50%"),
            "case_7" => Dict("constraint"=> "Region Islanding 75%"),
            "case_8" => Dict("constraint"=> "Region Islanding 100%"),)
    filtered_cases = Dict{String, Dict{Any, Any}}()
    for case_label in keys(all_data)
        case_data = all_data[case_label]
        case_number = split(case_label, '_')[2]
        case_data["RE"] = string(split(case_label, '_')[3], "%")
        if haskey(cases_constraint, "case_$case_number")
            case_data["constraint"] = cases_constraint["case_$case_number"]["constraint"]
            case_data["number"] = case_number
            filtered_cases[case_label] = case_data
        end
    end
    return filtered_cases
end

cases = load_cases(all_data)

function build_gen_dataframe(all_data, cases)
    gen_df = DataFrame()
    for c in keys(cases)
        gen_c = all_data[c]["results"]["gen"]
        gen_c_df = DataFrame()
        for k in keys(gen_c)
            gen_c[k] = rename(gen_c[k], :value => k)
            if ncol(gen_c_df) == 0
                gen_c_df = gen_c[k]
                gen_c_df[!, "constraint"] .= cases[c]["constraint"]
                gen_c_df[!, "RE"] .= cases[c]["RE"]
                gen_c_df[!, "number"] .= cases[c]["number"]

            else
                gen_c_df = innerjoin(gen_c_df, gen_c[k], on= :id)
            end
        end
        gen_df = vcat(gen_df, gen_c_df)
    end
    return gen_df
end


function calculate_weighted_pg_sums(all_data, cases)
    pg_sums = DataFrame(id = Any[], constraint = String[], RE = String[], pg_sum = Float64[])
    for (case_label, case_data) in all_data
        if haskey(case_data["results"], "gen_t") && haskey(case_data["results"], "weight")
            gen_t_data = case_data["results"]["gen_t"]["pg"]
            weights = case_data["results"]["weight"][!, :weight]
            for col in names(gen_t_data)[2:end]  
                pg_values = gen_t_data[:, col]
                weighted_sum = sum(skipmissing(pg_values .* weights))
                constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]
                id = parse(Int, col)
                is_string = string(id)
                push!(pg_sums, (id = is_string, constraint = constraint, RE = RE, pg_sum = weighted_sum))
            end
        end
    end
    return pg_sums
end


function calculate_weighted_cost_sums(all_data, cases)
    cost_sums = DataFrame(id = Any[], constraint = String[], RE = String[], 
                          sd_cost_sum = Float64[], su_cost_sum = Float64[], pg_cost_sum = Float64[])
    for (case_label, case_data) in all_data
        if haskey(case_data["results"], "gen_t") && haskey(case_data["results"], "weight")
            shutdown_cost_data = case_data["results"]["gen_t"]["shutdown_cost"]
            startup_cost_data = case_data["results"]["gen_t"]["startup_cost"]
            operational_cost_data = case_data["results"]["gen_t"]["pg_cost"]
            weights = case_data["results"]["weight"][!, :weight]
            for col in names(shutdown_cost_data)[2:end]
                sd_cost_values = shutdown_cost_data[:, col]
                su_cost_values = startup_cost_data[:, col]
                operat_cost_values = operational_cost_data[:, col]
                sd_cost_sum = sum(skipmissing(sd_cost_values .* weights))
                su_cost_sum = sum(skipmissing(su_cost_values .* weights))
                pg_cost_sum = sum(skipmissing(operat_cost_values .* weights))
                constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]
                id = parse(Int, col)
                is_string = string(id)
                push!(cost_sums, (id = is_string, constraint = constraint, RE = RE, 
                                  sd_cost_sum = sd_cost_sum, su_cost_sum = su_cost_sum, pg_cost_sum = pg_cost_sum))
            end
        end
    end
    return cost_sums
end


function create_results_gen_df()
    all_data = load_all_dataframes(PATH)
    cases = load_cases(all_data)

    gen_df = build_gen_dataframe(all_data, cases)
    carrier_names = CSV.read(CARRIERS_FILE, DataFrame, header=0)

    gen_df[!, :carrier_names] = [carrier_names[v + 1, "Column1"] for v in gen_df[:, :carrier]]
    gen_df[!, :pExp] = (gen_df[!, :nE] - gen_df[!, :n0]) .* gen_df[!, :P_b_nom] / 1e3
    gen_df[!, :pIns] = gen_df[!, :n0] .* gen_df[!, :P_b_nom] / 1e3
    pg_sums = calculate_weighted_pg_sums(all_data, cases)
    cost_sums = calculate_weighted_cost_sums(all_data, cases)
    gen_df = leftjoin(gen_df, cost_sums, on = [:id, :constraint, :RE])

    gen_df = leftjoin(gen_df, pg_sums, on = [:id, :constraint, :RE])
    gen_df[!, :carrier_group] = map(gen_df[!, :carrier]) do carrier
        if carrier == 6
            "PV"
        elseif carrier == 4
            "Onwind"
        elseif carrier == 2
            "Offwind"
        elseif carrier == 14
            "PHS_pump"
        elseif carrier in [0, 1]
            "Gas"
        else 
            "Others"
        end
    end    
    carrier_order = Dict("PV" => 1, "Onwind" => 2, "Offwind" => 3, "PHS_pump" => 4, "Gas" => 5, "Others" => 6)
    gen_df[!, :carrier_order] = map(v -> carrier_order[v], gen_df[!, :carrier_group])


    total_load = 3.21353e7

    gen_df[!, :percentage_pExp] = gen_df[!, :pExp] ./ gen_df[!, :pIns]
    g_gen_df_pExp = groupby(gen_df, [:constraint, :RE, :carrier_group])
    pExp_sums_df = combine(g_gen_df_pExp, :percentage_pExp => sum => :percentage_pExp_group)
    merged_gen_df = leftjoin(gen_df, pExp_sums_df, on = [:constraint, :RE, :carrier_group])

    g_gen_df = groupby(merged_gen_df, [:constraint, :RE, :carrier_group])
    pg_sums_df = combine(g_gen_df, :pg_sum => sum => :pg_sum_group)
    pg_sums_df[!, :percentage_pg_group] = pg_sums_df[!, :pg_sum_group] / total_load
    merged_gen_df = leftjoin(merged_gen_df, pg_sums_df, on = [:constraint, :RE, :carrier_group])

    return merged_gen_df
end

results_gen_df = create_results_gen_df()

