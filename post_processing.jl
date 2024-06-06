# Script to post-process the results of the optimization problem

using Plots
using Plots: text
using JLD2
using CSV
using DataFrames
using Colors
using StatsPlots
using Measures
using GraphRecipes
using GraphPlot
using Graphs
using GraphPlot
using LightGraphs
using StatsPlots
using Polynomials
using Statistics
using LaTeXStrings
using CategoricalArrays
# using PowerPlots
# using Glob
# using Random


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


function build_branch_dataframe(all_data, cases)
    branch_df = DataFrame()
    for c in keys(cases)
        branch_c = all_data[c]["results"]["branch"]
        branch_c_df = DataFrame()
        for k in keys(branch_c)
            branch_c[k] = rename(branch_c[k], :value => k)

            if ncol(branch_c_df) == 0
                branch_c_df = branch_c[k]
                branch_c_df[!, "constraint"] .= cases[c]["constraint"]
                branch_c_df[!, "RE"] .= cases[c]["RE"]
                branch_c_df[!, "number"] .= cases[c]["number"]
            else
                branch_c_df = innerjoin(branch_c_df, branch_c[k], on= :id)
            end
        end
        branch_df = vcat(branch_df, branch_c_df)
    end
    return branch_df
end



function calculate_pf_pt_sums(all_data, cases)
    pt_sums = DataFrame(id = Any[], constraint = String[], RE = String[], pt_sum = Float64[])
    pf_sums = DataFrame(id = Any[], constraint = String[], RE = String[], pf_sum = Float64[])

    for (case_label, case_data) in all_data
        pf_data = case_data["results"]["branch_t"]["pf"]
        pt_data = case_data["results"]["branch_t"]["pt"]

        for col in names(pf_data)[2:end]
            pf_values = pf_data[:, col]
            weighted_sum_pf = sum(skipmissing(pf_values))
            constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]
            id = parse(Int64, col)
            id_str = string(id)
            push!(pf_sums, (id = id_str, constraint = constraint, RE = RE, pf_sum = weighted_sum_pf))
        end

        for col in names(pt_data)[2:end]
            pt_values = pt_data[:, col]
            weighted_sum_pt = sum(skipmissing(pt_values))
            constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]
            id = parse(Int64, col)
            id_str = string(id)
            push!(pt_sums, (id = id_str, constraint = constraint, RE = RE, pt_sum = weighted_sum_pt))
        end
    end
    
    return pt_sums, pf_sums
end
pf_sums_test = calculate_pf_pt_sums(all_data, cases)[2]


function calculate_weighted_pg_sums(all_data, cases)
    pg_sums = DataFrame(id = Any[], constraint = String[], RE = String[], pg_sum = Float64[])
    for (case_label, case_data) in all_data
        if haskey(case_data["results"], "gen_t") && haskey(case_data["results"], "weight")
            gen_t_data = case_data["results"]["gen_t"]["pg"]
            weights = case_data["results"]["weight"][!, :weight]
            for col in names(gen_t_data)[2:end]  
                pg_values = gen_t_data[:, col]
                pg_weighted_sum = sum(skipmissing(pg_values .* weights))
                constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]
                id = parse(Int, col)
                is_string = string(id)
                push!(pg_sums, (id = is_string, constraint = constraint, RE = RE, pg_sum = pg_weighted_sum))
            end
        end
    end
    return pg_sums
end

pg_sums_test = calculate_weighted_pg_sums(all_data, cases)


function calculate_weighted_cost_sums(all_data, cases)
    cost_sums = DataFrame(id = Any[], constraint = String[], RE = String[], invest_cost_sum = Float64[],
                          sd_cost_sum = Float64[], su_cost_sum = Float64[], pg_cost_sum = Float64[])
    for (case_label, case_data) in all_data
        if haskey(case_data["results"], "gen_t") && haskey(case_data["results"], "weight")

            invest_cost_data = case_data["results"]["gen_t"]["investment_cost"]
            shutdown_cost_data = case_data["results"]["gen_t"]["shutdown_cost"]
            startup_cost_data = case_data["results"]["gen_t"]["startup_cost"]
            pg_cost_data = case_data["results"]["gen_t"]["pg_cost"]



            for col in names(shutdown_cost_data)[2:end]
                inv_cost_values = invest_cost_data[:, col]
                sd_cost_values = shutdown_cost_data[:, col]
                su_cost_values = startup_cost_data[:, col]
                pg_cost_values = pg_cost_data[:, col]

                invest_cost_sum = sum(skipmissing(inv_cost_values))
                sd_cost_sum = sum(skipmissing(sd_cost_values))
                su_cost_sum = sum(skipmissing(su_cost_values))
                pg_cost_sum = sum(skipmissing(pg_cost_values))

                constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]
                id = parse(Int, col)
                push!(cost_sums, (id = id, constraint = constraint, RE = RE,  invest_cost_sum = invest_cost_sum,
                                  sd_cost_sum = sd_cost_sum, su_cost_sum = su_cost_sum, pg_cost_sum = pg_cost_sum))
            end
        end
    end
    return cost_sums
end

test_cost_sums = calculate_weighted_cost_sums(all_data, cases)

function calculate_total_costs(all_data, cases)

    total_costs = DataFrame(constraint = String[], RE = String[], number = Int[], total_cost = Float64[])

    cost_sums = calculate_weighted_cost_sums(all_data, cases)

    for (case_label, case_data) in all_data
        if haskey(case_data["results"], "general_var")
            total_cost = case_data["results"]["general_var"][1, :objective]
            constraint, RE = cases[case_label]["constraint"], cases[case_label]["RE"]

            number = parse(Int, cases[case_label]["number"])  
            push!(total_costs, (constraint = constraint, RE = RE, number = number, total_cost = total_cost))
        end
    end

    cost_sums_grouped = combine(groupby(cost_sums, [:constraint, :RE]), 
                                :invest_cost_sum => sum => :invest_cost_sum_case,
                                :sd_cost_sum => sum => :sd_cost_sum_case, 
                                :su_cost_sum => sum => :su_cost_sum_case, 
                                :pg_cost_sum => sum => :pg_cost_sum_case)
    total_costs = innerjoin(total_costs, cost_sums_grouped, on = [:constraint, :RE])

    total_costs[!, :invest_cost_percentage] = total_costs[!, :invest_cost_sum_case] ./ total_costs[!, :total_cost] * 100
    total_costs[!, :sd_cost_percentage] = total_costs[!, :sd_cost_sum_case] ./ total_costs[!, :total_cost] * 100
    total_costs[!, :su_cost_percentage] = total_costs[!, :su_cost_sum_case] ./ total_costs[!, :total_cost] * 100
    total_costs[!, :pg_cost_percentage] = total_costs[!, :pg_cost_sum_case] ./ total_costs[!, :total_cost] * 100

    return total_costs
end


test_total_costs = calculate_total_costs(all_data, cases)

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
        elseif carrier in [0, 1]
            "Gas"
        else 
            "Others"
        end
    end    
    carrier_order = Dict("PV" => 1, "Onwind" => 2, "Offwind" => 3, "Gas" => 5, "Others" => 6)
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

    return merged_gen_df, pg_sums_df
end

results_gen_df, pg_sums_df = create_results_gen_df()


function create_results_branch_df()
    all_data = load_all_dataframes(PATH)
    cases = load_cases(all_data)

    branch_df = build_branch_dataframe(all_data, cases)
    g_branch_df = groupby(branch_df, [:constraint, :RE])
    pt_sums, pf_sums = calculate_pf_pt_sums(all_data, cases)
    merged_branch_df = leftjoin(branch_df, pt_sums, on = [:id, :constraint, :RE])
    merged_branch_df = leftjoin(merged_branch_df, pf_sums, on = [:id, :constraint, :RE])

    return merged_branch_df
end

results_branch_df = create_results_branch_df()

###################################################################### plot pg stacked ######################################################################


function plot_stacked_generation(results_gen_df)
    plots = []
    re_levels = sort(unique(results_gen_df.RE))

    for (i, re) in enumerate(re_levels)
        re_df = filter(row -> row.RE == re, results_gen_df)
        if isempty(re_df)
            println("Keine Daten für RE = $re gefunden!")
        else
            re_df[!, :percentage_pg_group] = re_df[!, :percentage_pg_group] .* 100
            re_df[!, :carrier_group] = categorical(re_df[!, :carrier_group])
            
            levels!(re_df[!, :carrier_group], ["Others", "Gas", "Offwind", "Onwind", "PV"] )  # ["PV", "Onwind", "Offwind", "Gas", "Others"]
            legend_setting = i == 1 ? :outertopright : :none

            sort!(re_df, :number)
            subplot = groupedbar(
                re_df[!, :constraint],
                re_df[!, :percentage_pg_group],
                group=re_df[!, :carrier_group],
                bar_width=0.95,
                bar_position=:stack,
                xlabel=i == 1 ? "" : "Constraint",
                xaxis = :flip,
                ylabel="Energy (%)",
                title="RE = $re",
                legend=legend_setting,
                titlefontsize=14,
                xguidefontsize=12,
                yguidefontsize=12,
                xtickfontsize=10,
                ytickfontsize=10,
                legendfontsize=12,
                xrotation=20
            )
            
            push!(plots, subplot)
        end
    end

    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 800), margin=6mm)
    return final_plot
end


plot_stacked_pg = plot_stacked_generation(results_gen_df)
savefig(plot_stacked_pg, "$PATH/plot_stacked_generation.pdf")
display(plot_stacked_pg)





###################################################################### plot pExp ######################################################################






function create_df_for_plot_Exp()

    results_gen_df = create_results_gen_df();
    all_data = load_all_dataframes(PATH);
    carriers = CSV.read(CARRIERS_FILE, DataFrame, header=0)
    rows = NamedTuple{(:id, :constraint, :RE, :carrier, :Pmax_Pb_nom, :pg_Pb_nom, :nE, :carrier_names), 
                      Tuple{Any, String, String, Any, Float64, Float64, Any, String}}[]
    cases = load_cases(all_data)

    for (case_key, case_details) in cases
        gen_t_data = all_data[case_key]["results"]["gen_t"]
        pmax_df = gen_t_data["pmax"]
        pg_df = gen_t_data["pg"]
        println(case_key)

        for gen_id in unique(results_gen_df.id)
            carrier = results_gen_df[results_gen_df.id .== gen_id, :carrier][1]
            Pb_nom = results_gen_df[results_gen_df.id .== gen_id, :P_b_nom][1]
            nE_value = results_gen_df[results_gen_df.id .== gen_id, :nE][1]
            carrier_name = results_gen_df[results_gen_df.id .== gen_id, :carrier_names][1]


            for nw in 1:size(pmax_df, 1)
                pmax_value = pmax_df[nw, gen_id]
                pg_value = pg_df[nw, gen_id]
                Pmax_Pb_nom = pmax_value * 100 / Pb_nom
                pg_Pb_nom = pg_value * 100 / (Pb_nom * nE_value)

                push!(rows, (id=gen_id, constraint=case_details["constraint"], RE=case_details["RE"], 
                             carrier=carrier, Pmax_Pb_nom=Pmax_Pb_nom, pg_Pb_nom=pg_Pb_nom, nE=nE_value, carrier_names=carrier_name))
            end
        end
    end

    df_Exp = DataFrame(rows)
    df_Exp = df_Exp[df_Exp.nE .!= 0, :] 
    df_Exp = df_Exp[in.(df_Exp.carrier_names, Ref(["solar", "onwind", "offwind"])), :] 


    return df_Exp
end

test_df_exp = create_df_for_plot_Exp()

function plot_carrier_data(carrier_name)
    df = create_df_for_plot_Exp()
    df_carrier = df[df[!, :carrier_names] .== carrier_name, :]

    melted_df = stack(df_carrier, [:Pmax_Pb_nom, :pg_Pb_nom], [:constraint])
    rename!(melted_df, :variable => :Measurement, :value => :Value)
    grouped = groupby(melted_df, [:constraint, :Measurement])

    p = plot(title=" ",
             ylabel=" ", xlabel="cases", legend=:outertopright)

    for g in grouped
    constraint_label = unique(g[!, :constraint])[1]
    measurement_label = unique(g[!, :Measurement])[1]
    values = g[!, :Value]
    
    violin!(p, [String(constraint_label)], values, label=String(measurement_label))
end

    display(p)
end

plot_carrier_data("solar")


function plot_carrier_data(carrier_name)
    df = create_df_for_plot_Exp()
    df_carrier = df[df.carrier_names .== carrier_name, :]

    melted_df = stack(df_carrier, [:Pmax_Pb_nom, :pg_Pb_nom], [:constraint])
    rename!(melted_df, :variable => :Measurement, :value => :Value)

    grouped = groupby(melted_df, [:constraint, :Measurement])

    p = plot(title=" ",
             ylabel=" ", xlabel="cases", legend=:topright,
             xrotation=20, size= (1200, 800), margin=5mm,
             titlefontsize=14, labelfontsize=12,
             xtickfontsize=10, legendfontsize=12)

    unique_measurements = Set()
    
    for g in grouped
        measurement_label = first(g).Measurement
        if !(measurement_label in unique_measurements)
            push!(unique_measurements, measurement_label)
            label = measurement_label
        else
            label = "" 
        end
        
        constraint_label = unique(g[!, :constraint])[1]
        values = g[!, :Value]
        
        violin!(p, [constraint_label * " " * label], values, label=label)
    end

    display(p)
end

plot_carrier_data("solar")


function plot_carrier_data(carrier_name)
    df = create_df_for_plot_Exp()
    df_carrier = df[df.carrier_names .== carrier_name, :]

    melted_df = stack(df_carrier, [:Pmax_Pb_nom, :pg_Pb_nom], [:constraint])
    rename!(melted_df, :variable => :Measurement, :value => :Value)

    grouped = groupby(melted_df, [:constraint])

    p = plot(title="Vergleich von Pmax/Pb_nom und Pg/Pb_nom",
             ylabel="Werte", xlabel="Fälle", legend=:outertopright,
             xrotation=45, size=(1200, 800), titlefontsize=14,
             labelfontsize=12, xtickfontsize=10, legendfontsize=12)

    color_map = Dict("Pmax_Pb_nom" => :blue, "pg_Pb_nom" => :red)
    labels_added = Dict("Pmax_Pb_nom" => false, "pg_Pb_nom" => false)

    for (i, g) in enumerate(grouped)
        for (j, measurement) in enumerate(["Pmax_Pb_nom", "pg_Pb_nom"])
            subdata = g[g.Measurement .== measurement, :]
            positions = [i + (j-1)*0.3]

            label = measurement 
            if !labels_added[measurement] 
            else 
                label =""
            labels_added[measurement] = true

            violin!(p, positions, subdata.Value, label=label, color=color_map[measurement])
        end
    end

    xticks!(p, [(i + 0.15) for i in 1:length(grouped)], [String(g[1, :constraint]) for g in grouped])

    display(p)
end


plot_carrier_data("solar")



###################################################################### plot satcked cost ######################################################################

function plot_stacked_costs()
    costs_df = calculate_total_costs(all_data, cases)
    plots = []
    re_levels = sort(unique(costs_df.RE))
    # desired_order = ["Unit Commitment", "System Inertia", "SD Weighted Load", "Region Islanding 50%", "Region Islanding 75%", "Region Islanding 25%", "SD Weighted Equal", "Region Islanding 100%"]

    for (i, re) in enumerate(re_levels)
        re_df = filter(row -> row.RE == re, costs_df)

        if isempty(re_df)
            println("No data for RE = $re ")
        else

            re_df = sort!(re_df, :number)

            constraints = re_df.constraint
            data = vcat(
                re_df.sd_cost_percentage,
                re_df.su_cost_percentage,
                re_df.pg_cost_percentage,
                re_df.invest_cost_percentage
            )

            groups = categorical(repeat(["Shutdown", "Startup", "Power", "Investment"], inner = nrow(re_df)))
            levels!(groups, ["Shutdown", "Startup", "Power", "Investment"])
            repeated_constraints = repeat(constraints, outer = 4)
            legend_setting = i == 1 ? :outertopright : :none

            subplot = groupedbar(
                repeated_constraints,
                data,
                group = categorical(groups),
                xlabel = i == 1 ? "" : "Cases",
                xaxis = :flip,
                ylabel = "Cost (%)",
                title = "RE = $re",
                bar_width = 0.95,
                bar_position = :stack,
                legend = legend_setting,
                titlefontsize = 14,
                xguidefontsize = 12,
                yguidefontsize = 12,
                xtickfontsize = 10,
                ytickfontsize = 10,
                legendfontsize = 12,
                xrotation = 20
            )
            
            push!(plots, subplot)
        end
    end

    final_plot = plot(plots..., layout = (length(plots), 1), size = (1200, 800), margin = 5mm)
    return final_plot
end

plot_stack_costs = plot_stacked_costs()
savefig(plot_stack_costs, "$PATH/plot_stacked_costs.pdf")
display(plot_stack_costs)


###################################################################### plot cost ######################################################################

function plot_cost_relative()
    all_data = load_all_dataframes(PATH)
    cases = load_cases(all_data)
    plots = []

    RE_levels = sort(unique([v["RE"] for v in values(cases)]))
    
    for (index, RE) in enumerate(RE_levels)
        filtered_cases = [k for (k, v) in cases if v["RE"] == RE]
        RE = replace(RE, r"\%" => "")

        base_case = "case_1_$RE"

        reference_cost = all_data[base_case]["results"]["general_var"].objective[1]
        objectives = []
        constraints = []
        numbers = []

        for case_label in filtered_cases
            if haskey(all_data, case_label) && hasproperty(all_data[case_label]["results"]["general_var"], :objective)
                objective = all_data[case_label]["results"]["general_var"].objective[1]
                percent_relative = (objective / reference_cost) * 100
                push!(constraints, cases[case_label]["constraint"])
                push!(objectives, percent_relative)
                push!(numbers, parse(Int, cases[case_label]["number"]))
            end
        end

        sorted_indices = sortperm(numbers)
        sorted_constraints = [constraints[i] for i in sorted_indices]
        sorted_objectives = [objectives[i] for i in sorted_indices]

        if !isempty(sorted_constraints)
            xlabel_text = index == 1 ? "" : "Cases"
            p = bar(sorted_constraints, sorted_objectives, title="RE = $RE%", ylabel="Reference-related costs(%)", xlabel=xlabel_text, legend=false, margin=5mm)
            plot!(p, legend=false, titlefontsize=14, xguidefontsize=12, yguidefontsize=12, xtickfontsize=10, ytickfontsize=10, xticks=(1:length(sorted_constraints), sorted_constraints), xrotation=20)
            push!(plots, p)
        end
    end

    if isempty(plots)
        println("No data available for plots")
    else
        final_plot = plot(plots..., layout = @layout([a{0.47h}; b{0.47h}]), size = (1200, 800))
        return final_plot
    end
end

plot_total_cost = plot_cost_relative()
savefig(plot_total_cost, "$PATH/plot_total_cost.pdf")
display(plot_total_cost)


###################################################################### plot time ######################################################################

function plot_time_relative()
    all_data = load_all_dataframes(PATH)
    cases = load_cases(all_data)
    plots = []

    RE_levels = sort(unique([v["RE"] for v in values(cases)]))
    
    for (index, RE) in enumerate(RE_levels)
        filtered_cases = [k for (k, v) in cases if v["RE"] == RE]
        RE = replace(RE, r"\%" => "")

        base_case = "case_1_$RE"

        if haskey(all_data, base_case) && hasproperty(all_data[base_case]["results"]["general_var"], :solve_time)
            reference_time = all_data[base_case]["results"]["general_var"].solve_time[1]
            times = []
            constraints = []
            numbers = []

            for case_label in filtered_cases
                if haskey(all_data, case_label) && hasproperty(all_data[case_label]["results"]["general_var"], :solve_time)
                    time = all_data[case_label]["results"]["general_var"].solve_time[1]
                    percent_relative = (time / reference_time) * 100
                    push!(constraints, cases[case_label]["constraint"])
                    push!(times, percent_relative)
                    push!(numbers, parse(Int, cases[case_label]["number"]))
                end
            end

            sorted_indices = sortperm(numbers)
            sorted_constraints = [constraints[i] for i in sorted_indices]
            sorted_times = [times[i] for i in sorted_indices]

            xlabel_text = index == 1 ? "" : "Cases"
            p = bar(sorted_constraints, sorted_times, title="RE = $RE%", ylabel="Reference-related time (%)", xlabel=xlabel_text, legend=false, margin=5mm)
            plot!(p, legend=false, titlefontsize=14, xguidefontsize=12, yguidefontsize=12, xtickfontsize=10, ytickfontsize=10, xticks=(1:length(sorted_constraints), sorted_constraints), xrotation=20)
            annotate!(p, [(sorted_constraints[1], sorted_times[1], text(string(round(reference_time, digits=0), " s"), 14, :black, :center, :bottom))])
            push!(plots, p)
        end
    end

    if isempty(plots)
        println("No data available for plots")
    else
        final_plot = plot(plots..., layout =@layout([a{0.47h}; b{0.47h}]), size = (1200, 800))
        return final_plot
    end
end

plot_total_time = plot_time_relative()
savefig(plot_total_time, "$PATH/plot_total_time.pdf")
display(plot_total_time)




################################################################### branch pf ###################################################################

function plot_normalized_pf_whisker(all_data, cases)
    plots = []
    RE_levels = sort(unique([v["RE"] for v in values(cases)]))

    for RE in RE_levels
        all_pf_values = []
        labels = []
        numbers = [] 

        filtered_cases = [k for (k, v) in cases if v["RE"] == RE]
        for case_label in filtered_cases

            branch_data = all_data[case_label]["results"]["branch_t"]["pf"]
            rate_a_data = all_data[case_label]["results"]["branch"]["rate_a"][:, 2]
            constraint = all_data[case_label]["constraint"]
            number = parse(Int, all_data[case_label]["number"])

            if !isempty(branch_data) && !isempty(rate_a_data)
                pf_values = vec(Matrix(branch_data[:, 2:end]))
                rate_a_repeated = repeat(rate_a_data, outer=size(branch_data, 1))
                normalized_pf = pf_values ./ rate_a_repeated

                push!(all_pf_values, normalized_pf)
                push!(labels, constraint)
                push!(numbers, number)
            else
                println("Daten fehlen für $case_label")
            end
        end
        
        sorted_indices = sortperm(numbers)
        sorted_labels = [labels[i] for i in sorted_indices]
        sorted_pf_values = [all_pf_values[i] for i in sorted_indices]

        p = boxplot(sorted_labels, sorted_pf_values, title="RE = $RE", legend=false, color=:lightblue,
                    whisker_color=:darkgrey, median_color=:black, alpha=0.6)
        
        plot!(p, titlefontsize=14, guidefontsize=12, tickfontsize=10, legendfontsize=10, ylabel="Branch utilization (%)", xrotation=20)

        push!(plots, p)
    end
    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 800), margin=5mm)
    return final_plot
end





function plot_normalized_pf_whisker(all_data, cases)
    plots = []
    RE_levels = sort(unique([v["RE"] for v in values(cases)]))

    for RE in RE_levels
        all_pf_values = []
        labels = []
        numbers = []

        filtered_cases = [k for (k, v) in cases if v["RE"] == RE]
        for case_label in filtered_cases
            branch_data = all_data[case_label]["results"]["branch_t"]["pf"]
            rate_a_data = all_data[case_label]["results"]["branch"]["rate_a"]
            
            rate_a_dict = Dict(string(id) => value for (id, value) in zip(rate_a_data[:, 1], rate_a_data[:, 2]))

            ids = names(branch_data)[2:end]
            rate_a_values = [rate_a_dict[string(id)] for id in ids]

            normalized_pf = [Vector(branch_data[row, 2:end]) ./ rate_a_values for row in 1:nrow(branch_data)]
            
            constraint = all_data[case_label]["constraint"]
            number = parse(Int, all_data[case_label]["number"])

            if !isempty(normalized_pf)
                push!(all_pf_values, vec(normalized_pf))
                push!(labels, constraint)
                push!(numbers, number)
            else
                println("Daten fehlen für $case_label")
            end
        end

        sorted_indices = sortperm(numbers)
        sorted_labels = [labels[i] for i in sorted_indices]
        sorted_pf_values = [all_pf_values[i] for i in sorted_indices]

        p = boxplot(sorted_labels, sorted_pf_values, title="RE = $RE", legend=false, color=:lightblue,
                    whisker_color=:darkgrey, median_color=:black, alpha=0.6, xaxis=:flip)
        
        plot!(p, titlefontsize=14, guidefontsize=12, tickfontsize=10, legendfontsize=10, ylabel="Branch utilization (%)", xrotation=20)

        push!(plots, p)
    end
    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 800), margin=5mm)
    return final_plot
end

normalized_pf_plot = plot_normalized_pf_whisker(all_data, cases)
savefig(normalized_pf_plot, "$PATH/normalized_pf_custom_plot.pdf")
display(normalized_pf_plot)


















normalized_pf_plot = plot_normalized_pf_whisker(all_data, cases)
savefig(normalized_pf_plot, "$PATH/normalized_pf_custom_plot.pdf")
display(normalized_pf_plot)


function plot_normalized_pf_bar(all_data, cases)
    plots = []
    RE_levels = sort(unique([v["RE"] for v in values(cases)]))

    for RE in RE_levels
        all_pf_values = []
        labels = []
        numbers = []

        filtered_cases = [k for (k, v) in cases if v["RE"] == RE]
        for case_label in filtered_cases
            branch_data = all_data[case_label]["results"]["branch_t"]["pf"]
            rate_a_data = all_data[case_label]["results"]["branch"]["rate_a"][:, 2]
            constraint = all_data[case_label]["constraint"]
            number = parse(Int, all_data[case_label]["number"])

            if !isempty(branch_data) && !isempty(rate_a_data)
                pf_values = vec(Matrix(branch_data[:, 2:end]))
                rate_a_repeated = repeat(rate_a_data, outer=size(branch_data, 1))
                normalized_pf = pf_values ./ rate_a_repeated

                push!(all_pf_values, normalized_pf)
                push!(labels, constraint)
                push!(numbers, number)
            else
                println("Daten fehlen für $case_label")
            end
        end
        
        sorted_indices = sortperm(numbers)
        sorted_labels = [labels[i] for i in sorted_indices]
        sorted_pf_values = [all_pf_values[i] for i in sorted_indices]

        means = mean.(sorted_pf_values)
        std_devs = std.(sorted_pf_values)
        p = bar(sorted_labels, means, yerr=std_devs, title="RE = $RE", color=:blue, 
                legend=false, ylabel="Branch utilization (%)")

        plot!(p, titlefontsize=14, guidefontsize=12, tickfontsize=10, legendfontsize=10, xrotation=20)
        push!(plots, p)
    end

    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 800), margin=5mm)
    return final_plot
end

normalized_pf_plot = plot_normalized_pf_bar(all_data, cases)
savefig(normalized_pf_plot, "$PATH/normalized_pf_custom_plot.pdf")
display(normalized_pf_plot)


function plot_normalized_pf_violin(all_data, cases)
    plots = []
    RE_levels = sort(unique([v["RE"] for v in values(cases)]))

    for RE in RE_levels
        all_pf_values = []
        labels = []
        numbers = []

        filtered_cases = [k for (k, v) in cases if v["RE"] == RE]
        for case_label in filtered_cases
            branch_data = all_data[case_label]["results"]["branch_t"]["pf"]
            rate_a_data = all_data[case_label]["results"]["branch"]["rate_a"][:, 2]
            constraint = all_data[case_label]["constraint"]
            number = parse(Int, all_data[case_label]["number"])

            if !isempty(branch_data) && !isempty(rate_a_data)
                pf_values = vec(Matrix(branch_data[:, 2:end]))
                rate_a_repeated = repeat(rate_a_data, outer=size(branch_data, 1))
                normalized_pf = pf_values ./ rate_a_repeated

                push!(all_pf_values, normalized_pf)
                push!(labels, constraint)
                push!(numbers, number)
            else
                println("Daten fehlen für $case_label")
            end
        end

        sorted_indices = sortperm(numbers)
        sorted_labels = [labels[i] for i in sorted_indices]
        sorted_pf_values = [all_pf_values[i] for i in sorted_indices]

        p = StatsPlots.violin(sorted_labels, sorted_pf_values, title="RE = $RE", legend=false,
                   color=:blue, ylabel="Branch utilization (%)")

        plot!(p, titlefontsize=14, guidefontsize=12, tickfontsize=10, legendfontsize=10, xrotation=20)
        push!(plots, p)
    end

    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 800), margin=5mm)
    return final_plot
end

normalized_pf_plot = plot_normalized_pf_violin(all_data, cases)
savefig(normalized_pf_plot, "$PATH/normalized_pf_violin_plot.pdf")
display(normalized_pf_plot)


############################################################# E_I_area scatter   #####################################################################################

function create_inertia_scatter_plots(cases, selected_cases, selected_RE, c=25)
    plots = []

    for RE in selected_RE
        for case_id in selected_cases
            filtered_cases = [k for (k, v) in cases if v["RE"] == RE && parse(Int, split(k, "_")[2]) == case_id]
            for case_key in filtered_cases
                p = scatter(title="RE = $RE", xlabel="Power difference in the Regions (p.u.)", ylabel="Region Inertia (p.u. s)",
                            titlefontsize=14, xlabelfontsize=12, ylabelfontsize=12,
                            guidefontsize=10, tickfontsize=10, legendfontsize=12,
                            legend=:bottomright, size=(1200, 800), margin=5mm)

                delta_p = DataFrame(cases[case_key]["results"]["bus_t"]["delta_p_area"])
                e_i_area = DataFrame(cases[case_key]["results"]["bus_t"]["E_I_area"])

                delta_p_values = vec(delta_p[:, 2])
                e_i_area_values = vec(e_i_area[:, 2])

                dp = maximum(abs.(delta_p_values))
                xlims!(p, -dp, dp)

                scatter!(p, delta_p_values, e_i_area_values, label=cases[case_key]["constraint"],
                         markersize=6)

                plot!(p, [0, dp], [-dp*c, dp*c], label=false, color=:black, linewidth=3)
                plot!(p, [-dp, 0], [dp*c, -dp*c], label=false, color=:black, linewidth=3)

                push!(plots, p)
            end
        end
    end

    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 400 * length(plots)))
    return final_plot
end

selected_cases = [6, 8]
selected_RE = ["45%"]
Inertia_scatter_plot_neu = create_inertia_scatter_plots(cases, selected_cases, selected_RE)
savefig(Inertia_scatter_plot_neu, "$PATH/inertia_scatter_plot_neu.pdf")
display(Inertia_scatter_plot_neu)



function create_scatter_plots(cases)
    plots = []
    RE_levels = sort(unique([v["RE"] for v in values(cases)]))
    colors = ["red", "blue"]
    markers = [:circle, :square]

    for RE in RE_levels
        p = scatter(title="RE $RE", xlabel="Power deference in the Regions (p.u.)", ylabel="Region Inertia (p.u. s)",
                    titlefontsize=14, xlabelfontsize=12, ylabelfontsize=12,
                    guidefontsize=10, tickfontsize=10, legendfontsize=12,
                    legend=:topright, size=(1200, 800), margin=5mm)

        filtered_cases = [k for (k, v) in cases if v["RE"] == RE && (v["number"] in ["6"])]
        sorted_filtered_cases = sort(filtered_cases, by = x -> parse(Int, split(x, "_")[2]))

        for (case_index, case_key) in enumerate(sorted_filtered_cases)
            constraint = cases[case_key]["constraint"]
            delta_p = DataFrame(cases[case_key]["results"]["bus_t"]["delta_p_area"])
            e_i_area = DataFrame(cases[case_key]["results"]["bus_t"]["E_I_area"])

            delta_p_values = vec(delta_p[:, 2])
            e_i_area_values = vec(e_i_area[:, 2])

            scatter!(p, delta_p_values, e_i_area_values, label="$constraint",
                     markersize=8, color=colors[mod1(case_index, length(colors))],
                     marker=markers[mod1(case_index, length(markers))])

            fit_result = Polynomials.fit(Polynomial, delta_p_values, e_i_area_values, 1)
            x_fit = range(minimum(delta_p_values), maximum(delta_p_values), length=100)
            y_fit = fit_result.(x_fit)
            plot!(p, x_fit, y_fit, label=false, linestyle=:dash, color=colors[mod1(case_index, length(colors))])
        end

        push!(plots, p)
    end

    final_plot = plot(plots..., layout=(length(plots), 1), size=(1200, 400 * length(RE_levels)))
    return final_plot
end

Inertia_scatter_plot = create_scatter_plots(cases)
savefig(Inertia_scatter_plot, "$PATH/ineria_scatter_plot.pdf")
display(Inertia_scatter_plot)

##################################################################################################################################################





function create_adjacency_matrix(df, from_col, to_col, value_col)
    nodes = sort(unique(vcat(df[!, from_col], df[!, to_col])))
    mat = zeros(Float64, length(nodes), length(nodes))

    for row in eachrow(df)
        from_idx = findfirst(==(row[from_col]), nodes)
        to_idx = findfirst(==(row[to_col]), nodes)
        value = coalesce(row[value_col], 0.0)
        mat[from_idx, to_idx] += value
    end

    return mat, nodes
end

function filter_and_create_matrix(results_df, constraint, RE)
    # Daten filtern
    filtered_df = results_df[(results_df[!, :constraint] .== constraint) .& (results_df[!, :RE] .== RE), :]
    adj_matrix, nodes = create_adjacency_matrix(filtered_df, :f_bus, :t_bus, :pf_sum)

    return adj_matrix, nodes
end

constraint_to_use = "Region Islanding 100%" # "System Inertia", "Region Islanding 50%", "Region Islanding 100%", "Unit Commitment"
RE_to_use = "50%"
adj_matrix, nodes = filter_and_create_matrix(results_branch_df, constraint_to_use, RE_to_use)

heatmap(adj_matrix, color = :viridis, aspect_ratio = 1, # color = :blues, :reds, :greens, :grays, :purples, :oranges, :viridis
        xticks = (1:length(nodes), nodes), yticks = (1:length(nodes), nodes),
        xlabel = "Bus", ylabel = "Bus",
        size = (1000, 800),
        xlims = (0.5, length(nodes) + 0.5),
        ylims = (0.5, length(nodes) + 0.5),
        xguidefontsize = 12, yguidefontsize = 12, xtickfontsize = 10, ytickfontsize = 10,
        
)       
savefig("$PATH/$(constraint_to_use)_$(RE_to_use)_plot_heatmap_pf_sum.pdf")



function calculate_and_plot_branch_load(test_case::String, all_data::Dict)

    pf_data_original = all_data[test_case]["results"]["branch_t"]["pf"]
    rate_a_data = all_data[test_case]["results"]["branch"]["rate_a"]
    rate_a_dict = Dict(rate_a_data.id .=> rate_a_data.value)
    pf_data = copy(pf_data_original)
    
    column_names = sort(names(pf_data)[2:end], by = x -> parse(Int, x))
    column_symbols = [:nw; Symbol.(column_names)]
    pf_data = pf_data[:, column_symbols]
    sort!(pf_data, :nw, by = x -> parse(Int, string(x)))

    for col in column_names
        if rate_a_dict[col] > 0
            pf_data[:, col] .= abs.(pf_data[:, col] ./ rate_a_dict[col]) * 100
        else
            println("Warnung: Ungültiger rate_a Wert für Spalte $(col): $(rate_a_dict[col])")
        end
    end
    
    heatmap(column_names,pf_data.nw, Matrix(pf_data[:, 2:end]), 
            c=:cividis, 
            xlabel="Bumber of Bus", ylabel="Day of the year",
            title=" ",
            size=(1200, 800),
            left_margin=5mm,
            xguidefontsize=12, yguidefontsize=12,
            xtickfontsize=10, ytickfontsize=10)
end


test_case = "case_1_50" # "case_1_35", "case_1_50", "case_2_35", "case_2_50", "case_6_35", "case_6_50", "case_8_35", "case_8_50"
plot_branch_laod = calculate_and_plot_branch_load(test_case, all_data)
display(plot_branch_laod)
savefig(plot_branch_laod, "$PATH/plot_branch_load_$test_case.pdf")


#===========================================================================================================================




#=
global gen_df = DataFrame()

for c in keys(cases)
    println(c)
    gen_c = all_data[c]["results"]["gen"]
    gen_c_df = DataFrame()
    kk = 0
    for k in keys(gen_c)
        println(k)
        gen_c[k] = rename(gen_c[k], :value => k)
        if kk == 0
            gen_c_df[!, :id] = gen_c[k][!, :id]
            gen_c_df[!, "constraint"] .= cases[c]["constraint"]
            gen_c_df[!, "RE"] .= cases[c]["RE"]
            gen_c_df[!, k] = gen_c[k][!, k]
            kk = 1
        else
            gen_c_df = innerjoin(gen_c_df, gen_c[k], on= :id)
        end
    end 
    global gen_df = vcat(gen_df, gen_c_df)
end

carriers_file = ".\\test\\data\\matpower\\multi_nw\\bus_37\\pypsa-eur_multinetwork_37\\carrier_indices.csv"
carrier_names = CSV.read(carriers_file, DataFrame, header =0)
gen_df[!, :carrier_names] = [carrier_names[v + 1, "Column1"] for v in gen_df[:, :carrier]]


gen_df[!, :pExp] = (gen_df[!, :nE] - gen_df[!, :n0]) .* gen_df[!, :P_b_nom] / 1e3
g_gen_df = groupby(gen_df, [:constraint, :RE, :carrier_names])
g_exp_df = filter(row -> row.pExp > 0, gen_df)
exp_df = combine(g_gen_df, :pExp => sum)
exp_df = filter(row -> row.pExp_sum > 0, exp_df)

g_exp_df35 = filter(row->row.RE=="35%", g_exp_df)
exp_df35 = filter(row->row.RE=="35%",exp_df)



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

pg_sums = calculate_weighted_pg_sums(all_data, cases)

merged_df = leftjoin(gen_df, pg_sums, on = [:id, :constraint, :RE])

gen_df = merged_df



=#


function plot_cost(all_data, plot_type="bar", save_path=".")
    p1 = plot(title = "35% RE", xlabel = " ", ylabel = "system cost (%)", legend = false, xrotation = 45)
    p2 = plot(title = "70% RE", xlabel = " ", ylabel = " ", legend = false, xrotation = 45)

    xticks1 = Int64[]
    xticks2 = Int64[]
    xlabels1 = String[]
    xlabels2 = String[]
    i1 = 1
    i2 = 1

    sorted_data = sort(collect(all_data), by = x -> parse(Int, split(x[1], '_')[2]))

    reference_costs = Dict()

    case_legend = Dict(
        "case_1" => "unit commitment",
        "case_2" => "system inertia",
        "case_3" => "small disturbance with weighted load area",
        "case_4" => "small disturbance with weighted equal area",
        "case_5" => "large disturbance with β = 0.25",
        "case_6" => "large disturbance with β = 0.50",
        "case_7" => "large disturbance with β = 0.75",
        "case_8" => "large disturbance with β = 1.00"
    )

    legends_added = Set()
    labels = []

    for (case_label, case_data) in sorted_data
        if haskey(case_data, "results") && hasproperty(case_data["results"]["general_var"], :objective)
            cost = case_data["results"]["general_var"].objective[1]
            subcase_label = split(case_label, '_')[2]
            println("subcase_label: ", subcase_label)
            if startswith(case_label, "case_1")
                reference_costs[subcase_label] = cost
            end
            cost_percentage = (cost / reference_costs[subcase_label]) * 100
            if endswith(case_label, "35")
                p = p1
                i = i1
                xticks = xticks1
                xlabels = xlabels1
            else
                p = p2
                i = i2
                xticks = xticks2
                xlabels = xlabels2
            end
            label_key = join(split(case_label, '_')[1:end-1], '_')
            label = in(label_key, legends_added) ? "" : case_legend[label_key]
            
            if plot_type == "scatter"
                scatter!(p, [i], [cost_percentage], marker = :circle, label = label)
            elseif plot_type == "bar"
                bar!(p, [i], [cost_percentage], label = label)
            end
            if !in(case_label, legends_added)
                push!(labels, label)
            end
            push!(legends_added, case_label)
            # push!(xticks, i)
            push!(xlabels, replace(subcase_label, r"35|70" => ""))
            if endswith(case_label, "35")
                i1 += 1
            else
                i2 += 1
            end
        else
            println("Missing 'objective' in 'results' for case '$case_label'.")
        end
    end
    xticks!(p1, Float64[], String[])
    xticks!(p2, Float64[], String[])

    combined_plot = plot(p1, p2, layout = (1, 2), legend = false)
    plot!(combined_plot, legend = :bottomleft, labels = labels) # bottomright, topleft, topright, bottomleft
    display(combined_plot)
    savefig(joinpath(save_path, "cost_plot.pdf"))
end


plot_cost(all_data, "bar", "C:/Users/Taiseer/Desktop/Master NEE/MASTERARBEIT/LaTex_schriftliche_Ausarbeitung_MA/plot_results")  # or "scatter" for scatter plot




function calculate_and_plot_pg(all_data)
    carriers = Dict(
        0 => "CCGT",
        1 => "OCGT",
        2 => "offwind",
        3 => "oil",
        4 => "onwind",
        5 => "ror",
        6 => "solar",
        7 => "lignite",
        8 => "biomass",
        9 => "coal",
        10 => "nuclear",
        11 => "geothermal",
        12 => "hydro",
        13 => "PHS_turb",
        14 => "PHS_pump"
    )

    tech_production = DataFrame(Tech = String[], Production = Float64[])

    for (key, case) in all_data
        pg = case["results"]["gen_t"]["pg"]
        weights = case["results"]["weight"]
        total_production_all_tech = sum([sum(pg[:, col]) for col in 2:ncol(pg)]) * weights[1, :weight][1]
        carrier = case["results"]["gen"]["carrier"]
        carrier[!, :tech] = map(x -> get(carriers, x, "Unknown"), carrier[!, :value])

        for nw in unique(pg[:, :nw])
            current_weight = weights[weights[:, :nw] .== nw, :weight][1]
            for col in 2:ncol(pg)
                tech = carrier[!, :tech][col - 1]
                column_name = names(pg)[col]
                production = sum(pg[pg[:, :nw] .== nw, column_name] * current_weight)
                production_percent = (production / total_production_all_tech) * 100
                push!(tech_production, (tech, production_percent))
            end
        end
    end

    # Aggregation der Produktionswerte nach Technologie
    tech_grouped = groupby(tech_production, :Tech)
    final_production = combine(tech_grouped, :Production => sum => :TotalProduction)

    # Daten für das Plotten vorbereiten
    tech_labels = unique(tech_production.Tech)
    production_values = [final_production[final_production.Tech .== tech, :TotalProduction][1] for tech in tech_labels]
    production_array = hcat(production_values...)
    p = groupedbar(production_array, bar_position=:stack, bar_width=0.7, labels=false)
    ylabel!("Energy in %")
    display(p)
    savefig(p, "energy_production_stack_plot.pdf")
end
calculate_and_plot_pg(all_data)




function calculate_and_plot_pg(all_data)
    carriers = Dict(
        0 => "CCGT",
        1 => "OCGT",
        2 => "offwind",
        3 => "oil",
        4 => "onwind",
        5 => "ror",
        6 => "solar",
        7 => "lignite",
        8 => "biomass",
        9 => "coal",
        10 => "nuclear",
        11 => "geothermal",
        12 => "hydro",
        13 => "PHS_turb",
        14 => "PHS_pump"
    )

    all_techs = values(carriers) |> unique |> sort
    all_case_productions = [Float64[] for _ in all_techs]

    for (key, case) in all_data
        tech_production = DataFrame(Tech = String[], Production = Float64[])

        pg = case["results"]["gen_t"]["pg"]
        weights = case["results"]["weight"]
        carrier = case["results"]["gen"]["carrier"]
        carrier[!, :tech] = map(x -> get(carriers, x, "Unknown"), carrier[!, :value])

        for nw in unique(pg[:, :nw])
            current_weight = weights[weights[:, :nw] .== nw, :weight][1]
            for col in 2:ncol(pg)
                tech = carrier[!, :tech][col - 1]
                column_name = names(pg)[col]
                production = sum(pg[pg[:, :nw] .== nw, column_name] * current_weight)
                push!(tech_production, (tech, production))
            end
        end

        total_production = sum(tech_production[:, :Production])
        for (i, tech) in enumerate(all_techs)
            tech_total = sum(tech_production[tech_production.Tech .== tech, :Production])
            push!(all_case_productions[i], tech_total / total_production * 100)
        end
    end

    p = groupedbar(all_case_productions, bar_position=:stack, bar_width=0.7, labels=all_techs)
    ylabel!(p, "Energy in %")
    xlabel!(p, "Cases")
    display(p)
    savefig(p, "energy_production_stack_plot.pdf")
end

calculate_and_plot_pg(all_data)

=################################################################################################################################












#=


function load_all_results(path)
    closeall()

    # path = "results\\results_bus_$(bus_system)\\relax_$(relax_integrality)"
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
    
    return all_data , mn_data_main
end

all_data, mn_data_main = load_all_results(path);

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

# all_sol_data = add_missing_values_to_sol_data(all_data, mn_data_main)[1];
# mn_data_main = add_missing_values_to_sol_data(all_data, mn_data_main)[2];
all_sol_data, mn_data_main = add_missing_values_to_sol_data(all_data, mn_data_main);


function create_dataframes(all_sol_data::Dict{Any, Any})
    dfs = Dict()
    
    time_dependent_variables = Dict(
        "gen" => Set(["pg", "gen_status", "gen_startup", "gen_shutdown"]),
        "branch" => Set(["pf", "pt"]),
        "bus" => Set(["va"]),
        "dcline" => Set(["pfl", "ptl"])
    )
    
    for (case_label, results_v) in all_sol_data
        for (re_inj_int, result) in results_v
            case_re_label = "$(case_label)_$(re_inj_int)"
            dfs[case_re_label] = Dict()

            if result["results"]["termination_status"] == JuMP.OPTIMAL
                for comp in ["gen", "bus", "branch", "dcline"]
                    dfs[case_re_label]["$(comp)_t"] = Dict()
                    dfs[case_re_label][comp] = Dict()
                    
                    all_vars = Set()
                    for (_, comp_data) in result["results"]["solution"]["nw"]
                        if haskey(comp_data, comp)
                            for (_, details) in comp_data[comp]
                                union!(all_vars, keys(details))
                            end
                        end
                    end

                    non_time_dependent_vars = setdiff(all_vars, time_dependent_variables[comp])

                    for var in time_dependent_variables[comp]
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
                        sort!(df, :nw, by = x -> parse(Int, string(x)))
                        dfs[case_re_label]["$(comp)_t"][var] = df
                    end

                    for var in non_time_dependent_vars
                        unique_data = Dict()
                        for (_, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                for (id, details) in comp_data[comp]
                                    unique_data[id] = get(details, var, NaN)
                                end
                            end
                        end

                        if isempty(unique_data)
                            continue
                        end
                        
                        id_keys = sort(collect(keys(unique_data)), by = x -> parse(Int, x))
                        values_list = [unique_data[id] for id in id_keys]
                        df = DataFrame(id = id_keys, value = values_list)
                        dfs[case_re_label][comp][var] = df
                    end

                end
            end
        end
    end

    return dfs
end

dfs = create_dataframes(all_sol_data)

# save dsf to jld2 file
JLD2.save("results\\results_bus_$(bus_system)\\relax_$(relax_integrality)\\results_as_df/dfs.jld2", "dfs", dfs)

function plot_cost(all_data, plot_type="bar", save_path=".")
    p1 = plot(title = "35% RE", xlabel = " ", ylabel = "system cost (%)", legend = false, xrotation = 45)
    p2 = plot(title = "70% RE", xlabel = " ", ylabel = " ", legend = false, xrotation = 45)

    xticks1 = Int64[]
    xticks2 = Int64[]
    xlabels1 = String[]
    xlabels2 = String[]
    i1 = 1
    i2 = 1

    sorted_data = sort(collect(all_data), by = x -> parse(Int, split(x[1], '_')[2]))

    reference_costs = Dict()

    case_legend = Dict(
        "case_1" => "unit commitment",
        "case_2" => "system inertia",
        "case_3" => "small disturbance with weighted load area",
        "case_4" => "small disturbance with weighted equal area",
        "case_5" => "large disturbance with β = 0.25",
        "case_6" => "large disturbance with β = 0.50",
        "case_7" => "large disturbance with β = 0.75",
        "case_8" => "large disturbance with β = 1.00"
    )

    legends_added = Set()
    labels = []

    for (case_label, case_data) in sorted_data
        sorted_subcases = sort(collect(case_data), by = x -> parse(Int, x[1]))
        for (subcase_label, subcase_data) in sorted_subcases
            if haskey(subcase_data, "results") && haskey(subcase_data["results"], "objective")
                cost = subcase_data["results"]["objective"]
                if case_label == "case_1"
                    reference_costs[subcase_label] = cost
                end
                cost_percentage = (cost / reference_costs[subcase_label]) * 100
                if endswith(subcase_label, "35")
                    p = p1
                    i = i1
                    xticks = xticks1
                    xlabels = xlabels1
                else
                    p = p2
                    i = i2
                    xticks = xticks2
                    xlabels = xlabels2
                end
                label = in(case_label, legends_added) ? "" : case_legend[case_label]
                if plot_type == "scatter"
                    scatter!(p, [i], [cost_percentage], marker = :circle, label = label)
                elseif plot_type == "bar"
                    bar!(p, [i], [cost_percentage], label = label)
                end
                if !in(case_label, legends_added)
                    push!(labels, label)
                end
                push!(legends_added, case_label)
                push!(xticks, i)
                push!(xlabels, replace(subcase_label, r"35|70" => ""))
                if endswith(subcase_label, "35")
                    i1 += 1
                else
                    i2 += 1
                end
            else
                println("Missing 'objective' in 'results' for case '$(case_label)_$subcase_label'.")
            end
        end
    end
    xticks!(p1, xticks1, xlabels1)
    xticks!(p2, xticks2, xlabels2)

    combined_plot = plot(p1, p2, layout = (1, 2), legend = false)
    plot!(legend = :bottomleft) # bottomright, topleft, topright, bottomleft
    display(combined_plot)
    savefig(joinpath(save_path, "cost_plot.pdf"))
end

plot_cost(all_sol_data, "bar", "C:/Users/Taiseer/Desktop/Master NEE/MASTERARBEIT/LaTex_schriftliche_Ausarbeitung_MA/plot_results")  # or "scatter" for scatter plot


function plot_cost(all_data, plot_type="bar", save_path=".")
    p = plot(title = "Kostenvergleich", xlabel = " ", ylabel = "Kosten", legend = :outertopright, xrotation = 45)
    xticks = Int64[]
    xlabels = String[]
    i = 1

    sorted_data = sort(collect(all_data), by = x -> parse(Int, split(x[1], '_')[2]))

    reference_costs = Dict()

    for (case_label, case_data) in sorted_data
        sorted_subcases = sort(collect(case_data), by = x -> parse(Int, x[1]))
        for (subcase_label, subcase_data) in sorted_subcases
            if haskey(subcase_data, "results") && haskey(subcase_data["results"], "objective")
                cost = subcase_data["results"]["objective"]
                if case_label == "case_1"
                    reference_costs[subcase_label] = cost
                end
                if plot_type == "scatter"
                    scatter!(p, [i], [cost], label = "$(case_label)_$subcase_label", marker = :circle)
                elseif plot_type == "bar"
                    bar!(p, [i], [cost], label = "$(case_label)_$subcase_label")
                    if haskey(reference_costs, subcase_label)
                        cost_increase_percentage = ((cost - reference_costs[subcase_label]) / reference_costs[subcase_label]) * 100
                        annotate!(p, [(i, cost, text(string(round(cost_increase_percentage, digits=2), "%"), 8))])  # Adjust the position here
                    
                    end
                end
                push!(xticks, i)
                push!(xlabels, "$(case_label)_$subcase_label")
                i += 1
            else
                println("Missing 'objective' in 'results' for case '$(case_label)_$subcase_label'.")
            end
        end
    end
    xticks!(p, xticks, xlabels)
    display(p)
    savefig(p, joinpath(save_path, "cost_plot.pdf"))  # Save the plot as a PDF
end
plot_cost(all_sol_data, "bar", "C:/Users/Taiseer/Desktop/Master NEE/MASTERARBEIT/LaTex_schriftliche_Ausarbeitung_MA/plot_results")  # or "scatter" for scatter plot



function plot_solver_time(all_data, plot_type="bar", save_path=".")
    p = plot(title = "Vergleich der Lösungszeit", xlabel = " ", ylabel = "Zeit (s)", legend = :outertopright, xrotation = 45)
    xticks = Int64[]
    xlabels = String[]
    i = 1

    sorted_data = sort(collect(all_data), by = x -> parse(Int, split(x[1], '_')[2]))

    reference_times = Dict()

    for (case_label, case_data) in sorted_data
        sorted_subcases = sort(collect(case_data), by = x -> parse(Int, x[1]))
        for (subcase_label, subcase_data) in sorted_subcases
            if haskey(subcase_data, "results") && haskey(subcase_data["results"], "solve_time")
                solver_time = subcase_data["results"]["solve_time"]
                if case_label == "case_1"
                    reference_times[subcase_label] = solver_time
                end
                if plot_type == "scatter"
                    scatter!(p, [i], [solver_time], label = "$(case_label)_$subcase_label", marker = :circle)
                elseif plot_type == "bar"
                    bar!(p, [i], [solver_time], label = "$(case_label)_$subcase_label")
                    if haskey(reference_times, subcase_label)
                        time_increase_percentage = ((solver_time - reference_times[subcase_label]) / reference_times[subcase_label]) * 100
                        annotate!(p, [(i, solver_time, text(string(round(time_increase_percentage, digits=2), "%"), 8))])
                    end
                end
                push!(xticks, i)
                push!(xlabels, "$(case_label)_$subcase_label")
                i += 1
            else
                println("Missing 'solve_time' in 'results' for case '$(case_label)_$subcase_label'.")
            end
        end
    end
    xticks!(p, xticks, xlabels)
    display(p)
    savefig(p, joinpath(save_path, "solver_time_plot.pdf"))  # Save the plot as a PDF
end
plot_solver_time(all_sol_data, "bar", "C:\\Users\\Taiseer\\Desktop\\Master NEE\\MASTERARBEIT\\LaTex_schriftliche_Ausarbeitung_MA\\plot_results")



function plot_pg_data(dfs, pg_plot_options, save_path=".")
    cases = get(pg_plot_options, "cases", [])
    component = get(pg_plot_options, "component", "")
    variable = get(pg_plot_options, "variable", "")
    generator_ids = get(pg_plot_options, "generator_ids", [])
    nw_range = get(pg_plot_options, "nw_range", nothing)
    plot_type = get(pg_plot_options, "plot_type", :line)
    plot_name = get(pg_plot_options, "plot_name", "")
    p = plot() 

    if plot_name == "pg_pro_gen"
        xlabel = "Laststunden"
        ylabel = "Eruegerleistung in p.u."
        title = "Vergleich der Erzeugerleistung"
        p = plot(title = title, xlabel = xlabel, ylabel = ylabel, legend = :outertopright)

        for case_label in cases
            df = dfs[case_label][component][variable]

            nw_col = parse.(Int, df[!, :nw])
            
            if nw_range !== nothing
                indices = findall(x -> x in nw_range, nw_col)
                df = df[indices, :]
                nw_col = nw_col[indices]
            end

            if !all(id in names(df) for id in generator_ids)
                println("One or more generator_ids not found in DataFrame for case $case_label.")
                continue
            end

            filtered_df = df[:, generator_ids]
            for id in generator_ids
                variable_data = filtered_df[!, id]
                if isempty(variable_data)
                    println("No data available for generator ID $id in case $case_label.")
                    continue
                end

                if plot_type == :line
                    plot!(p, nw_col, variable_data, label = "$case_label, $id")
                elseif plot_type == :scatter
                    scatter!(p, nw_col, variable_data, label = "$case_label, $id")
                else
                    println("Invalid plot type '$plot_type'.")
                end
            end
        end
    elseif plot_name == "pg_pro_carrier"
        xlabel = ""
        ylabel = "Erzeugerleistung in p.u."
        title = "Vergleich der Erzeugerleistung"
        p = plot(title = title, xlabel = xlabel, ylabel = ylabel, legend = false)
    
        carrier_labels = Dict("0" => "CGT", "1" => "CGT", "2" => "BM", "3" => "WE", "4" => "WE", "5" => "PV", "6" => "HD", "7" => "PHS", "8" => "PHS", "9" => "COAL")
        
        all_carriers = sort(unique(values(carrier_labels)))
        carrier_colors = Dict(
            "CGT" => RGB(0.000, 0.502, 0.502),  # Dunkles Türkis
            "BM"  => RGB(0.690, 0.188, 0.376),  # Dunkles Rosa
            "WE"  => RGB(0.416, 0.353, 0.804),  # Dunkles Indigo
            "PV"  => RGB(0.855, 0.647, 0.125),  # Goldgelb
            "HD"  => RGB(0.988, 0.894, 0.423),  # Sonnengelb
            "PHS" => RGB(0.627, 0.322, 0.176),   # Siena
            "COAL" => RGB(0.000, 0.000, 0.000)  # Schwarz
            "NUCLEAR" =
        )
        
        case_offsets = Dict()
    
        for (case_index, case_label) in enumerate(cases)
            case_offset = case_index * length(all_carriers) * 1.3
            case_offsets[case_label] = case_offset
    
            df_gen = dfs[case_label][component][variable]
            df_carrier = dfs[case_label]["gen"]["carrier"]
            df_carrier_map = Dict(row.id => carrier_labels[string(row.value)] for row in eachrow(df_carrier))
    
            nw_col = parse.(Int, df_gen[!, :nw])
            indices = findall(x -> x in nw_range, nw_col)
    
            carrier_pg = Dict{String, Float64}()
            for gen_id_str in names(df_gen)
                if gen_id_str == "nw"
                    continue
                end
                gen_id = parse(Int, gen_id_str)
                if haskey(df_carrier_map, string(gen_id))
                    carrier_id = df_carrier_map[string(gen_id)]
                    summed_pg = sum(df_gen[indices, gen_id_str])
                    carrier_pg[carrier_id] = get(carrier_pg, carrier_id, 0.0) + summed_pg
                end
            end
    
            for (carrier_name, pg_value) in carrier_pg
                carrier_index = findfirst(x -> x == carrier_name, all_carriers)
                pos = case_offsets[case_label] + carrier_index
                bar!(p, [pos], [pg_value], width=0.8, color=carrier_colors[carrier_name], label="")
                annotate!(p, [(pos, pg_value + 0.5, text(carrier_name, :bottom, 10))])            
            end
        end
    
        xticks = [case_offset + (length(all_carriers) * 0.5) for case_offset in values(case_offsets)]
        xticks!(p, xticks, cases)
    
    end    
    display(p)
    save_path = joinpath(save_path, "$plot_name.pdf")
    savefig(p, save_path) 
end


pg_plot_options = Dict(
    "cases" => ["case_1_35", "case_1_70"],
    "component" => "gen_t",
    "variable" => "pg",
    "generator_ids" => [],
    "plot_name" => "pg_pro_carrier",  # "pg_pro_gen" or "pg_pro_carrier"
    "plot_type" => :bar, # line or scatter or bar
    "nw_range" => 0:47
)

plot_pg_data(dfs, pg_plot_options, "C:/Users/Taiseer/Desktop/Master NEE/MASTERARBEIT/LaTex_schriftliche_Ausarbeitung_MA/plot_results")



function plot_power_flows(pf_options, dfs)
    case_label = pf_options["case"]
    nw = pf_options["nw"]
    num_edges = size(dfs[case_label]["branch"]["f_bus"], 1)
    num_buses = size(dfs[case_label]["bus"]["bus_i"], 1)

    bus_data = dfs[case_label]["bus"]["bus_i"]
    branch_data = dfs[case_label]["branch"]
    pf = dfs[case_label]["branch_t"]["pf"]
    pt = dfs[case_label]["branch_t"]["pt"]
    rate_a = branch_data["rate_a"]

    g = Graphs.SimpleGraph(num_buses)

    for i in 1:num_edges
        f_bus = branch_data["f_bus"][i, :value]
        t_bus = branch_data["t_bus"][i, :value]
        Graphs.add_edge!(g, f_bus, t_bus)
    end

    edge_colors = Vector{RGB}(undef, num_edges)
    edge_labels = Vector{String}(undef, num_edges)

    for i in 1:num_edges
        pf_value = pf[nw, i + 1]
        rate_value = rate_a[i, :value]

        utilization = abs(pf_value) / rate_value
        color = utilization > 0.8 ? RGB(1, 0, 0) : utilization > 0.5 ? RGB(1, 1, 0) : RGB(0, 1, 0)

        edge_colors[i] = color
        edge_labels[i] = @sprintf("%.2f", pf_value)
    end

    offset = 0.1
    Random.seed!(35)
    layout = (x) -> spring_layout(g)
    graph_plot = gplot(g,
                       layout=layout,
                       edgelabel=edge_labels,
                       edgestrokec=edge_colors,
                       nodelabel=bus_data[:, :value],
                       nodefillc=colorant"purple",
                       edgelabeldistx=offset,
                       edgelabeldisty=offset)

    display(graph_plot)
    save_path = "C:/Users/Taiseer/Desktop/Master NEE/MASTERARBEIT/LaTex_schriftliche_Ausarbeitung_MA/plot_results"
    draw(PDF(joinpath(save_path, "pf_bus_$(num_buses)_$(case_label)_$(nw).pdf")), graph_plot)

end

pf_options = Dict(
    "case" => "case_1_35",
    "nw" => 20
)

plot_power_flows(pf_options, dfs)




#=
function create_dataframes(all_sol_data::Dict{Any, Any})
    dfs = Dict()
    
    # Define time-dependent variables for each component
    time_dependent_variables = Dict(
        "gen" => ["pg", "gen_status", "gen_startup", "gen_shutdown"],
        "branch" => ["pf", "pt"],
        "bus" => ["va"]
    )
    
    for (case_label, results_v) in all_sol_data
        for (re_inj_int, result) in results_v
            case_re_label = "$(case_label)_$(re_inj_int)"
            dfs[case_re_label] = Dict()

            if result["results"]["termination_status"] == JuMP.OPTIMAL

                for comp in ["gen", "bus", "branch"]
                    dfs[case_re_label]["$(comp)_t"] = Dict()
                    dfs[case_re_label][comp] = Dict()
                    
                    all_vars = Set()
                    for (_, comp_data) in result["results"]["solution"]["nw"]
                        if haskey(comp_data, comp)
                            for (_, details) in comp_data[comp]
                                all_vars = union(all_vars, Set(keys(details)))
                            end
                        end
                    end

                    non_time_dependent_vars = setdiff(all_vars, time_dependent_variables[comp])

                    for var in time_dependent_variables[comp]
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
                        
                        # Erstellen einer Liste von eindeutigen nw-Schlüsseln
                        nw_keys = sort(parse.(Int, collect(keys(data))))

                        # Erstellen einer Liste von eindeutigen id-Schlüsseln
                        id_keys = Set{Int}()
                        for (nw, comp_data) in result["results"]["solution"]["nw"]
                            if haskey(comp_data, comp)
                                id_keys = union(id_keys, parse.(Int, collect(keys(comp_data[comp]))))
                            end
                        end
                        id_keys = sort(collect(id_keys))

                        # Erstellen der Matrix aus den Daten
                        matrix = [get(get(data, string(nw), Dict()), string(id), NaN) for nw in nw_keys, id in id_keys]

                        # Erstellen des DataFrames
                        df = DataFrame(matrix, Symbol.(id_keys))
                        insertcols!(df, 1, :nw => nw_keys)

                        dfs[case_re_label]["$(comp)_t"][var] = df
                    end

                    # Non-time-dependent variables processing
                    for var in non_time_dependent_vars
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
                        sort!(df, :nw, by = x -> parse(Int, string(x)))
                        dfs[case_re_label][comp][var] = df
                    end
                end
            end
        end
    end

    return dfs
end

=#

