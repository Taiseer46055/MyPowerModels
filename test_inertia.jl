# Author: Taiseer Alhaj Ali
# Created: 2024-06-15


using MyPowerModels
using PowerModelsAnalytics
using JuMP
using Gurobi
using Plots
using Ipopt
using Juniper
using JLD
using JLD2
using GLPK
using SCIP
using Cbc

dir = ".\\test\\data\\matpower\\multi_nw"
case_name = "case2"


options = Dict( 
    "f" => Dict(
        "inertia_constraint"=> "true", # "true", "false"
        "system" => "false", # "true", "false"
        "disturbance" => "large", # "small", "large" 
        "weighted_area" => "load", # "load", "equal", "none"
        "area" => "true", # "true", "false"
        "bus" => "false", # "true", "false"
        "calc_delta_P" => [1,250], # "internal" or array of [gen_id, delta_P in MW]
        "alpha_factor" => 0.1, # [0, 1]
        "rocof" => 1.0,
    ), 
    "v" => Dict(
        "voltage_constraint" => "false",
        "reactive_power_limit" => "false",
        "max_rvc" => "false",    # rapid voltage change
        "max_pas" => "false",    # phase angle shift
        "area" => "false",
        "weighted_area" => "false", # Abhängig von den Leistungsflüssen und größe der verfügbaren Erzeugern in den betroffenen Region. Je grösser LF umso mehr Spannungsunterstützung erforderlich
    )
)

configurations = [
    ("Fall_1", Dict("inertia_constraint" => "false", "system" => "false", "disturbance" => "small", "weighted_area" => "none", "area" => "false", "bus" => "false")),
    ("Fall_2", Dict("inertia_constraint" => "true", "system" => "true", "disturbance" => "small", "weighted_area" => "none", "area" => "false", "bus" => "false")),
    ("Fall_3", Dict("inertia_constraint" => "true", "system" => "false", "disturbance" => "small", "weighted_area" => "load", "area" => "false", "bus" => "false")),
    ("Fall_4", Dict("inertia_constraint" => "true", "system" => "false", "disturbance" => "small", "weighted_area" => "equal", "area" => "false", "bus" => "false")),
    ("Fall_5", Dict("inertia_constraint" => "true", "system" => "false", "disturbance" => "large", "weighted_area" => "none", "area" => "true", "bus" => "false"))
    #("Fall_6", Dict("inertia_constraint" => "true", "system" => "false", "disturbance" => "large", "weighted_area" => "none", "area" => "false", "bus" => "true"))
]

# minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>2), "log_levels"=>[:all])
m = Model()
function update_options!(base_options, case_options)
    for (key, case_value) in case_options
        if haskey(base_options["f"], key)
            base_options["f"][key] = case_value
        end
    end
end

function extract_and_parse_number_from_filename(filename)
    numeric_part_match = match(r"(\d+).m$", filename)
    if numeric_part_match !== nothing
        parsed_number = tryparse(Int, numeric_part_match[1])
        if isnothing(parsed_number)
            println("Warning: The numeric part found could not be converted to an integer.")
        end
        return parsed_number
    else
        println("Warning: No numeric part found in the filename.")
        return nothing
    end
end


function load_multinetwork_data(case_name, dir)
    filenames = readdir(dir)
    pattern = "$(case_name)_"
    filtered_filenames = filter(f -> occursin(pattern, f) && endswith(f, ".m"), filenames)


    sorted_filenames = sort(filtered_filenames, by=f -> extract_and_parse_number_from_filename(f))

    mn_data = Dict("nw" => Dict(), "per_unit" => true, "multinetwork" => true)
    for filename in sorted_filenames
        full_path = joinpath(dir, filename)
        parsed_id = extract_and_parse_number_from_filename(filename) 
        if isfile(full_path) && parsed_id !== nothing
            data = MyPowerModels.parse_file(full_path)
            mn_data["nw"][parsed_id] = data
        else
            println("Datei $filename nicht gefunden oder ID ist ungültig.")
        end
    end

    return mn_data
end


function main()

    mn_data = load_multinetwork_data(case_name, dir)

    string_nw_keys = Dict(string(k) => v for (k, v) in pairs(mn_data["nw"]))
    mn_data["nw"] = string_nw_keys

    for (case_label, case_options) in configurations
        update_options!(options, case_options)

        result_mn = solve_mn_opf_with_inertia_and_generator_expansion(mn_data, DCPPowerModel, Gurobi.Optimizer, options, jump_model=Model(); multinetwork=true)
        results = result_mn["solution"]
        
        save_results(case_name, results, results_filename, data_filename, options_filename, mn_data)
        script_path = joinpath(pwd(), "post_processing.jl")
        println("Skriptpfad: ", script_path)

        if isfile(script_path)
            println("Running post_processing.jl...")
            include(script_path)
        else
            println("post_processing.jl not found. Skipping post-processing..")
        end
    end
end

# Save the results
results_filename = Dict()
data_filename = Dict()
options_filename = Dict()

function save_results(case_name, results, results_filename, data_filename, options_filename, mn_data)
    if "nw" in keys(results)
        println(" save the result for the case $case_name")

        results_filename[case_name] = "results\\multi_network_results\\results_$case_name.jld2"
        data_filename[case_name] = "results\\multi_network_results\\data_$case_name.jld2"
        options_filename[case_name] = "results\\multi_network_results\\options_$case_name.jld2"

        network_keys = keys(results["nw"])
        sorted_results = Dict(string(key) => results["nw"][string(key)] for key in network_keys)

        JLD2.save(results_filename[case_name], "results", sorted_results)
        JLD2.save(data_filename[case_name], "data", mn_data)
        JLD2.save(options_filename[case_name], "options", options)
    else

        results_filename[case_name] = "results\\single_network_results\\results_$case_name.jld2"
        data_filename[case_name] = "results\\single_network_results\\data_$case_name.jld2"
        options_filename[case_name] = "results\\single_network_results\\options_$case_name.jld2"

        JLD2.save(results_filename[case_name], "results", results)
        JLD2.save(data_filename[case_name], "data", data)
        JLD2.save(options_filename[case_name], "options", options)
    end
end


main()

# function update_load_data!(mn_data, last_profile)
#     network_keys = sort(collect(keys(mn_data["nw"])))
#     for (n, key) in enumerate(network_keys)
#         network = mn_data["nw"][key]
#         for load in values(network["load"])
#              load["pd"] *= last_profile[n]
#         end
#     end
# end
# update_load_data!(mn_data, last_profile)

    # case_name = "case2"
    # grid = ".\\test\\data\\matpower\\$case_name.m"
    # data = MyPowerModels.parse_file(grid)
    # mn_data = MyPowerModels.replicate(data, 24)
    # last_profile = [0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 1.0, 1.0, 0.8, 0.8, 1.0, 1.1, 1.1, 1.2, 1.3, 1.4, 1.2, 1.0, 1.0, 0.9, 0.9, 0.8, 0.8, 0.7]  # Beispiel für Lastprofile

    # update_load_data!(mn_data, last_profile)

# save_results(case_name, results, results_filename, data_filename, options_filename)

#=
case_name = "case2"
case_name
grid = ".\\test\\data\\matpower\\$case_name.m"
data = MyPowerModels.parse_file(grid)

mn_data = MyPowerModels.replicate(data, 24)
last_profile = [0.9, 1.0, 1.2]
#last_profile = [2.8, 2.9, 2.8]
last_profile = [0.7, 0.7, 0.7, 0.7, 0.7, 1.1, 1.3, 1.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2, 1.5, 1.5, 1.3, 1.3, 0.9, 0.9, 0.9, 0.9, 0.7]

m = Model()

result_sn = solve_opf_with_inertia(data, DCPPowerModel, Gurobi.Optimizer, options)
results = result_sn["solution"]
result_mn = solve_mn_opf_with_inertia(mn_data, DCPPowerModel, Gurobi.Optimizer , options, multinetwork=true)
results = result_mn["solution"]

result_mn = solve_mn_opf_with_inertia_and_generator_expansion(mn_data, DCPPowerModel,  Gurobi.Optimizer, options, jump_model=m; multinetwork=true)
results = result_mn["solution"]
println(m)
=#


#=
# define potential_generators for generator expansion planning
pot_gens = Dict(
    "Gen_1" => Dict("pg_min" => 50.0, "pg_max" => 5000.0, "qg_min" => -50.0, "qg_max" => 50.0, "investment_cost" => 2000000, "operating_cost" => 10,"H" => 8.0, "num_blocks" => 4),
    "Gen_2" => Dict("pg_min" => 20.0, "pg_max" => 2000.0, "qg_min" => -50.0, "qg_max" => 50.0, "investment_cost" => 1000000, "operating_cost" => 10,"H" => 5.0 , "num_blocks" => 5),

)

function add_gen_exp_data!(mn_data, pot_gens)
    network_keys = sort(collect(keys(mn_data["nw"])))
    for (n, key) in enumerate(network_keys)
        network = mn_data["nw"][key]
        # Create a new dictionary 'gen_exp' if it doesn't exist
        if !haskey(network, "gen_exp")
            network["gen_exp"] = Dict()
        end
        # Add pot_gens to 'gen_exp'
        for (gen_id, gen_attrs) in pot_gens
            if !haskey(network["gen_exp"], gen_id)
                network["gen_exp"][gen_id] = gen_attrs
            end
        end
    end
end

add_gen_exp_data!(mn_data, potential_generators) 
=#

#=
# Plot the network connections

graph = plot_network(data; 
    aggregate_extra_nodes=true,
    node_size_limits=[50, 100], 
    edge_width_limits=[10, 20], 
    label_nodes=true, 
    fontsize=10, 
    plot_size=(800,800),
    plot_dpi=100);


node_labels = Dict(node_id => "Bus: $(get_bus_number(node_id)), Gen: $(data["gen"][node_id])" for node_id in keys(data["bus"]) if haskey(data["gen"], node_id))
node_label_colors = Dict(node_id => "red" for node_id in keys(data["bus"]))

graph = plot_network(data; 
    aggregate_extra_nodes=true,
    node_size_limits=[50, 100], 
    edge_width_limits=[10, 20], 
    label_nodes=true, 
    node_labels=node_labels,
    node_label_colors=node_label_colors,
    fontsize=10, 
    plot_size=(500,500),
    plot_dpi=100);

=#
