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



#=
mn_data = Dict("nw" => Dict(), "multinetwork" => true, "per_unit" => true)

for i in 1:3
    filename = ".\\test\\data\\matpower\\case9_"  * lpad(i, 2, "0") * ".m"
    network_value = MyPowerModels.parse_file(filename)
    mn_data["nw"][string(i)] = network_value
end
=#

case_name = "case2"
case_name
netz = ".\\test\\data\\matpower\\$case_name.m"
data = MyPowerModels.parse_file(netz)

mn_data = MyPowerModels.replicate(data, 24)
#last_profile = [0.8, 1.5, 1.2]
last_profile = [0.7, 0.7, 0.7, 0.7, 0.7, 1.1, 1.3, 1.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2, 1.5, 1.5, 1.3, 1.3, 0.9, 0.9, 0.9, 0.9, 0.7]

function update_load_data!(mn_data, last_profile)
    network_keys = sort(collect(keys(mn_data["nw"])))
    for (n, key) in enumerate(network_keys)
        network = mn_data["nw"][key]
        for load in values(network["load"])
             load["pd"] *= last_profile[n]
        end
    end
end



update_load_data!(mn_data, last_profile)

minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>2), "log_levels"=>[:all])
options = Dict( 
    "f" => Dict(
        "inertia_constraint"=> "true", 
        "system" => "true",
        "disturbance" => "large", # "small", "large" 
        "weighted_area" => "none", # "load", "equal", "none"
        "area" => "false", 
        "bus" => "true", 
        "calc_delta_P" => [1,250], # "internal" or array of [gen_id, delta_P in MW]
        "alpha_factor" => 0.1, # [0, 1]
        "rocof" => 1.0
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

# result_sn = solve_ac_opf_with_inertia(data, DCPPowerModel, Gurobi.Optimizer, options)
# results = result_sn["solution"]
result_mn = solve_mn_ac_opf_with_inertia(mn_data, DCPPowerModel, Gurobi.Optimizer , options, multinetwork=true)
results = result_mn["solution"]



# Save the results
results_filename = Dict()
data_filename = Dict()
options_filename = Dict()


function save_results(case_name, results, results_filename, data_filename, options_filename)
    if "nw" in keys(results)
        println(" save the result for the case $case_name")

        results_filename[case_name] = "results\\multi_network_results\\results_$case_name.jld2"
        data_filename[case_name] = "results\\multi_network_results\\data_$case_name.jld2"
        options_filename[case_name] = "results\\multi_network_results\\options_$case_name.jld2"

        network_keys = sort([parse(Int, key) for key in keys(results["nw"])])
        sorted_results = Dict(string(key) => results["nw"][string(key)] for key in network_keys)

        JLD2.save(results_filename[case_name], "results", sorted_results)
        JLD2.save(data_filename[case_name], "data", mn_data)
        JLD2.save(options_filename[case_name], "options", options)
    else
        local gen_data = results["gen"]

        pg_values = Dict()
        for (gen_id, gen_attrs) in gen_data
            pg_values[gen_id] = gen_attrs["pg"]
        end
        println("Pg-Werte der Generatoren mit H_min:")
        for (gen_id, pg_value) in pg_values
            println("Generator $gen_id: Pg = $pg_value")
        end

        results_filename[case_name] = "results\\single_network_results\\results_$case_name.jld2"
        data_filename[case_name] = "results\\single_network_results\\data_$case_name.jld2"
        options_filename[case_name] = "results\\single_network_results\\options_$case_name.jld2"

        JLD2.save(results_filename[case_name], "results", results)
        JLD2.save(data_filename[case_name], "data", data)
        JLD2.save(options_filename[case_name], "options", options)
    end
end

save_results(case_name, results, results_filename, data_filename, options_filename)
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
