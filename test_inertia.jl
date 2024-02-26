using MyPowerModels
using PowerModelsAnalytics
using JuMP
using Gurobi
using Plots
using Ipopt
using Juniper
using JLD


mn_data = Dict("nw" => Dict(), "multinetwork" => true, "per_unit" => true)

for i in 1:1
    filename = ".\\test\\data\\matpower\\case9_"  * lpad(i, 2, "0") * ".m"
    network_value = MyPowerModels.parse_file(filename)
    mn_data["nw"][string(i)] = network_value
end

netz = ".\\test\\data\\matpower\\case2.m"
data = MyPowerModels.parse_file(netz)


minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>2), "log_levels"=>[:all])
options = Dict( 
    "f" => Dict(
        "inertia_constraint"=> "false", 
        "system" => "true",
        "disturbance" => "small", # "small", "large" 
        "weighted_area" => "none", # "load", "equal", "none"
        "area" => "false", 
        "bus" => "false", 
        "calc_delta_P" => [4,200], # "internal" or array of [gen_id, delta_P in MW]
        "alpha_factor" => 0.1, # [0, 1]
        "rocof" => 1.0
    ), 
    "v" => Dict(
        "voltage_constraint" => "true",
        "reactive_power_limit" => "true",
        "max_rvc" => "false",    # rapid voltage change
        "max_pas" => "false",    # phase angle shift
        "area" => "false",
        "weighted_area" => "false", # Abhängig von den Leistungsflüssen und größe der verfügbaren Erzeugern in den betroffenen Region. Je grösser LF umso mehr Spannungsunterstützung erforderlich
    )
)

result_sn = solve_ac_opf_with_inertia(data, ACPPowerModel, minlp_solver, options)
#result_mn = solve_mn_ac_opf_with_inertia(mn_data, ACPPowerModel, minlp_solver, options, multinetwork=true)
results = result_sn["solution"]



function save_results(results)
    if "nw" in keys(results)
        # Speichern Sie die Ergebnisse für Multi Network
        for n in keys(results["nw"])
            println("Ergebnisse für Netzwerk $n:")

            local gen_data = results["nw"][n]["gen"]
            local bus_data = results["nw"][n]["bus"]
            local branch_data = results["nw"][n]["branch"]

            results_filename[n] = "results_$n.jld"
            data_filename[n] = "data_$n.jld"
            options_filename[n] = "options_$n.jld"

            JLD.save(results_filename[n], "results", results["nw"][n])
            JLD.save(data_filename[n], "data", data["nw"][n])
            JLD.save(options_filename[n], "options", options)
        end
    else
        # Speichern Sie die Ergebnisse für Single Network
        local gen_data = results["gen"]
        local bus_data = results["bus"]
        local branch_data = results["branch"]

        pg_values = Dict()
        for (gen_id, gen_attrs) in gen_data
            pg_values[gen_id] = gen_attrs["pg"]
        end
        println("Pg-Werte der Generatoren mit H_min:")
        for (gen_id, pg_value) in pg_values
            println("Generator $gen_id: Pg = $pg_value")
        end


        results_filename = "results.jld"
        data_filename = "data.jld"
        options_filename = "options.jld"

        JLD.save(results_filename, "results", results)
        JLD.save(data_filename, "data", data)
        JLD.save(options_filename, "options", options)
    end
end


save_results(results)


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
