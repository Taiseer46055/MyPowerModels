using MyPowerModels
using PowerModelsAnalytics
using JuMP
using Gurobi
using Plots
using Ipopt
using Juniper
using JLD

data = Dict("nw" => Dict(), "multinetwork" => true, "per_unit" => true)

for i in 1:3
    filename = ".\\test\\data\\matpower\\case9_" * lpad(i, 2, "0") * ".m"
    network_value = MyPowerModels.parse_file(filename)
    data["nw"][string(i)] = network_value
end

# netz = ".\\test\\data\\matpower\\case9.m"
# data = MyPowerModels.parse_file(netz)
# mn_data = MyPowerModels.parse_file(netz)
# multi_network_data = MyPowerModels.replicate(mn_data, 3)


minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>2), "log_levels"=>[:all])
options = Dict( 
    "f" => Dict(
        "inertia_constraint"=> "true", 
        "system" => "true",
        "disturbance" => "small", # "small", "large" 
        "weighted_area" => "none", # "load", "equal", "none"
        "area" => "true", 
        "bus" => "false", 
        "calc_delta_P" => [1,100], # "internal" or array of [gen_id, delta_P in MW]
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

# result_1 = solve_ac_opf_with_inertia(sn_data, ACPPowerModel, minlp_solver, options)
result_1 = solve_mn_ac_opf_with_inertia(data, ACPPowerModel, minlp_solver, options, multinetwork=true)
results = result_1["solution"]
if "multiinfrastructure" in keys(results["nw"])
    gen_data = results["nw"]["multiinfrastructure"]["gen"]
    bus_data = results["nw"]["multiinfrastructure"]["bus"]
    branch_data = results["nw"]["multiinfrastructure"]["branch"]
end

results_filename = Dict()
data_filename = Dict()
options_filename = Dict()
for n in keys(results["nw"])
    println("Ergebnisse für Netzwerk $n:")

    gen_data = results["nw"][n]["gen"]
    bus_data = results["nw"][n]["bus"]
    branch_data = results["nw"][n]["branch"]

    pg_values = Dict()
    for (gen_id, gen_attrs) in gen_data
        pg_values[gen_id] = gen_attrs["pg"]
    end
    println("Pg-Werte der Generatoren mit H_min:")
    for (gen_id, pg_value) in pg_values
        println("Generator $gen_id: Pg = $pg_value")
    end

    results_filename[n] = "results_$n.jld"
    data_filename[n] = "data_$n.jld"
    options_filename[n] = "options_$n.jld"

    JLD.save(results_filename[n], "results", results["nw"][n])
    JLD.save(data_filename[n], "data", data["nw"][n])
    JLD.save(options_filename[n], "options", options)
end
#=
result_1["solution"]["gen"]

gen_data = result_1["solution"]["gen"]
bus_data = result_1["solution"]["bus"]
branch_data = result_1["solution"]["branch"]

pg_values = Dict()
for (gen_id, gen_attrs) in gen_data
    pg_values[gen_id] = gen_attrs["pg"]
end
println("Pg-Werte der Generatoren mit H_min:")
for (gen_id, pg_value) in pg_values
    println("Generator $gen_id: Pg = $pg_value")
end

results = result_1["solution"]

results_filename = "results.jld"
data_filename = "data.jld"
options_filename = "options.jld"

JLD.save(results_filename, "results", results)
JLD.save(data_filename, "data", data)
JLD.save(options_filename, "options", options)


=#


# Plot the network connections


#=

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






######  print_level: Controls the output of the optimizer #######

print_level = 0: Keine Ausgabe
print_level = 1: Nur Fehlermeldungen
print_level = 2: Nur geringfügige Ausgabe bei jedem Iterationsschritt
print_level = 3: Mehr detaillierte Ausgabe
print_level = 4: Noch detailliertere Ausgabe
print_level = 5: Maximale Ausgabe

logLevel = 0: Keine Ausgabe
logLevel = 1: Nur Fehlermeldungen
logLevel = 2: Nur geringfügige Ausgabe bei jedem Iterationsschritt
logLevel = 3: Mehr detaillierte Ausgabe
logLevel = "all" Alle verfügbaren Ausgaben



#minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-4, "print_level"=>0), "log_levels"=>[])

#print(result_1["solution"]["gen"])

gen_data = result_1["solution"]["gen"]

pg_values = Dict()
for (gen_id, gen_attrs) in gen_data
    pg_values[gen_id] = gen_attrs["pg"]
end
println("Pg-Werte der Generatoren mit H_min:")
for (gen_id, pg_value) in pg_values
    println("Generator $gen_id: Pg = $pg_value")
end



# Plot the network connections

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



result = solve_ac_opf(data, minlp_solver) #Ipopt.Optimizer)
data = MyPowerModels.parse_file(netz)
result_1 = solve_ac_opf_H_min(data, minlp_solver, 1,0.0, 1.0)



gen_data = result["solution"]["gen"]


pg_values = Dict()
for (gen_id, gen_attrs) in gen_data
    pg_values[gen_id] = gen_attrs["pg"]
end
println("Pg-Werte der Generatoren:")
for (gen_id, pg_value) in pg_values
    println("Generator $gen_id: Pg = $pg_value")
end


gen_data = result_1["solution"]["gen"]

pg_values = Dict()
for (gen_id, gen_attrs) in gen_data
    pg_values[gen_id] = gen_attrs["pg"]
end
println("Pg-Werte der Generatoren mit H_imn:")
for (gen_id, pg_value) in pg_values
    println("Generator $gen_id: Pg = $pg_value")
end

gen_data
=#
