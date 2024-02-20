using MyPowerModels
using PowerModelsAnalytics
using JuMP
using Gurobi
using Plots
using Ipopt
using Juniper
using JLD


netz = ".\\test\\data\\matpower\\case2.m"
#netz = "C:\\Users\\Taiseer\\Projekt\\pmProjekt\\MyPM\\MyPowerModels\\test\\data\\matpower\\case24.m"
data = MyPowerModels.parse_file(netz)


minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>2), "log_levels"=>[:all])
m = Model()

options = Dict( 
    "f" => Dict(
        "inertia_constraint"=> "true", 
        "system" => "true",
        "disturbance" => "small", # "small", "large" 
        "weighted_area" => "load", # "load", "equal", "none"
        "area" => "false", 
        "bus" => "false", 
        "calc_delta_P" => [1, 100], # "internal" or array of [gen_id, delta_P in MW]
        "alpha_factor" => 0.5, # [0, 1]
        "rocof" => 1.0
    ), 
    "v" => Dict()
)
result_1 = solve_ac_opf_H_min(data, minlp_solver, options)

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

JLD.save("results.jld", "results", results)
JLD.save("data.jld", "data", data)
JLD.save("options.jld", "options", options)





#=

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
=#


#minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-4, "print_level"=>0), "log_levels"=>[])

#print(result_1["solution"]["gen"])
#=
gen_data = result_1["solution"]["gen"]

pg_values = Dict()
for (gen_id, gen_attrs) in gen_data
    pg_values[gen_id] = gen_attrs["pg"]
end
println("Pg-Werte der Generatoren mit H_min:")
for (gen_id, pg_value) in pg_values
    println("Generator $gen_id: Pg = $pg_value")
end

=#

# Plot the network connections
#=
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


