# Script to post-process the results of the optimization problem

using MyPowerModels
import MyPowerModels: solve_ac_opf_H_min
using PowerModelsAnalytics
using JuMP
using Gurobi
using Plots
using Ipopt
using Juniper
using JLD


data = JLD.load("data.jld")["data"]
options = JLD.load("options.jld")["options"]
results = JLD.load("results.jld")["results"]

gen_data = data["gen"]
bus_data = data["bus"]

branch_data = data["branch"]
load_data = data["load"]
baseMVA = data["baseMVA"]

sol_gen_data = results["gen"]
sol_bus_data = results["bus"]

f_options = options["f"]
alpha = f_options["alpha_factor"]
calc_delta_P = f_options["calc_delta_P"]
rocof = f_options["rocof"]
bus = f_options["bus"]
area = f_options["area"]
weighted_area = f_options["weighted_area"]

for (gen_id, gen) in gen_data
    H = data["gen"][gen_id]["H"]
    pmin = data["gen"][gen_id]["pmin"]
    pmax = data["gen"][gen_id]["pmax"]
    gen_bus = data["gen"][gen_id]["gen_bus"]

    if !haskey(sol_gen_data, gen_id)
        sol_gen_data[gen_id] = Dict()
    end

    sol_gen_data[gen_id]["H"] = H
    sol_gen_data[gen_id]["pmin"] = pmin
    sol_gen_data[gen_id]["pmax"] = pmax
    sol_gen_data[gen_id]["gen_bus"] = gen_bus
end

for (bus_id, bus) in bus_data
   
    area = data["bus"][bus_id]["area"]
    sol_bus_data[bus_id]["area"] = area

end

for i in keys(sol_gen_data)
    if sol_gen_data[i]["gen_status"] < 1e-2
        sol_gen_data[i]["gen_status"] = 0
    end
end

H_sys = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in eachindex(sol_gen_data)) / sum(sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in eachindex(sol_gen_data))
println("H_sys = ", H_sys)

delta_P = data["delta_P"]
println("delta_P = ", delta_P)

H_min = data["H_min"]
println("H_min = ", H_min)

# Initialisierung der Dictionaries für Berechnungen
f0 = 50
H_min_bus = Dict()
H_bus = Dict()
P_gen_bus = Dict()
P_load_bus = Dict()
delta_p_bus = Dict()
#=
# berechnung für weighted_area
if f_options["weighted_area"] == "load"

    areas = unique([bus_data[j]["area"] for j in keys(bus_data)])
    sum_H_area_weighted = 0
    for area in areas

        W_v = Dict()
        H_area = Dict()
        load_sum_area = Dict()
        gens_in_area = [i for i in eachindex(gen_data) if sol_bus_data[gen_data[i]["gen_bus"]]["area"] == area]
        loads_in_area = [i for i in eachindex(load_data) if bus_data[load_data[i]["load_bus"]]["area"] == area]
        load_sum_area[area] = sum(load_data[i]["pd"] for i in loads_in_area; init=0)
        H_area[area] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0) / sum(gen_data[i]["pmax"] * z[i] for i in gens_in_area; init=0)
        W_v[area] = load_sum_area[area] / P_load
        println("W_v von area $area: ", W_v[area])
        sum_H_area_weighted += W_v[area] * H_area[area]
        println("sum_H_area_weighted: ", sum_H_area_weighted)

    end
end
=#
if f_options["bus"] == "true" && f_options["disturbance"] == "large"
    for j in keys(sol_bus_data)

        gens_at_bus = [i for i in keys(sol_gen_data) if string(sol_gen_data[i]["gen_bus"]) == j]
        P_gen_bus[j] = sum(sol_gen_data[i]["pg"] for i in gens_at_bus; init=0)
        P_load_bus[j] = sum(load_data[i]["pd"] for i in keys(load_data) if string(load_data[i]["load_bus"]) == j; init=0)
        delta_p_bus[j] = abs(P_gen_bus[j] - P_load_bus[j])
        H_min_bus[j] = (delta_p_bus[j] * f0) / (P_load_bus[j] * 2 * rocof)
        global H_bus = sum(sol_gen_data[i]["H"] * sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_at_bus; init=0) / sum(sol_gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_at_bus; init=0)

        println("Generatoren am Bus $j: ", gens_at_bus)
        println("P_load_bus $j = ", P_load_bus[j])
        println("P_gen_bus $j = ", P_gen_bus[j])
        println("delta_p_bus $j = ", delta_p_bus[j])
        println("H_min_bus $j = ", H_min_bus[j])
        println("H_bus $j = ", H_bus)

    end
end

# Berechnung für jede Area
H_min_area = Dict()
H_area = Dict()
P_gen_area = Dict()
P_load_area = Dict()
delta_p_area = Dict()

global areas = unique([bus_data[j]["area"] for j in keys(bus_data)])

if f_options["area"] == "true" && f_options["disturbance"] == "large"
    for r in areas

        gens_in_area = [i for i in keys(sol_gen_data) if haskey(sol_bus_data, string(sol_gen_data[i]["gen_bus"])) && sol_bus_data[string(sol_gen_data[i]["gen_bus"])]["area"] == r]
        buses_in_area = [i for i in keys(bus_data) if bus_data[i]["area"] == r]

        P_load_area[r] = sum(load_data[i]["pd"] for i in keys(load_data) if string(load_data[i]["load_bus"]) in buses_in_area; init=0)
        P_gen_area[r] = sum(sol_gen_data[i]["pg"] for i in gens_in_area; init=0)
        delta_p_area[r] = abs(P_gen_area[r] - P_load_area[r])
        H_min_area[r] = (delta_p_area[r] * f0) / (P_load_area[r] * 2 * rocof)
        H_area[r] = sum(gen_data[i]["H"] * gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0) / sum(gen_data[i]["pmax"] * sol_gen_data[i]["gen_status"] for i in gens_in_area; init=0)

        println("Gens in area $r: ", gens_in_area)
        println("Buses in area $r: ", buses_in_area)
        println("P_load_area $r = ", P_load_area[r])
        println("P_gen_area $r = ", P_gen_area[r])
        println("delta_p_area $r = ", delta_p_area[r])
        println("H_area_min $r = ", H_min_area[r])
        println("H_area $r = ", H_area[r])

    end
end

# Berechnung für gewichtete buses
buses = collect(keys(H_min_bus))
h_min_values = [H_min_bus[bus] for bus in buses]

bar(buses, h_min_values, title="H_min für verschiedene Busse", xlabel="Bus", ylabel="H_min", legend=false)

#=
# Berechnung für gewichtete areas
areas = collect(keys(H_min_area))
h_min_values = [H_min_area[area] for area in areas]

bar(areas, h_min_values, title="H_min für verschiedene Areas", xlabel="Area", ylabel="H_min", legend=false)

=#