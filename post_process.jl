using JLD2
using DataFrames
using CSV
using StatsPlots
gr(size=(400,300))

cases = Dict( "case_2_35" => Dict("constraint"=> "System Inertia", "RE"=> "35%"),
"case_8_35" => Dict("constraint"=> "Region Islanding 100%", "RE"=> "35%"),
"case_1_35" => Dict("constraint"=> "Unit Commitment", "RE"=> "35%"),
"case_8_50" => Dict("constraint"=> "Region Islanding 100%", "RE"=> "50%"),
"case_6_50" => Dict("constraint"=> "Region Islanding 50%", "RE"=> "50%"),
"case_1_50"=> Dict("constraint"=> "Unit Commitment", "RE"=> "50%"),
"case_6_35"=> Dict("constraint"=> "Region Islanding 50%", "RE"=> "35%"),
"case_2_50"=> Dict("constraint"=> "System Inertia", "RE"=> "50%"), ) 



r_df2 = JLD2.load("results_dataframes/case_2_35.jld2")
r_df1 = JLD2.load("results_dataframes/case_1_35.jld2")
r_df6 = JLD2.load("results_dataframes/case_6_35.jld2")

r_df2_50 = JLD2.load("results_dataframes/case_2_50.jld2")
r_df1_50 = JLD2.load("results_dataframes/case_1_50.jld2")

r_df6_50 = JLD2.load("results_dataframes/case_6_50.jld2")


carrier_names = CSV.read("carrier_indices.csv",DataFrame)
# r_df = r_df["dfs"]
global gen_df = DataFrame()
println(gen_df)
for c in keys(cases)
    println( c)
    r_df = JLD2.load("results_dataframes/$c.jld2")
    # println(r_df)
    gen_c = r_df["results"]["gen"]
    gen_t_c_df=DataFrame()
    gen_c_df = DataFrame()
    kk=0
    CSV.write("pg_t/$c.csv",(r_df["results"]["gen_t"]["pg"]))
    for k in keys(gen_c)
        println(k)
        gen_c[k] = rename(gen_c[k], :value =>k)
        
        if kk==0
            gen_c_df[!,:id]=gen_c[k][!,:id]
            gen_c_df[!,"constraint"] .= cases[c]["constraint"]
            gen_c_df[!,"RE"] .= cases[c]["RE"]
            gen_c_df[!,k] = gen_c[k][!,k]
            kk = 1
        else
            gen_c_df = innerjoin(gen_c_df, gen_c[k], on= :id)
        end

        
    end 
    # println(gen_df)
    # return gen_c_df
    global gen_df =vcat(gen_df,gen_c_df)
    # println(gen_df)
end








gen_df[!,:carrier_names]=combine(gen_df,:carrier => v -> carrier_names[v.+1,3]).carrier_function;
gen_df[!,:carrier_names_h]=combine(gen_df,:carrier => v -> carrier_names[v.+1,1]).carrier_function;
gen_df[!,:pExp]= (gen_df[!,:nE]-gen_df[!,:n0]).*gen_df[!,:P_b_nom]/1e3
g_gen_df = groupby(gen_df, [:constraint, :RE, :carrier_names_h])
g_exp_df = filter(row->row.pExp>0, gen_df)
exp_df = combine(g_gen_df, :pExp => sum)
exp_df = filter(row->row.pExp_sum>0,exp_df)

CSV.write("gen_df.csv", gen_df)
CSV.write("exp_df.csv", exp_df)

# g_exp_df35 = filter(row->row.RE=="35%", g_exp_df)
# g_exp_df70 = filter(row->row.RE=="50%", g_exp_df)

# exp_df35 = filter(row->row.RE=="35%",exp_df)


# @df exp_df35 plot(x= :RE, y= :pExp_sum, color=:carrier_names)

# exp_df70 = filter(row->row.RE=="50%",exp_df)


# @df exp_df70 plot(x= :RE, y= :pExp_sum, color=:carrier_names)

# # return gen_df
# # gen = r_df["gen"]
# # gen_df = innerjoin(innerjoin(rename(gen["nE"],[:id, :nE]),rename(gen["n0"],[:id, :n0]),on= :id),rename(gen["P_b_nom"],[:id,:pb]),on= :id)