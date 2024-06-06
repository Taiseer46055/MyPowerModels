# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# %%

# os.path.join(os.getcwd(), "pg_t\case1_35.csv") 

# %%
gen_t35_1 =pd.read_csv(os.path.join(os.getcwd(), "pg_t\case_1_35.csv"))
gen_t35_2 =pd.read_csv("./pg_t/case_2_35.csv")
w =pd.read_csv("./pg_t/weights.csv",delimiter=";")
c_df =pd.read_csv("carrier_indices.csv",delimiter=",")
gen_t35_6 =pd.read_csv("./pg_t/case_6_35.csv")
gen_t35_8 =pd.read_csv("./pg_t/case_8_35.csv")

# gen_t70_1 =pd.read_csv("./pg_t/case1_70.csv")
# gen_t70_2 =pd.read_csv("./pg_t/case2_70.csv")
# gen_t70_6 =pd.read_csv("./pg_t/case6_70.csv")
# gen_t70_8 =pd.read_csv("./pg_t/case8_70.csv")
gen_t35_1 = gen_t35_1.multiply(w.w,axis=0)
gen_t35_2 = gen_t35_2.multiply(w.w,axis=0)
gen_t35_6 = gen_t35_6.multiply(w.w,axis=0)
gen_t35_8 = gen_t35_8.multiply(w.w,axis=0)
gen_t35_1 = gen_t35_1.drop(columns="nw")
gen_t35_2 = gen_t35_2.drop(columns="nw")
gen_t35_6 = gen_t35_6.drop(columns="nw")
gen_t35_8 = gen_t35_8.drop(columns="nw")

gen_t35_1.columns = gen_t35_1.columns.astype(int)
gen_t35_2.columns = gen_t35_2.columns.astype(int)
gen_t35_6.columns = gen_t35_6.columns.astype(int)
gen_t35_8.columns = gen_t35_8.columns.astype(int)

gen_df = pd.read_csv("gen_df.csv",delimiter=",")
gen_df ['agg_tech'] = gen_df.carrier.apply(lambda x: c_df.loc[int(x),'agg_tech'])
# gen_df_red = gen_df.loc[(gen_df.constraint=="System Inertia") & (gen_df.RE=="35%")]
# gen_df_red['pMax'] = gen_df_red["P_b_nom"] * gen_df_red["nE"]
# gen_df_red['pSum1'] =gen_df_red.id.apply(lambda x: gen_t35_1.loc[:,x].sum())
# gen_df_red['pSum2'] =gen_df_red.id.apply(lambda x: gen_t35_2.loc[:,x].sum())
# gen_df_red['pSum6'] =gen_df_red.id.apply(lambda x: gen_t35_6.loc[:,x].sum())
# gen_df_red['pSum8'] =gen_df_red.id.apply(lambda x: gen_t35_8.loc[:,x].sum())

# gen_t_agg = gen_df_red.groupby("carrier")[['pSum1','pSum2', 'pSum6', 'pSum8','pMax']].sum()*1e-3


# gen_df_red = gen_df.loc[(gen_df.constraint=="System Inertia") & (gen_df.RE=="35%")]
gen_df['pMax'] = gen_df["P_b_nom"] * gen_df["nE"]
gen_df['pMax0'] = gen_df["P_b_nom"] * gen_df["n0"]
gen_df['pSum'] = 0
gen_df.loc[(gen_df.constraint=="Unit Commitment") & (gen_df.RE=="35%"),'pSum'] =gen_df.loc[(gen_df.constraint=="Unit Commitment") & (gen_df.RE=="35%"),'id'].apply(lambda x: gen_t35_1.loc[:,x].sum())
gen_df.loc[(gen_df.constraint=="System Inertia") & (gen_df.RE=="35%"),'pSum'] =gen_df.loc[(gen_df.constraint=="System Inertia") & (gen_df.RE=="35%"),'id'].apply(lambda x: gen_t35_2.loc[:,x].sum())

gen_df.loc[(gen_df.constraint=="Region Islanding 50%") & (gen_df.RE=="35%"),'pSum'] =gen_df.loc[(gen_df.constraint=="Region Islanding 50%") & (gen_df.RE=="35%"),'id'].apply(lambda x: gen_t35_6.loc[:,x].sum())

gen_df.loc[(gen_df.constraint=="Region Islanding 100%") & (gen_df.RE=="35%"),'pSum'] =gen_df.loc[(gen_df.constraint=="Region Islanding 100%") & (gen_df.RE=="35%"),'id'].apply(lambda x: gen_t35_8.loc[:,x].sum())


gen_df_sum = gen_df.groupby(['RE','constraint','agg_tech'])[['pMax0','pMax','pSum']].sum()

df_pMax0 = gen_df_sum.loc['35%','pMax0'].unstack()


df_pMax = gen_df_sum.loc['35%','pMax'].unstack()

df_Mix = gen_df_sum.loc['35%','pSum'].unstack()


df_pExp = ((df_pMax-df_pMax0)*1e-3)

df_pExp_rel = df_pExp.div(df_pExp.sum(axis=1),axis=0)*100

df_Mix_rel = df_Mix.div(df_Mix.sum(axis=1),axis=0)*100


# gen_df_red['pSum2'] =gen_df_red.id.apply(lambda x: gen_t35_2.loc[:,x].sum())
# gen_df_red['pSum6'] =gen_df_red.id.apply(lambda x: gen_t35_6.loc[:,x].sum())
# gen_df_red['pSum8'] =gen_df_red.id.apply(lambda x: gen_t35_8.loc[:,x].sum())

# gen_t_agg = gen_df_red.groupby("carrier")[['pSum1','pSum2', 'pSum6', 'pSum8','pMax']].sum()*1e-3


## gen_t_agg = gen_df_red.groupby("carrier_names")[['pSum1','pSum2', 'pSum6', 'pSum8','pMax']].sum()*1e-3

gen_t_agg = gen_t_agg.rename(columns={'pSum1':'UC','pSum2':'System Inertia','pSum6': '50% Islanding','pSum8':'100% Islanding','pMax':'Installed Capacity'})

gen_t_agg ["tech"] =c_df.tech 
gen_t_agg ["agg_tech"] =  c_df.agg_tech

gen_t_agg_red = gen_t_agg.groupby("agg_tech").sum()
# gen_t_agg["UC"] = gen_t_agg["UC"]/gen_t_agg["pMax"]
# gen_t_agg["System Inertia"] = gen_t_agg["System Inertia"]/gen_t_agg["UC"]*100-100
gen_t_agg_red.loc[["Gas","Coal","PV","Wind"],["UC","System Inertia","50% Islanding","100% Islanding"]].map(round, ndigits=1)
# gen_t_agg["50% Islanding"] = gen_t_agg["50% Islanding"]/gen_t_agg["UC"]*100-100
# gen_t_agg["100% Islanding"] = gen_t_agg["100% Islanding"]/gen_t_agg["UC"]*100-100


# %%
gen_exp = pd.read_csv("exp_df2.csv",delimiter=";")
# gen_df = pd.read_csv("gen_df.csv",delimiter=",")
# gen_df_red = gen_df.loc[(gen_df.constraint=="System Inertia") & (gen_df.RE=="35%")]


# %%
gen_t35_1.sum()

# %%
gen_exp35=gen_exp.loc[gen_exp.RE=="35%"]
g_piv35= gen_exp35.pivot(index=["RE","constraint"],columns="carrier_names",values=["pExp_sum"])
# p_piv.rename
g_piv35 = g_piv35.droplevel(0)
g_piv35 = g_piv35.droplevel(0,axis=1)
g_piv35
# g_piv35= g_piv35/g_piv35.loc["Unit Commitment",:]
g_piv35

gen_exp70=gen_exp.loc[gen_exp.RE=="70%"]
g_piv70= gen_exp70.pivot(index=["RE","constraint"],columns="carrier_names",values=["pExp_sum"])
# p_piv.rename
g_piv70 = g_piv70.droplevel(0)
g_piv70 = g_piv70.droplevel(0,axis=1)
# g_piv70= g_piv70/g_piv70.loc["Unit Commitment",:]
g_piv35.sort_index(ascending=False,inplace=True)
g_piv35.sort_index(axis=1,ascending=False,inplace=True)

g_piv70.sort_index(ascending=False,inplace=True)
g_piv70.sort_index(axis=1,ascending=False,inplace=True)
g_piv35.columns.name=''
g_piv70.columns.name=''
print(g_piv35)
g_piv35.to_latex()
# g_piv70

# %%
g_piv35.columns.name=''

# %%
g_piv35.plot.bar(stacked=True,subplots=False,logy=False,rot=45,legend=["Gas", "PV","Wind"],ylabel="Capacity Expansion in GW")

# %%
g_piv35.plot.bar(stacked=False,subplots=False,logy=False,rot=45,legend=["Gas", "PV","Wind"])

# %%



