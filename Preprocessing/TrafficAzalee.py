"""
Ce script calcule la répartiton du trafic avant et après la mise en place de la rue a sens unique.
Les données de trafic Geomobility, le résau de rue et la répartiton spatioal du trafic (Emis_rues.dat) doivent être copier dans le fichier Preprocessing

"""

mainfile = "Preprocessing/"

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#Import des données de trafic de Geomobility
bx15 = pd.read_csv(mainfile+"/Brussels #15 (Per Time Window).csv", sep=";")
bx16 = pd.read_csv(mainfile+"/Brussels #16 Chazal_Rogier (Per Time Window).csv", sep=";")

for df in [bx15,bx16]:
    df["Site No."] = df["Site No."].str.strip()
    df["total"] = df.loc[:, "M/C":"PSV"].sum(axis=1)
    df["time"] = pd.to_datetime(df["Date"] + " " +df["Start Time"])
    df["hour"] = df["time"].dt.hour


all = pd.concat([bx15,bx16])
all["Campagne"] = "Avant"
all.loc[all.time >= "2022-08-01","Campagne"] = "Apres"
all["dayofweek"] = all.time.dt.dayofweek
all = all.reset_index(drop = True)

#liste des rues et groupement de celles-ci. Les numeros sont les ID des segments de rue
id_per_road = {"E. Cambier (1)": [3323,3324,3325,3326], "E. Cambier (2)":[3121,3122,3327,3328,3329,22327,22328], "E. Cambier (3)":[22083,22084,22085,22086],
               "Pavots (1)":[22324,22325,22326], "Pavots (2)":[22105,22106],
               "Chardons (1)":[3298,3299,3300,3301], "Chardons (2)":[3302,3303,3304],
               "P. Devigne (1)":[3291,3292,3293,3294],"P. Devigne (2)":[3295,3296,3297],"P. Devigne (3)":[3183,3184,3185],
               "Chazal (1)":[3186,3187,3188],"Chazal (2)":[3260,3261,3262],"Chazal (3)":[3266,3267,3268],
               "G. Eisenhower (1)":[2874,2875,2876,2903,2904,2905,2906,2956,2957,2958],"G. Eisenhower (2)":[2817,2818,2819,2820],"G. Eisenhower (3)":[3313,3314,3315,3316,3317,3318,3319,3320,3321,3322,3330,3331,3332],
               "H. Stacquet":[3208,3209,3210,3231,3232,3233,3305,3306,3307],
               "Vandenbussche (1)":[3180,3181,3182], "Vandenbussche (2)": [3279,3280,3281],
               "J. Strobbaerts (1)":[2880, 2881, 2882, 3239, 3240, 3241], "J. Strobbaerts (2)":[3195,3196,3197,3257,3258,3259],
               "J. Impens (1)":[2861,2862,2863,2864,2821,2822,2823,2824],"J. Impens (2)":[2888,2889,2890],
               "F. Binje" :[2842,2843,2844,2858,2859,2860],
               "J. Devreese": [2830,2831,2832,2833,2834,2835],
               "Paquerettes": [2825,2826,2827,2828,2829,2836,2837,2838,2918,2919,2920,2921],
               "F d'amour":[2855,2856,2857,2865,2866,2867],
               "Azalee":[2911,2912,2913,2914,2915,2916,2917,2921,2922,2923,2924,2925,2926,2927,2928,2936,2937,2938,2939],
               "Rogier (1)":[2710,2711,2845,2846,2847,2891,2892,2893,2965,2966,2967,2978,2979,2981,2982,2983,3030,3173,3174,3175,3198,3199,3200],
               "Rogier (2)":[3228,3229,3230,3234,3235,3236,3237,3238,3245,3246,3247,3285,3286,3287,3288,3289,3290,3308,3309,3310,3311,3312,22121,22122,22123]}

road_per_id = {id: road for road, ids in id_per_road.items() for id in ids}
roads_id =  list(road_per_id.keys())

#Liste des des rues par groupe de rues d'évolution avant/après
road_per_camp = {"E. Cambier (1)":  "Chazal (Site 1) - E", "E. Cambier (2)":"Chazal (Site 1) - E", "E. Cambier (3)":"Chazal (Site 1) - E",
               "Pavots (1)":"Chazal (Site 1) - C", "Pavots (2)":"Chazal (Site 1) - C",
               "Chardons (1)":"Chazal (Site 1) - C", "Chardons (2)":"Chazal (Site 1) - C",
               "P. Devigne (1)":"Chazal (Site 1) - C","P. Devigne (2)":"Chazal (Site 1) - C","P. Devigne (3)":"Chazal (Site 1) - C",
               "Chazal (1)":"Chazal (Site 1) - A", "Chazal (2)":"Chazal (Site 1) - D", "Chazal (3)":"Chazal (Site 1) - D",
               "G. Eisenhower (1)":"Chazal (Site 1) - B","G. Eisenhower (2)":"Chazal (Site 1) - B","G. Eisenhower (3)":"Chazal (Site 1) - B",
               "H. Stacquet":"Chazal (Site 1) - C",
               "Vandenbussche (1)":"Chazal (Site 1) - C", "Vandenbussche (2)": "Chazal (Site 1) - C",
               "J. Strobbaerts (1)":"Chazal (Site 1) - C", "J. Strobbaerts (2)":"Chazal (Site 1) - C",
               "J. Impens (1)":"Chazal (Site 1) - C","J. Impens (2)":"Chazal (Site 1) - C",
               "F. Binje" :"Chazal (Site 1) - C",
               "J. Devreese": "Chazal (Site 1) - C",
               "Paquerettes": "Chazal (Site 1) - C",
               "F d'amour":"Chazal (Site 1) - C",
               "Azalee":"Chazal (Site 1) - B",
               "Rogier (1)":"Link Arm B - A", "Rogier (2)":"Link Arm E - A"}


Map_azalee = gpd.read_file(Preprocessing+"Reseau_rues-SIRANE.shp")
Map_azalee.index= Map_azalee.ID.astype(int)
roads_id_toharmonize = list(Map_azalee[~Map_azalee.Rue.isna()].index)
#plt.show()
rue_to_count = {'Chazal (Site 1) - A': 8, 'Chazal (Site 1) - B':10, 'Chazal (Site 1) - C':9, 'Chazal (Site 1) - B+C' : 19,
       'Chazal (Site 1) - D' : 7, 'Chazal (Site 1) - E' : 6, 'Link Arm B - A' : 17 ,'Link Arm E - B': 18,
       'A - 1' : 15, 'A - 2' : 16, 'A - 3' : 12, 'A - 4':14, 'A - 5':13, 'A - 4+5' : 20,'B - 1':11,
       'C - 2' : 4, 'C - 3' : 5, 'D - 2' : 2, 'D - 3' : 1, 'D - 4': 3}
count_to_rue= {val : key for key,val in rue_to_count.items()}

Map_azalee["Groupe rues"] = Map_azalee.Rue.map(count_to_rue)
ax = Map_azalee.loc[Map_azalee.Rue.isna(),:].clip(Map_azalee.loc[roads_id_toharmonize,:].buffer(100).total_bounds).plot(color="lightgrey")
Map_azalee.loc[roads_id_toharmonize,:].plot("Groupe rues", cmap="tab20b",categorical=True, legend=True, ax=ax)
plt.show()

Camp2 = pd.read_excel(Preprocessing+"/Campagne 13mars.xlsx", sheet_name="Feuil1").rename(columns={"Croisement":"Site No.", "Intersection":"Origin"})
Camp2["Campagne"] = "Apres"
data =     (all.pivot_table(index=["Site No.","Origin"     ], columns = "Campagne", values = "total", aggfunc='sum')+\
            all.pivot_table(index=["Site No.","Destination"], columns = "Campagne", values = "total", aggfunc='sum'))/4


Camp2.loc[Camp2["Site No."] == "B","Véhicules"] = Camp2.loc[Camp2["Site No."] == "B","Véhicules"] / Camp2.loc[(Camp2["Site No."] == "B") & (Camp2.Origin == 2),"Véhicules"].values * data.loc[["Chazal (Site 1)"], ["B"] ,:]["Apres"].values
Camp2.loc[Camp2["Site No."] == "A","Véhicules"] = Camp2.loc[Camp2["Site No."] == "A","Véhicules"] / Camp2.loc[(Camp2["Site No."] == "A") & (Camp2.Origin == 3),"Véhicules"].values * Camp2.loc[(Camp2["Site No."] == "B") & (Camp2.Origin == 3),"Véhicules"].values
Camp2.loc[Camp2["Site No."] == "C","Véhicules"] = Camp2.loc[Camp2["Site No."] == "C","Véhicules"] / Camp2.loc[(Camp2["Site No."] == "C") & (Camp2.Origin == 1),"Véhicules"].values * data.loc[["Chazal (Site 1)"], ["E"] ,:]["Apres"].values
Camp2.loc[Camp2["Site No."] == "D","Véhicules"] = Camp2.loc[Camp2["Site No."] == "D","Véhicules"] / Camp2.loc[(Camp2["Site No."] == "D") & (Camp2.Origin == 1),"Véhicules"].values * Camp2.loc[(Camp2["Site No."] == "C") & (Camp2.Origin == 2),"Véhicules"].values
all_camp = pd.concat([data,Camp2.pivot_table(index=["Site No.","Origin"], columns = "Campagne", values="Véhicules", aggfunc='sum')])
all_camp.index = all_camp.reset_index()['Site No.']+ ' - ' + all_camp.reset_index()["Origin"].astype(str)
all_camp.loc["Chazal (Site 1) - B+C",:] = all_camp.loc["Chazal (Site 1) - B",:] + all_camp.loc["Chazal (Site 1) - C",:]
all_camp.loc["A - 4+5",:] = all_camp.loc["A - 4",:] + all_camp.loc["A - 5",:]

Map_azalee["Apres"] = Map_azalee.Rue.map(count_to_rue).map(all_camp.Apres)
Map_azalee["Avant"] = Map_azalee["Apres"] / Map_azalee.ID.astype(int).map(road_per_id).map(evol_per_road)
Map_azalee["Apres/Avant"] = Map_azalee["Apres"]/Map_azalee["Avant"]
Map_azalee.loc[Map_azalee["Apres/Avant"]==1,["Apres/Avant"]] = np.nan

vmin, vmax = (Map_azalee[["Apres","Avant"]].min().min(), Map_azalee[["Apres","Avant"]].max().max())
ax = Map_azalee.loc[Map_azalee.Rue.isna(),:].clip(Map_azalee.loc[roads_id_toharmonize,:].buffer(100).total_bounds).plot(color="lightgrey")
Map_azalee.loc[roads_id_toharmonize,:].plot("Apres", legend=True, ax=ax, vmin=vmin, vmax=vmax)

ax = Map_azalee.loc[Map_azalee.Rue.isna(),:].clip(Map_azalee.loc[roads_id_toharmonize,:].buffer(100).total_bounds).plot(color="lightgrey")
Map_azalee.loc[roads_id_toharmonize,:].plot("Avant", legend=True, ax=ax, vmin=vmin, vmax=vmax)
Map_azalee.describe()

vmin, vmax = (Map_azalee[["Apres/Avant"]].min().min(), Map_azalee[["Apres/Avant"]].max().max())
ax = Map_azalee.loc[Map_azalee["Apres/Avant"].isna(),:].clip(Map_azalee.loc[roads_id_toharmonize,:].buffer(100).total_bounds).plot(color="lightgrey")
Map_azalee.plot("Apres/Avant", legend=True, ax=ax, vmin=vmin, vmax=vmax)
plt.show()
#plt.show()

em_rue = pd.read_csv(Preprocessing+"Emis_rues.dat", na_values=-9999, sep="\t")
em_rue["evol"] = em_rue.Id.astype(int).map(road_per_id).map(evol_per_road)
em_rue["roadcount"] = em_rue.Id.astype(int).map(road_per_id).map(road_per_camp)
em_rue.loc[em_rue["evol"].isna(),"evol"] = 1
em_rue.describe()
em_rue["Q"] = em_rue.Id.astype(int).map(Map_azalee["Apres"])
em_rue["L"] = em_rue.Id.astype(int).map(Map_azalee["length_2"])
E_total_Azalee = em_rue.loc[em_rue.Id.apply(lambda x : x in roads_id_toharmonize), ['NO',"NO2","O3"]].sum()
em_rue_New = em_rue.copy()
for esp in ['NO',"NO2","O3"]:
    em_rue_New.loc[em_rue_New.Id.apply(lambda x : x in roads_id_toharmonize), esp] = ((em_rue.Q * em_rue.L) * E_total_Azalee[esp]/(em_rue.Q * em_rue.L).sum())

em_rue_New[["Id",'NO',"NO2","O3"]].to_csv(Preprocessing+"Emis_rues-modified_Apres.dat", na_rep=-9999, sep="\t", index=False)
em_rue_New.describe()
for element in ["NO","NO2",'O3']:
    em_rue_New[element] = em_rue_New[element] / em_rue_New["evol"]

em_rue_New[["Id",'NO',"NO2","O3"]].to_csv(Preprocessing+"Emis_rues-modified_Avant.dat", na_rep=-9999, sep="\t", index=False)
