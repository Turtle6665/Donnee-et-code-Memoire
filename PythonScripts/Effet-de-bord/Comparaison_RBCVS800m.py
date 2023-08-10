###
#Doit être executé après la partie validation !
#Comparaison des données modélisée sur tout BX avec les données modélisée finale
###

mainfile = "PythonScripts/Effet-de-bord/"

import PySiWrap.PySiWrap as ps
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload as r
import os
import rasterio as rio
r(ps)

resultfiles = "Resultats/"+mainfile.split("/")[-2]+"/"
os.makedirs(resultfiles,exist_ok=True)

extends = pd.read_csv(f"{mainfile}/effet-de-bord_Extends.txt", sep="\t", header=[0,1], index_col=0)


def models(name, extends_id, dateDeb = "13/05", dateFin="19/05"):
    mod = ps.Model(name, pathToSirane="SIRANE")
    #Change results file
    mod.inputs.FICH_DIR_RESUL = f"RESULT_Effet-de-bord/{name}"
    #Change NO to NO2_equivalent
    mod.inputs.B_EMIS_NOEQNO2 = 1
    #change dates
    mod.inputs.DATE_DEB = f"{dateDeb}/2018 00:00:00"
    mod.inputs.DATE_FIN = f"{dateFin}/2018 23:00:00"
    #Change dispertions parameter
    mod.site_disp.LATITUDE = 50.833
    mod.site_disp.Z0D = 1.4
    mod.site_disp.ZDISPL = 10
    mod.site_disp.ALBEDO = 0.15
    mod.site_disp.EMISSIVITE = 0.94
    mod.site_disp.PRIESTLEY_TAYLOR = 0.75
    #Change meteo parameter
    mod.site_meteo.LATITUDE = 50.797
    mod.site_meteo.ALTITUDE = 27
    mod.site_meteo.Z0D = 1.5
    mod.site_meteo.ZDISPL = 15
    mod.site_meteo.ALBEDO = 0.15
    mod.site_meteo.EMISSIVITE = 0.94
    mod.site_meteo.PRIESTLEY_TAYLOR = 0.75
    #Ne modeliser que les NOx et O3
    Esp = ps.read_dat(mod.inputs.FICH_ESPECES, mod_input = mod)
    Esp.loc[["PM",'PM25'],["Flag"]] = 0
    mod.save_dat(Esp, mod.inputs.FICH_ESPECES)
    #ajout d'un recepteur au dessus de celui d'Ixelles (1m au dessus de la rue canyon)
    recept = ps.read_dat(mod.inputs.FICH_RECEPT, mod_input=mod)
    recept.loc["IXL-Canope",:] = {"name":"Ixelles Rooftop", "X": 151115,  "Y": 168294,  "Z": 1,   "Type": 1, "Fichier":"RECEPTEURS/points_vides.dat"}
    mod.save_dat(recept, mod.inputs.FICH_RECEPT)
    #changer l'extend
    mod.save_dat(extends["info-grid-meteo"].T[[extends_id]].T, mod.inputs.FICH_GRD_MET, index=False)
    mod.save_dat(extends["info-grid-sortie"].T[[extends_id]].T, mod.inputs.FICH_GRD_SORTIE, index=False)
    #desactive un les stats:
    mod.inputs.B_CALC_STAT = 0
    mod.inputs.B_CALC_CRITERES = 0
    return mod

mod = models("S1-50m","50")
Receptors = ps.read_dat(mod.inputs.FICH_RECEPT,mod_input=mod)
Receptors["FullName"] = Receptors.name + " (" +Receptors.index+ ")"
#mod.run()
#ps.force_all_files_to_def(mod) #pour si un modèle s'est arreté précédement et a laisser des fichier "-def".

mod_list = []
extends_ls = list(extends.index)
for extends_id in extends_ls:
    for semaine,dates in enumerate([["13/05","19/05"],["23/05","29/05"]]):
        print(f"S{semaine+1}-{extends_id}m")
        mod_list.append(models(f"S{semaine+1}-{extends_id}m", extends_id, dates[0], dates[1]))
        if os.path.isfile(f"{mod_list[-1].pathToSirane}/{mod_list[-1].inputs.FICH_DIR_RESUL}/log.txt"):
            print("Already runned")
        else:
            mod_list[-1].run(silent=True)

out = {} #recuperation des deux semaines modélisées sur tout Bruxelles
for mod in mod_list:
    if "BXm" in mod.name:
        out[mod.name] = ps.Output(mod)

data_FullRBC = pd.concat([out[o].alldata for o in out])
filter_dates = list(set(data_FullRBC.index))
# récuperations des données de chaque recepteurs
extends_local = pd.read_csv("All_Extends.csv", header=[0,1], index_col=[0,1])
data_local_ls = []
for ids in extends_local.index:
    name = ids[0].replace(" ","_")
    print(name)
    mod = ps.Model(name, pathToSirane="Final-Sirane")
    mod.inputs.FICH_DIR_RESUL = f"RESULT_Validation/{name}"
    list_recept = ids[1].split("/")
    out_mod = ps.Output(mod)
    filter_recept = out_mod.alldata.Recept.apply(lambda x : x in list_recept)
    data_local_ls.append(out_mod.alldata[filter_recept])#.loc[filter_dates,:])

data_local = pd.concat(data_local_ls)
a = data_local.groupby("Recept").min()
a[[f"{i}_{m}" for i in ["NO2","NO","O3"] for m in ["Mod","Mes","Res"]]]
all_recept = list(set(data_local.Recept))
#Elever les faux recepteurs du datafull
data_FullRBC = data_FullRBC[data_FullRBC.Recept.apply(lambda x: x in all_recept)]
#elever les index et trier pour avoir le même ordre
data_FullRBC = data_FullRBC.reset_index().sort_values(["Date","Recept"]).set_index(["Date","Recept"])
data_local = data_local.reset_index().sort_values(["Date","Recept"]).set_index(["Date","Recept"])

#comparaison des données
compared_NO2 = pd.concat({"RBC":data_FullRBC.NO2_Mod,"800m":data_local.NO2_Mod}, axis=1)
#compared_NO2 = pd.concat({"RBC":data_FullRBC.NO_Mod,"800m":data_local.NO_Mod}, axis=1)
#compared_NO2 = pd.concat({"RBC":data_FullRBC.O3_Mod,"800m":data_local.O3_Mod}, axis=1)

Type_Recept = {"41R012": "Urbain avec très faible influence du trafic",
               "41B011": "Urbain avec très faible influence du trafic",
               "41MEU1":"Urbain avec faible influence du trafic",
               "41B006":"Urbain avec faible influence du trafic",
               "41R001":"Urbain avec influence modérée du trafic",
               "41B004":"Urbain avec influence modérée du trafic",
               "41WOL1":"Urbain avec influence modérée du trafic",
               "41N043":"Industriel avec influence modérée du trafic",
               "41R002":"Urbain avec forte influence du trafic",
               "41B008":"Urbain avec très forte influence du trafic",
               "41B001":"Urbain avec très forte influence du trafic"}
type_list = list({i:1 for i in Type_Recept.values()}.keys())
Recept_order = list(Type_Recept.keys())

f,ax = plt.subplots(1,layout="constrained")
sns.boxplot(compared_NO2.dropna().set_index(["Name","Recept","Date"]).stack().reset_index().loc[lambda x: x.level_3 != "Ratio"].rename(columns={"level_3":"Zone",0:"NO2_Mod"}), y="Name",x="NO2_Mod",order=Receptors.FullName[Recept_order].values,hue="Zone")
#sns.boxplot(compared_NO2, y="Name",x="Ratio", order=Receptors.FullName[Recept_order].values)
ax.set(xlabel="Concentration en NO₂ modélisés [-]",ylabel="")
plt.savefig(resultfiles+"Boxplot_RBC_VS_800m.png")
plt.close("all")
