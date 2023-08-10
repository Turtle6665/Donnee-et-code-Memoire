"""
Ce script regarde l'incertitude liée aux paramtre de site météo et de dispertion.

#"""

mainfile = "PythonScripts/Incertitude/"

import PYSIRANE as ps
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload as r
import os
r(ps)

resultfiles = "Resultats/"+mainfile.split("/")[-2]+"/"
os.makedirs(resultfiles,exist_ok=True)

extends = pd.read_csv("effet-de-bord_Extends.txt", sep="\t", header=[0,1], index_col=0)

def models(name, Z0D_disp, ZDISPL_disp, ALBEDO, EMISSIVITE, PRIESTLEY_TAYLOR, addon = ""):
    extends_id = "800"
    mod = ps.Model(name, pathToSirane="SIRANE/")
    #Change results file
    mod.inputs.FICH_DIR_RESUL = f"RESULT_Analyse-sensibilite{addon}/{name}"
    #change dates
    mod.inputs.DATE_DEB = f"01/08/2018 00:00:00"
    mod.inputs.DATE_FIN = f"31/12/2018 23:00:00"
    #Change NO to NO2_equivalent
    mod.inputs.B_EMIS_NOEQNO2 = 1
    #Change dispertions parameter to default
    mod.site_disp.LATITUDE = 50.833
    mod.site_disp.Z0D = Z0D_disp#1.4
    mod.site_disp.ZDISPL = ZDISPL_disp#10
    mod.site_disp.ALBEDO = ALBEDO#0.15
    mod.site_disp.EMISSIVITE = EMISSIVITE#0.94
    mod.site_disp.PRIESTLEY_TAYLOR = PRIESTLEY_TAYLOR#0.75
    #Change meteo parameter to default
    mod.site_meteo.LATITUDE = 50.797
    mod.site_meteo.ALTITUDE = 27
    mod.site_meteo.Z0D = 1.5
    mod.site_meteo.ZDISPL = 15
    mod.site_meteo.ALBEDO = ALBEDO#0.15
    mod.site_meteo.EMISSIVITE = EMISSIVITE#0.94
    mod.site_meteo.PRIESTLEY_TAYLOR = PRIESTLEY_TAYLOR#0.75
    #Ne modeliser que les NOx et O3
    Esp = ps.read_dat(mod.inputs.FICH_ESPECES, mod_input = mod)
    Esp.loc[["PM",'PM25'],["Flag"]] = 0
    mod.save_dat(Esp, mod.inputs.FICH_ESPECES)
    #ajout d'un recepteur au dessus de celui d'Ixelles (1m au dessus de la rue canyon)
    recept = ps.read_dat(mod.inputs.FICH_RECEPT, mod_input=mod)
    recept.loc["IXL-Canope",:] = {"name":"Ixelles Rooftop", "X": 151115,  "Y": 168294,  "Z": 1,   "Type": 1, "Fichier":"RECEPTEURS/points_vides.dat"}
    mod.save_dat(recept, mod.inputs.FICH_RECEPT)
    #changer l'extend pour celui d'IXL avec un buffer de Xm
    mod.save_dat(extends["info-grid-meteo"].T[[extends_id]].T, mod.inputs.FICH_GRD_MET, index=False)
    mod.save_dat(extends["info-grid-sortie"].T[[extends_id]].T, mod.inputs.FICH_GRD_SORTIE, index=False)
    #desactive les stats:
    mod.inputs.B_CALC_STAT = 0
    mod.inputs.B_CALC_CRITERES = 0
    return mod

#Default parameters
Dparams = {"Z0D_disp" : 1.4, "ZDISPL_disp": 10, "ALBEDO": 0.15, "EMISSIVITE": 0.94, "PRIESTLEY_TAYLOR": 0.75}
#min and max typical values
MINparams = {"Z0D_disp" : 0.8, "ZDISPL_disp": 7, "ALBEDO": 0.1, "EMISSIVITE": 0.92, "PRIESTLEY_TAYLOR": 0.5}
MAXparams = {"Z0D_disp" : 1.5, "ZDISPL_disp": 15, "ALBEDO": 0.27, "EMISSIVITE": 0.961, "PRIESTLEY_TAYLOR": 1}
extend_Param = {"Min" : MINparams, "Def" : Dparams,"Max" : MAXparams}

#Creation des modèles
mod_list = []
mod_list.append(models("Def", **Dparams))
for MinMax in extend_Param:
    if MinMax != "Def":
        for param in Dparams:
            print(MinMax, ":", param)
            NEWparams = Dparams.copy()
            NEWparams[param] = extend_Param[MinMax][param]
            mod_list.append(models(f"{param}-{MinMax}", **NEWparams))

modandparams = pd.DataFrame({i.name : [i.site_disp.Z0D, i.site_disp.ZDISPL, i.site_disp.ALBEDO, i.site_disp.EMISSIVITE, i.site_disp.PRIESTLEY_TAYLOR]   for i in mod_list},
                            index=Dparams.keys()).T

#lancer les modèles
#ps.force_all_files_to_def(mod_list[-1]) #pour si un modèle s'est arreté précédement et a laisser des fichier "-def".
import datetime
for mod in mod_list:
    print(mod.name, "Started at", datetime.datetime.now())
    if os.path.isfile(f"{mod.pathToSirane}/{mod.inputs.FICH_DIR_RESUL}/log.txt"):
        None
    else:
        mod.run(silent=True)

#Récuperation des données et calcules des profiles
Fond = ps.read_dat("FOND/Concentration_fond.dat", mod_input=mod)#mod_list["41R002"])
filters = Fond[["NO","NO2",'O3']].isna().sum(axis=1) > 0
o = {}
for mod in mod_list:
    print(mod.name)
    mod2 = models(mod.name, **NEWparams,addon="-1erpartie")
    o[mod.name] = pd.concat([ps.Output(mod2).alldata.loc[~filters,["Recept","NO2_Mod"]].reset_index().set_index(["Recept","Date"]),ps.Output(mod).alldata.loc[~filters,["Recept","NO2_Mod"]].reset_index().set_index(["Recept","Date"])])

data = pd.concat(o, axis=1)
data = data.swaplevel(0,1, axis=1).NO2_Mod.reset_index().set_index("Date")
data["Dir"] = ps.Output(mod).alldata.loc[lambda x : x["Recept"] =="41R002","Dir"]
data["NO2_Mes"] = ps.Output(mod).alldata.loc[lambda x : x["Recept"] =="41R002","NO2_Mes"]
data["Dir"] = data.Dir.round(-1)/180 * np.pi#data.index.dayofweek + data.index.hour/24
meanData = data.groupby(["Recept","Dir"]).mean().reset_index()
meanData = meanData[meanData.Recept=="41R002"]
fig, axs= plt.subplots(2,3, layout="constrained", subplot_kw={'projection': 'polar', "theta_offset":np.pi/2, "theta_direction":-1})
axs = axs.flatten()
for i,element in enumerate(MINparams):
    print(element,i)
    axs[i].plot(meanData["Dir"], meanData[f"{element}-Min"], label="_Min")
    axs[i].plot(meanData["Dir"], meanData["Def"], label="_Def")
    axs[i].plot(meanData["Dir"], meanData[f"{element}-Max"], label="_Max")
    axs[i].plot(meanData["Dir"], meanData[f"NO2_Mes"], label="_Mes")

    axs[i].set(title=element, ylim=(20,120))

i+=1
axs[i].plot([0,0],[1,1], label="Min")
axs[i].plot([0,0],[1,1], label="Def")
axs[i].plot([0,0],[1,1], label="Max")
fig.legend(loc="lower right")
fig.delaxes(axs[i])
plt.savefig(resultfiles+"plotenFonctionVent.png")
plt.close("all")

#recuperations des données et calcules des quantiles
NO2_quant =[]
NO2_all = []
Fond = ps.read_dat("FOND/Concentration_fond.dat", mod_input=mod)#mod_list["41R002"])
filters = Fond[["NO","NO2",'O3']].isna().sum(axis=1) > 0
for mod in mod_list:
    print(mod.name)
    mod2 = models(mod.name, **NEWparams,addon="-1erpartie")
    o = pd.concat([ps.Output(mod2).alldata.loc[~filters,["Recept","NO2_Mod"]],ps.Output(mod).alldata.loc[~filters,["Recept","NO2_Mod"]]])
    o = o[o.Recept.apply(lambda x: x in ["41R002","IXL-Canope"])]
    o["param"] = mod.name.split("-")[0]
    o["MinMax"] = mod.name.split("-")[-1]
    o["mod_name"] = mod.name
    NO2_all.append(o.reset_index())
    q = o.loc[:,["Recept","NO2_Mod"]].groupby("Recept").quantile([0.05,0.5,0.95]).NO2_Mod.reset_index()
    q["param"] = mod.name.split("-")[0]
    q["MinMax"] = mod.name.split("-")[-1]
    q["mod_name"] = mod.name
    NO2_quant.append(q)

allNO2 = pd.concat(NO2_all)
allquants = pd.concat(NO2_quant)


allNO2.loc[allNO2.mod_name.isna(), "mod_name"] = "Def"
for param in modandparams:
    allNO2[param] = allNO2.mod_name.map(modandparams[param])
    allquants[param] = allquants.mod_name.map(modandparams[param])

#quantplot
import matplotlib as mpl
n = allquants.copy()
default = n[n.param =="Def"].loc[:,['Recept',"level_1", "NO2_Mod"]]
default.index = default.Recept + default.level_1.astype(str)
n.NO2_Mod = (n.NO2_Mod/(n.Recept+n.level_1.astype(str)).map(default.NO2_Mod)-1)*100
ParamName = {'Z0D_disp': "Rugosité aérodynamique [m]", 'ZDISPL_disp': "Epaisseur de déplacement [m]", 'ALBEDO': "Albédo [-]", 'EMISSIVITE': "Emissivité [-]", 'PRIESTLEY_TAYLOR': "Coefficient de Priesley-Talor [-]"}
n["Recept"] = n["Recept"].map({"IXL-Canope":"IXL-Canopée","41R002":"41R002"})
mpl.rcParams['font.size'] = 15
fig, axs = plt.subplots(2,3,sharey=True, layout="constrained")
axs = axs.flatten()
for i,param in enumerate(ParamName.keys()):
    sns.lineplot(n[(n.param == param)|(n.param == "Def")].sort_values(param).rename(columns={"Recept":"Récepteur", "level_1":"Quantile"}), x=param, y="NO2_Mod", style="Quantile", hue="Récepteur", sort=False, ax=axs[i],linewidth = 2.5)
    axs[i].set(xlabel = ParamName[param], ylabel="Différence relative [%]")

fig.delaxes(axs[-1])
plt.savefig(resultfiles+"difference relatve.png")
plt.close("all")
mpl.rcParams['font.size'] = 10

#Box plot
fig, axs_2D = plt.subplots(2,3, sharey=True)
axs= axs_2D.flatten()
for i,param in enumerate([i for i in set(allNO2.param) if i != "Def"]):
    sns.boxplot(data = allNO2[(allNO2.param == param)|(allNO2.param == "Def")],x=param, y="NO2_Mod", hue="Recept", ax=axs[i])

#pd.DataFrame(p)
plt.savefig(resultfiles+"Box plots.png")
plt.close("all")
