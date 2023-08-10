"""
conda activate PyGIS
python
#"""

mainfile = "PythonScripts/Validation/"

import PYSIRANE as ps
import pandas as pd
import datetime
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload as r
import seaborn as sns
import os
r(ps)
from tkinter import Tk
import time


resultfiles = "Resultats/"+mainfile.split("/")[-2]+"/"
os.makedirs(resultfiles,exist_ok=True)

extends = pd.read_csv("All_Extends.csv", header=[0,1], index_col=[0,1])
recept_order=np.array(["41R012", "41B011", "41MEU1","41B006","41R001","41B004","41WOL1", "41N043", "41R002", "41B008", "41B001"])

def clipboard(text):
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(text)
    r.update()
    time.sleep(1)
    r.update()
    r.destroy()

def models(name, extends_id, dateDeb = "01/01", dateFin="31/12"):
    mod = ps.Model(name, pathToSirane="SIRANE")
    mod.inputs.FICH_RESEAU = "RESEAU/Reseau_rues-SIRANE"
    #Change results file
    mod.inputs.FICH_DIR_RESUL = f"RESULT_Validation/{name}"
    #NO expressed in NO
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
    #desactive les stats et l'output sur la grille:
    mod.inputs.B_CALC_STAT = 0
    mod.inputs.B_CALC_CRITERES = 0
    mod.inputs.CALC_GRID = 0
    return mod

#Crée et lancer touts les modèles
#ps.force_all_files_to_def(models("..", extends.index[0]))
extends.iloc[[3],:]
mod_list = {}
for extends_id in extends.index:
    print(extends_id[0], ": started at", datetime.datetime.now())
    mod_list[extends_id[1]] = models(extends_id[0].replace(' ','_'),extends_id)
    if os.path.isfile(f"{mod_list[extends_id[1]].pathToSirane}/{mod_list[extends_id[1]].inputs.FICH_DIR_RESUL}/log.txt"):
        None
    else:
        mod_list[extends_id[1]].run(silent=True)


Receptors = ps.read_dat(mod_list[list(mod_list.keys())[0]].inputs.FICH_RECEPT,mod_input=mod_list[list(mod_list.keys())[0]])

#######
# Analyse
#######

#Creation des class output

out_list = {}
Fond = ps.read_dat("FOND/Concentration_fond.dat", mod_input=mod_list["41R002"])
filters = Fond[["NO","NO2",'O3']].isna().sum(axis=1) > 0
filters["2018":"2018"].mean()*100
r(ps)
for mod_id in mod_list:
    out_list[mod_id] = ps.Output(mod_list[mod_id], AddHourOfDay = True, Relative_Residual=False, extended_meteo=False)
    #Enelevr les données recepteurs qui ne sont pas ceux prévu de la zone
    list_recept = mod_id.split("/")
    filter_recept = out_list[mod_id].alldata.Recept.apply(lambda x : x in list_recept)
    out_list[mod_id].alldata = out_list[mod_id].alldata[filter_recept]
    out_list[mod_id].list_of_Recept = list_recept
    for esp in ["NO","NO2",'O3']: #Enlever les données modélisées si l'un des trois est absent sur la station de fond
        out_list[mod_id].alldata.loc[filters,[f"{esp}_Mod",f"{esp}_Res"]] = np.nan
    out_list[mod_id].alldata.loc[(out_list[mod_id].alldata.Recept=="41R002")&(out_list[mod_id].alldata.index =="2018-08-01 00:00:00"), ["NO2_Mes","NO_Mes"]] = np.nan

#cree un Output combiné de touts les sous-moddèles.
r(ps)
out_list["All"] = ps.Output(models(extends_id[0].replace(' ','_'),extends_id), AddHourOfDay = True, Relative_Residual=False, extended_meteo=False)
out_list["All"].alldata = pd.concat([out_list[id].alldata for id in mod_list])
#out_list["All"].alldata.Recept = "All"
out_list["All"].list_of_Recept = [i for id8 in mod_list for i in out_list[id8].list_of_Recept]#["All"]#
out_list["All"].list_of_species = ["NO","NO2",'O3']
out_list["All"].mod.name="All"

#Calcule des indicateurs statistiques
inds = pd.concat(out_list["All"].indicateurs(yearlyMQI=True))
inds["order"] = [(np.array(["NO","NO2","O3"])==esp).argmax()*100 + (recept_order==recp).argmax() for esp, recp in inds.index]
inds = inds.sort_values("order")

inds[["FB","NMSE","R","FAC2","NAD","MQI","MQI_y"]].to_csv(resultfiles+"indicateurs_stats.csv")


def analise_indicateurs(all_ind):
    Rural = all_ind[["FB", "NMSE", "FAC2", "NAD"]].dropna(axis=0).copy()
    Rural["FB"] = all_ind.FB.abs() <= 0.3
    #Rural["MG"] = (all_ind.MG <= 1.3)|(all_ind.MG >= 0.7)
    Rural["NMSE"] = all_ind.NMSE <= 3
    #Rural["VG"] = all_ind.VG <= 11
    Rural["FAC2"] = all_ind.FAC2 >= 0.5
    Rural["NAD"] = all_ind.NAD <=0.3
    Rural["Respect"] = Rural.sum(axis=1) >= 2
    Rural_O = Rural.reset_index().groupby("level_0")["Respect"].mean() >= 0.5
    Urbain = all_ind[["FB", "NMSE", "FAC2", "NAD"]].dropna(axis=0).copy()
    Urbain["FB"] = all_ind.FB.abs() <= 0.67
    Urbain["NMSE"] = all_ind.NMSE <= 6
    Urbain["FAC2"] = all_ind.FAC2 >= 0.3
    Urbain["NAD"] = all_ind.NAD <=0.5
    Urbain["Respect"] = Urbain.sum(axis=1) >= 2
    Urbain_O = Urbain.reset_index().groupby("level_0")["Respect"].mean() >= 0.5
    def MQO_test(Serie):
        Serie = Serie[Serie>=0]
        MQI = Serie.dropna().sort_values().values
        if len(MQI) == 0:
            return np.nan
        #print(MQI, Serie)
        S90 = int(len(MQI)*0.9)-1
        dist = len(MQI)*0.9 - S90 - 1
        MQI90th = MQI[S90] + (MQI[S90+1] - MQI[S90])*dist
        MQO = MQI90th
        return MQO
    MQO = all_ind.reset_index().groupby("level_0")["MQI"].apply(MQO_test)
    MQO_y = all_ind.reset_index().groupby("level_0")["MQI_y"].apply(MQO_test)
    return {"table" : pd.concat({"Rural":Rural, "Urbain":Urbain}, axis=1), "Rural": Rural_O, "Urbain":Urbain_O,"MQO":MQO, "MQO_y":MQO_y}

print(analise_indicateurs(inds))


#ajout des NOx et NO2_relatif
out_list["All"].alldata["NOx_Mes"] = out_list["All"].alldata["NO2_Mes"] + out_list["All"].alldata["NO_Mes"] * 46/30
out_list["All"].alldata["NOx_Mod"] = out_list["All"].alldata["NO2_Mod"] + out_list["All"].alldata["NO_Mod"] * 46/30
out_list["All"].alldata["NO2R_Mes"] = out_list["All"].alldata["NO_Mes"]/  out_list["All"].alldata["NOx_Mes"]
out_list["All"].alldata["NO2R_Mod"] = out_list["All"].alldata["NO_Mod"]/  out_list["All"].alldata["NOx_Mod"]

roll = out_list["All"].alldata.groupby("Recept").apply(lambda x : x.reset_index().set_index(["Recept","Date"]).rolling(24*30*3, min_periods = 24*20*3).mean().reset_index().set_index("Date"))
RollNO2 = roll.reset_index().pivot_table(values=["NO_Mod","NO_Mes"],index="Date",columns="Recept").swaplevel(1,0,axis=1)
f,axs = plt.subplots(3,4,layout="constrained")
axs=axs.flatten()
Receptors["FullName"] = Receptors.name + " (" + Receptors.index +")"
for i,recept in enumerate(recept_order):
    (RollNO2[recept]/RollNO2[recept].dropna().mean()).plot(ax=axs[i])
    axs[i].set_title(Receptors["FullName"][recept])

plt.savefig(resultfiles+"rollingavrages_NO2overmeanNO2.png")
plt.close("all")


out_list[mod_id].alldata["Recept_name"] = out_list[mod_id].alldata.Recept.map(extends.reset_index().set_index("Id")["Name"])

# Graphiques pour toutes les stations de mesures
plt.rcParams['figure.constrained_layout.use'] = True
r(ps)
for mod_id in mod_list:
    out_list[mod_id] = ps.Output(mod_list[mod_id], AddHourOfDay = True, Relative_Residual=False, extended_meteo=False)
    #Enelevr les données recepteurs qui ne sont pas ceux prévu de la zone
    list_recept = mod_id.split("/")
    filter_recept = out_list[mod_id].alldata.Recept.apply(lambda x : x in list_recept)
    out_list[mod_id].alldata = out_list[mod_id].alldata[filter_recept]
    out_list[mod_id].list_of_Recept = list_recept
    for esp in ["NO","NO2",'O3']: #Enlever les données modélisées si l'un des trois est absent sur la station de fond
        out_list[mod_id].alldata.loc[filters,[f"{esp}_Mod",f"{esp}_Res"]] = np.nan
    out_list[mod_id].alldata.loc[(out_list[mod_id].alldata.Recept=="41R002")&(out_list[mod_id].alldata.index =="2018-08-01 00:00:00"), ["NO2_Mes","NO_Mes"]] = np.nan
    out_list[mod_id].alldata["NOx_Mes"] = out_list[mod_id].alldata["NO2_Mes"] + out_list[mod_id].alldata["NO_Mes"] * 46/30
    out_list[mod_id].alldata["NOx_Mod"] = out_list[mod_id].alldata["NO2_Mod"] + out_list[mod_id].alldata["NO_Mod"] * 46/30
    out_list[mod_id].alldata["NOx_Res"] = out_list[mod_id].alldata["NOx_Mod"]-out_list[mod_id].alldata["NOx_Mes"]
    out_list[mod_id].list_of_species = np.append(out_list[mod_id].list_of_species,["NOx"])
    for esp in ["NO2",'O3']:
        out_list[mod_id].alldata.loc[:,f"{esp}_Res"] = out_list[mod_id].alldata.loc[:,f"{esp}_Res"]/U(out_list[mod_id].alldata[f"{esp}_Mes"],esp=esp)


for mod_id in out_list:#["All"]:#out_list:#
    if mod_id=="All":
        break
    for Receptor in mod_id.split("/"):
        plt.close("all")
        fig, axs = out_list[mod_id].scatterplots(list_of_Recept=[Receptor])
        if fig is not None:
            fig.suptitle(f"{Receptors['name'][Receptor]} ({Receptor})",fontsize=15)
            fig.set_size_inches(5* len(axs),5)
            os.makedirs(resultfiles+f"/scatterplots", exist_ok=True)
            fig.savefig(resultfiles+f"/scatterplots/{Receptor}.png")
            fig.savefig(resultfiles+f"/scatterplots/{Receptor}.svg")
        plt.close("all")
        fig, axs = out_list[mod_id].qQplots(list_of_Recept=[Receptor])
        if fig is not None:
            fig.suptitle(f"{Receptors['name'][Receptor]} ({Receptor})",fontsize=15)
            fig.set_size_inches(5* axs.shape[1],5*axs.shape[0])
            os.makedirs(resultfiles+f"/qQplots", exist_ok=True)
            fig.savefig(resultfiles+f"/qQplots/{Receptor}.png")
            fig.savefig(resultfiles+f"/qQplots/{Receptor}.svg")
        plt.close("all")
        fig, axs = out_list[mod_id].residualplots(list_of_Recept=[Receptor])
        if fig is not None:
            fig.suptitle(f"{Receptors['name'][Receptor]} ({Receptor})")
            fig.set_size_inches(axs.shape[1] * 7,30)
            #[ax.set_yscale("log") for ax in axs.flat]
            os.makedirs(resultfiles+f"/residualplots", exist_ok=True)
            fig.savefig(resultfiles+f"/residualplots/{Receptor}-abs.png")
            fig.savefig(resultfiles+f"/residualplots/{Receptor}-abs.svg")
        plt.close("all")

####Plot de serie temporelles


out_list["All"].alldata = pd.concat([out_list[id].alldata for id in mod_list if id != "All"])

fig, axs = plt.subplots(len(set(out_list["All"].alldata.Recept)),1, sharex=True)
for i,recept in enumerate(set(out_list["All"].alldata.Recept)):
    out_list["All"].alldata.loc[out_list["All"].alldata.Recept==recept,["NO2_Mod","NO2_Mes"]].plot(ax=axs[i])
    axs[i].set(ylabel=recept)
plt.show(block=False)

for id in ["All"]:#out_list
    for Receptor in id.split("/"):
        plt.close("all")
        print(id, Receptor)
        fig, axs = out_list[id].residualplots(list_of_VarMeteo=[f"{s}_Mes" for s in out_list[id].list_of_species],
                                   list_of_Recept=[Receptor])
        if fig is not None:
            fig.suptitle(f"{Receptors['name'][Receptor]} ({Receptor})")
            fig.savefig(resultfiles+f"/residualplots/{Receptor}-Rel_over_Mes.png")
            fig.savefig(resultfiles+f"/residualplots/{Receptor}-Rel_over_Mes.svg")

plt.close("all")

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
out_list["All"].alldata["Type"] = out_list["All"].alldata.Recept.map(Type_Recept)

#Profils moyen en fonction de la direction du vent
import pytz
SIRANE_tz = pytz.FixedOffset(+1*60) #GMT+1 : le fuseau horaire de sirane
Data_tz = "CET"


out_list["All"].alldata["Date_Local"] = out_list["All"].alldata.index.tz_localize(SIRANE_tz)#.tz_convert(Data_tz)
out_list["All"].alldata["Heure"] = out_list["All"].alldata["Date_Local"].dt.hour

out_list["All"].alldata["typeJour"] = (out_list["All"].alldata.Jour_semaine-4).clip(0,10)
out_list["All"].alldata["Jours"] = out_list["All"].alldata.Heure/24 + out_list["All"].alldata.Jour_semaine
out_list["All"].alldata["Dir_Group"] = pd.cut(out_list["All"].alldata.Dir/180*np.pi, 35)

data = out_list["All"].alldata[out_list["All"].alldata.isna()[["NO2_Mod","NO2_Mes"]].sum(axis=1)==0].groupby(["Dir_Group","Day","Recept"])[["NO2_Mod","NO2_Mes"]].mean().rename(columns={"NO2_Mod":"Modèle", "NO2_Mes":"Mesure"}).unstack()
data = data.swaplevel(0,1, axis=1)
data = data.reset_index().set_index("Dir_Group")
data.index = pd.IntervalIndex(data.index).mid
#data = data.reset_index().set_index(["index","Day"])

Receptors["FullName"] = Receptors.name + " (" + Receptors.index +")"
from matplotlib.offsetbox import AnchoredText
f, axs = plt.subplots(3,4, figsize=(13,10), layout="constrained")
axs = axs.flatten()
miny,maxy = data.min().min()*0.5, data.max().max()*1.1
for i,recept in enumerate(recept_order):
    axs[i].plot(data[recept].index, data[recept]["Modèle"], label="_Modèle")
    axs[i].plot(data[recept].index, data[recept]["Mesure"], label="_Mesure")
    axs[i].set_title(Receptors["FullName"][recept])
    axs[i].set(xlabel="", ylabel="")
    axs[i].tick_params(which="both", bottom=True,left=False,labelbottom=True, labelleft=True)

i+=1
axs[i].plot(data[recept].reset_index()["index"], data[recept]["Modèle"], label="Modèle")
axs[i].plot(data[recept].reset_index()["index"], data[recept]["Mesure"], label="Mesure")
f.legend(loc="lower right")
f.delaxes(axs[-1])
plt.savefig(resultfiles+"profils moyen du vent.png")
plt.close("all")

#profils hebdomadaires
data = out_list["All"].alldata[out_list["All"].alldata.isna()[["NO2_Mod","NO2_Mes"]].sum(axis=1)==0].groupby(["Jours","Day","Recept"])[["NO2_Mod","NO2_Mes"]].mean().rename(columns={"NO2_Mod":"Modèle", "NO2_Mes":"Mesure"}).unstack()
data = data.swaplevel(0,1, axis=1)
data = data.reset_index().set_index("Jours")

f, axs = plt.subplots(3,4, figsize=(13,10), layout="constrained")
axs = axs.flatten()
miny,maxy = data.min().min()*0.5, data.max().max()*1.1
for i,recept in enumerate(recept_order):
    axs[i].plot(data[recept].index, data[recept]["Modèle"], label="_Modèle")
    axs[i].plot(data[recept].index, data[recept]["Mesure"], label="_Mesure")
    axs[i].set_title(Receptors["FullName"][recept])
    axs[i].set(xlabel="", ylabel="")
    axs[i].tick_params(which="both", bottom=True,left=False,labelbottom=True, labelleft=True)

i+=1
axs[i].plot(data[recept].reset_index()["index"], data[recept]["Modèle"], label="Modèle")
axs[i].plot(data[recept].reset_index()["index"], data[recept]["Mesure"], label="Mesure")
f.legend(loc="lower right")
f.delaxes(axs[-1])
plt.savefig(resultfiles+"profils hebdomadaires.png")
plt.close("all")


#Histogrames
def plothists(data, esp, ax=None,hide="_"):
    x = data.loc[~data[[f"{esp}_Mes",f"{esp}_Mod"]].isna().any(axis=1),f"{esp}_Mes"]
    y = data.loc[~data[[f"{esp}_Mes",f"{esp}_Mod"]].isna().any(axis=1),f"{esp}_Mod"]
    maxval = max(np.max(x),np.max(y))
    minval = min(np.min(x),np.min(y))
    bins = np.linspace(minval,maxval,50)
    if ax is None:
      ax = plt.gca()
    ax.hist(x,histtype="step",bins=bins, label=hide+"Mesure", density=True)
    ax.hist(y,histtype="step",bins=bins, label=hide+"Modèle", density=True)
    return ax

fig, axs = plt.subplots(4,3,layout="constrained")
axs=axs.flatten()
esp = "NO2"
for mod_id in out_list:#["All"]:#out_list:#
    if mod_id=="All":
        break
    for Receptor in mod_id.split("/"):
        i = np.argmax(recept_order==Receptor)
        plothists(out_list[mod_id].alldata[out_list[mod_id].alldata.Recept == Receptor], esp, ax=axs[i])
        axs[i].set_title(f"{Receptors['name'][Receptor]} ({Receptor})")
        #axs[i].set_yscale("log")
        axs[i].set(xlim=(1,out_list[mod_id].alldata.loc[out_list[mod_id].alldata.Recept == Receptor,[f"{esp}_Mes",f"{esp}_Mod"]].max().max()))
        #axs[i].set_xscale("log")

plothists(out_list[mod_id].alldata[out_list[mod_id].alldata.Recept == Receptor], esp, ax=axs[-1],hide="")
fig.legend(loc="lower right")
fig.delaxes(axs[-1])
plt.savefig(resultfiles+"HistogramesDeNO2.png")
plt.close("all")
