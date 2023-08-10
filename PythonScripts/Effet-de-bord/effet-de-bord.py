"""
Le script execute SIRANE pour avec différentes taille de modélisaiton de 50 à 1.6km pour deux semaine de 2018

#"""

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


##########
# graphique de l'effet de la réduction de la zone
##########
[i.name for i in mod_list]
NO2_quant =[]
NO2_all = []
for mod in mod_list:
    o = ps.Output(mod).alldata.loc[:,["Recept","NO2_Mod","NO_Mod"]]
    o = o[o.Recept.apply(lambda x: x in ["41R002","IXL-Canope"])]
    o["semaine"] = int(mod.name.split("-")[0][1])
    o["buffer"] = mod.name.split("-")[1][:-1]
    NO2_all.append(o.reset_index())
    q = o.loc[:,["Recept","NO_Mod"]].groupby("Recept").quantile([0.05,0.5,0.95]).NO_Mod.reset_index()
    q["semaine"] = int(mod.name.split("-")[0][1])
    q["buffer"] = mod.name.split("-")[1][:-1]
    NO2_quant.append(q)

allNO2 = pd.concat(NO2_all)

allquants = pd.concat(NO2_quant)
allquants.columns = ["Recepteur", "Quantile", "C_NO2", "Semaine", "Taille zone tampon"]
ax2 = sns.lineplot(data=allquants.loc[(allquants.Semaine == 1)], x="Taille zone tampon", y="C_NO2", hue="Recepteur",style="Quantile")
plt.savefig(resultfiles+"EvolQuantilesNO2.png")
plt.close("all")

allquants_pivot = allquants.pivot(values = "C_NO2", columns = "Taille zone tampon", index=["Recepteur", "Quantile", "Semaine"])
for col in allquants_pivot:
    allquants_pivot[col] = allquants_pivot[col] / allquants_pivot["BX"]

allquants_rel = allquants_pivot[allquants_pivot.columns[:-1]].stack()._set_name("C_NO2").reset_index()
allquants_rel.loc[:,["Quantile", "C_NO2", "Semaine", "Taille zone tampon"]] = allquants_rel[["Quantile", "C_NO2", "Semaine", "Taille zone tampon"]].astype(float)
allquants_rel["Zones"] = "Zone ouverte"
allquants_rel.loc[allquants_rel.Recepteur == "41R002",["Zones"]] = "Rue canyon"
fig, axs = plt.subplots(1,2,sharey=True)
sns.lineplot(ax = axs[0],data=allquants_rel.loc[(allquants_rel.Semaine == 1)], x="Taille zone tampon", y="C_NO2", hue="Zones",style="Quantile")
sns.lineplot(ax = axs[1],data=allquants_rel.loc[(allquants_rel.Semaine == 2)], x="Taille zone tampon", y="C_NO2", hue="Zones",style="Quantile")
axs[0].set(ylim=(0,1), ylabel="Concentration de NO2 modélisée relative \n à la modélisation sur tout Bruxelles [-]", xlabel="Taille de la zone tampon [m]")
axs[1].set(xlabel="Taille de la zone tampon [m]")
axs[0].set_title("Du 13/05 au 19/05")
axs[1].set_title("Du 23/05 au 29/05")
plt.savefig(resultfiles+"EvolQuantilesNO2_Relatif.png")
plt.close("all")

allNO2 = allNO2.reset_index(drop=True)
allNO2["Name"] = allNO2.Date.astype(str)+allNO2.Recept
allNO2["RatioNO2"] = allNO2.NO2_Mod/allNO2.Name.map(allNO2.loc[allNO2.buffer == "BX", ["NO2_Mod","Name"]].set_index("Name").NO2_Mod)
fig, axs = plt.subplots(1,2, layout="constrained")
sns.boxplot(data=allNO2[(allNO2.semaine==1)&(allNO2.buffer != "BX")], x="buffer", y="RatioNO2", hue="Recept", ax=axs[0])
sns.boxplot(data=allNO2[(allNO2.semaine==2)&(allNO2.buffer != "BX")], x="buffer", y="RatioNO2", hue="Recept", ax=axs[1])
axs[0].set_title("Du 13/05 au 19/05")
axs[1].set_title("Du 23/05 au 29/05")
axs[0].set(xlabel="", ylabel="", ylim=(0,1))
axs[1].set(xlabel="", ylabel="", ylim=(0,1))
plt.savefig(resultfiles+"BoxplotNO2_Relative.png")
plt.close("all")

fig, axs = plt.subplots(2,1, layout="constrained")
sns.boxplot(data=allNO2[allNO2.semaine==1], x="buffer", y="NO2_Mod", hue="Recept", ax=axs[0])
sns.boxplot(data=allNO2[allNO2.semaine==2], x="buffer", y="NO2_Mod", hue="Recept", ax=axs[1])
axs[0].set_title("Du 13/05 au 19/05")
axs[1].set_title("Du 23/05 au 29/05")
axs[0].set(xlabel="", ylabel="")
axs[1].set(xlabel="", ylabel="")
plt.savefig(resultfiles+"BoxplotNO2.png")
plt.close("all")

############
# Creation des extends
############

#show bx and receptors
import geopandas as gpd
mod = out_list["S1-BXm"].mod
sources_ponct = ps.read_dat("EMISSIONS/EMIS_PONCT/sources-ponct.txt",mod_input=mod)
Receptors = ps.read_dat(mod.inputs.FICH_RECEPT,mod_input=mod)
Communes_RBC = gpd.read_file("shp_file/Communes_RBC.shp")
import rasterio
from rasterio.mask import mask
import rasterio.plot as rplot

#creation des bbox
import geopandas as gpd
from shapely.geometry import MultiPoint
def get_combinations(lst):
    lst = tuple(lst)
    output = set()
    n = len(lst)
    def combine(start, sublist):
        if start == n and sublist:
            output.add(tuple(sublist)) # convertir chaque sous-liste en tuple pour pouvoir l'ajouter à l'ensemble
        else:
            for i in range(start, n):
                combine(i+1, sublist+[lst[start:i+1]])
    combine(0, [])
    return list(output)

def find_Comb(zones):
    zone_combinations = get_combinations(zones.index)
    rectangles = []
    for possibleCom in zone_combinations:
        poss_rect = []
        for rect in possibleCom:
            union = gpd.GeoSeries(list(map(lambda x: zones.geometry[x],rect))).unary_union
            poss_rect.append(union.bounds)
        rectangles.append([MultiPoint((rec[:2],rec[2:])).envelope for rec in poss_rect])
    a = pd.concat({id: gpd.GeoDataFrame(geometry=rect).set_crs(zones.crs) for id,rect in enumerate(rectangles)})
    return a.loc[np.argmin(a.area.unstack().sum(axis=1)),:]

def create_rect(Receptors_gdf):
    #print(Receptors_gdf)
    union = gpd.GeoDataFrame(geometry=[Receptors_gdf.unary_union]).explode(index_parts=False).reset_index(drop=True).set_crs(Receptors_gdf.crs)
    #union["Grouped_Geom"] = union.geometry.copy()
    Receptors_Grouped = gpd.sjoin(union, Receptors_gdf, how='right')
    rectangles = np.array([])
    for groupID in union.index:
        rect_grouped = Receptors_Grouped[Receptors_Grouped["index_left"] == groupID].copy()
        if len(rect_grouped) == 1 :
            rectangles = np.append(rectangles,rect_grouped["geometry"][0].envelope)
            #rectangles[groupID] = rect_grouped["geometry"][0]
        else:
            #print(groupID,find_Comb(rect_grouped))
            rectangles = np.append(rectangles,find_Comb(rect_grouped))
    return rectangles

Receptors = ps.read_dat(mod.inputs.FICH_RECEPT,mod_input=mod)
Receptors = Receptors[(Receptors.index != "41REG1")&(Receptors.index != "41CHA1")]
Receptors["Bsize"] = 800
Receptors_gdf = gpd.GeoDataFrame(Receptors, geometry=gpd.points_from_xy(Receptors.X, Receptors.Y, crs = "EPSG:31370"))
Receptors_gdf = Receptors_gdf[Receptors_gdf.index != "IXL-Canope"]
Receptors_gdf.geometry = Receptors_gdf.buffer(Receptors_gdf.Bsize, cap_style=3)
rectangles = gpd.GeoSeries(create_rect(Receptors_gdf)).set_crs(Receptors_gdf.crs)
print(f"Reduction of {100 - rectangles.area.sum()/Receptors_gdf.envelope.area.sum()*100:.2f}%")

plt.close()
ax = Receptors_gdf.plot(edgecolor="k",facecolor=(0,0,1,0.3))
rectangles.plot(edgecolor="r", facecolor=(1,0,0,0.), ax=ax)
plt.savefig(resultfiles+"Extents graphs.png")
plt.close("all")




bboxs = rectangles.bounds
bboxs.index = pd.MultiIndex.from_frame(pd.DataFrame([["/".join(Receptors_gdf.name[Receptors_gdf.within(rect)]), "/".join(Receptors_gdf.index[Receptors_gdf.within(rect)])] for rect in rectangles]),
                                       names=['Name', 'Id'])
bboxs["Cell_size_meteo"] = 100
bboxs["Cell_size_sortie"] = 10
grids = {}
for type in ["meteo", "sortie"]:
    for x in ["x","y"]:
        bboxs[f"center{x}"] = round((bboxs[f"min{x}"] + bboxs[f"max{x}"])/2)
        delta = np.ceil((-bboxs[f"min{x}"] + bboxs[f"max{x}"])/bboxs[f"Cell_size_{type}"])*bboxs[f"Cell_size_{type}"]
        if type=="sortie":
            delta = np.ceil((-bboxs[f"min{x}"] + bboxs[f"max{x}"]-2*bboxs[f"Cell_size_{type}"])/bboxs[f"Cell_size_{type}"])*bboxs[f"Cell_size_{type}"]
        bboxs[f"{x}min"] =  bboxs[f"center{x}"] - delta/2
        bboxs[f"{x}max"] =  bboxs[f"center{x}"] + delta/2
        bboxs[f"N{x}"] = -(bboxs[f"{x}min"] - bboxs[f"{x}max"])/bboxs[f"Cell_size_{type}"]+1
        bboxs[[f"{x}min",f"{x}max",f"N{x}"]] = bboxs[[f"{x}min",f"{x}max",f"N{x}"]].astype(int)
    grids[f"info-grid-{type}"] = bboxs[['Nx','Ny','xmin','xmax','ymin','ymax']]

pd.concat(grids, axis=1).to_csv("PythonScipts/Validation/All_Extends.csv")
