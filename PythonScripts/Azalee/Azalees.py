"""
conda activate PyGIS2
python
#"""

mainfile = "PythonScripts/Azalee/"

import PYSIRANE as ps
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib_scalebar.scalebar import ScaleBar
import geopandas as gpd
import pandas as pd
import numpy as np
from scipy import stats
import glob,os
import rasterio as rio
#from importlib import reload as r
#r(ps)

resultfiles = "Resultats/"+mainfile.split("/")[-2]+"/"
os.makedirs(resultfiles,exist_ok=True)

extends = pd.read_csv(mainfile+"Azalee_Extend.csv", header=[0,1], index_col=[0])

def map_subplots(nrow,ncol, subplots_kws=None, Drawarrow=False, arrow_width=5, arrow_length=0.1):
    if subplots_kws is None:
        subplots_kws = {"layout":"constrained", "sharex":True, "sharey":True}
    fig,axs = plt.subplots(nrow,ncol, **subplots_kws)
    for ax in np.array(axs).flat:
        ax.tick_params(labelbottom=True, labelleft=True)
        ax.add_artist(ScaleBar(1))
        if Drawarrow:
            x, y, arrow_length, arrow_width = 0.9, 0.9, 0.1, 5
            ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                        arrowprops=dict(facecolor='black',edgecolor='w', width=arrow_width, headwidth=3*arrow_width),
                        ha='center', va='center', fontsize=4*arrow_width,
                        xycoords=ax.transAxes)
    return fig, axs

def model(name, dateDeb = "01/01", dateFin="31/12"):
    mod = ps.Model(name, pathToSirane="SIRANE")
    mod.inputs.FICH_RESEAU = "RESEAU/Reseau_rues-SIRANE"
    #Change results file
    mod.inputs.FICH_DIR_RESUL = f"RESULT_Azalees_big/{name}"
    #NO expressed in NO2eq
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
    #pas de calcule de stats et critères
    mod.inputs.B_CALC_STAT = 0
    mod.inputs.B_CALC_CRITERES = 0
    #change extend #[["Nx","Ny","xmin","xmax","ymin","ymax"]]
    mod.save_dat(extends["info-grid-meteo"].T[["Azalee_big"]].T[["Nx","Ny","xmin","xmax","ymin","ymax"]], mod.inputs.FICH_GRD_MET, index=False)
    mod.save_dat(extends["info-grid-sortie"].T[["Azalee_big"]].T[["Nx","Ny","xmin","xmax","ymin","ymax"]], mod.inputs.FICH_GRD_SORTIE, index=False)
    #ajout d'un recepteur au dessus de celui d'Ixelles (1m au dessus de la rue canyon)
    recept = ps.read_dat(mod.inputs.FICH_RECEPT, mod_input=mod)
    #recept.loc["IXL-Canope",:] = {"name":"Azalee", "X": 150703,  "Y": 171636,  "Z": 1,   "Type": 1, "Fichier":"RECEPTEURS/points_vides.dat"}
    mod.save_dat(recept, mod.inputs.FICH_RECEPT)
    return mod

#modélisation initiale
moddef = model("def")
#moddef.inputs.AFFICH=2
#ps.force_all_files_to_def(moddef)
#if not os.path.isfile(f"{moddef.pathToSirane}/{moddef.inputs.FICH_DIR_RESUL}/log.txt"):
#    moddef.run(silent=True)

#modélisation Après
mod_p = model("Apres")
data_Evol_Emis_Trafic = ps.read_dat("EMISSIONS/EMIS_LIN/Evol_Emis_Trafic.dat", mod_input=mod_p)
data_Evol_Emis_Trafic.Fich_Emis = "EMISSIONS/EMIS_LIN/Emis_rues-modified_Apres.dat"
mod_p.save_dat(data_Evol_Emis_Trafic, "EMISSIONS/EMIS_LIN/Evol_Emis_Trafic.dat")
mod_p.inputs.B_CALC_STAT = 1
stat = ps.read_dat(mod_p.inputs.FICH_STATS, mod_input=mod_p)
stat = pd.concat([stat,pd.DataFrame({"Espece":"NO2", "Calcul":"1", "Grille":1,"Recept":1,"Recept_Mes":1, "Rue":0,"Valeur":np.nan}, index=pd.Index(["Moy_Horaire"], name="Type"))])
mod_p.save_dat(stat, mod_p.inputs.FICH_STATS)
if not os.path.isfile(f"{mod_p.pathToSirane}/{mod_p.inputs.FICH_DIR_RESUL}/log.txt"):
    mod_p.run(silent=True)

#modélisation Avant
mod_v = model("Avant")
mod_v.inputs.B_CALC_STAT = 1
stat = ps.read_dat(mod_v.inputs.FICH_STATS, mod_input=mod_v)
stat = pd.concat([stat,pd.DataFrame({"Espece":"NO2", "Calcul":"1", "Grille":1,"Recept":1,"Recept_Mes":1, "Rue":0,"Valeur":np.nan}, index=pd.Index(["Moy_Horaire"], name="Type"))])
mod_v.save_dat(stat, mod_v.inputs.FICH_STATS)
data_Evol_Emis_Trafic = ps.read_dat("EMISSIONS/EMIS_LIN/Evol_Emis_Trafic.dat", mod_input=mod_v)
data_Evol_Emis_Trafic.Fich_Emis = "EMISSIONS/EMIS_LIN/Emis_rues-modified_Avant.dat"
mod_v.save_dat(data_Evol_Emis_Trafic, "EMISSIONS/EMIS_LIN/Evol_Emis_Trafic.dat")
if not os.path.isfile(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/log.txt"):
    mod_v.run(silent=True)



##### Enregistrement des concentrations avant et après
Fond = ps.read_dat("FOND/Concentration_fond.dat", mod_input=mod_p)#mod_list["41R002"])
filters = Fond[["NO","NO2",'O3']].isna().sum(axis=1) > 0
grids_path = glob.glob(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/GRILLE/Conc_NO2_*.grd")
grids_path = pd.DataFrame({p.split("_")[-1].split(".")[0]:{"Espece":p.split("_")[-2],"Fichier":p} for p in grids_path}).T
grids_path.index = pd.to_datetime(grids_path.index, format="%Y%m%d%H")
grids_path2 = glob.glob(f"{mod_p.pathToSirane}/{mod_p.inputs.FICH_DIR_RESUL}/GRILLE/Conc_NO2_*.grd")
grids_path2 = pd.DataFrame({p.split("_")[-1].split(".")[0]:{"Espece":p.split("_")[-2],"Fichier":p} for p in grids_path2}).T
grids_path2.index = pd.to_datetime(grids_path2.index, format="%Y%m%d%H")
grids_path = grids_path.loc[~filters,:] #bazarder les données pas bonnes
grids_path2 = grids_path2.loc[~filters,:] #bazarder les données pas bonnes
if not (os.path.isfile(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/GRILLE/All_NO2.npz") and os.path.isfile(f"{mod_p.pathToSirane}/{mod_p.inputs.FICH_DIR_RESUL}/GRILLE/All_NO2.npz")):
    l = len(grids_path)
    len(grids_path),len(grids_path2)
    ar, ext, t = ps.transph_grd(grids_path.Fichier.iloc[0])
    ultrabig_data = np.ndarray(shape = (ar.shape[0],ar.shape[1],l), dtype=np.float16)
    ultrabig_data[:,:,0] = ar[:,:]
    #
    ar, ext, t = ps.transph_grd(grids_path2.Fichier.iloc[0])
    ultrabig_data2 = np.ndarray(shape = (ar.shape[0],ar.shape[1],l), dtype=np.float16)
    ultrabig_data2[:,:,0] = ar[:,:]
    ar.shape
    for id, g_p in enumerate(grids_path.index[1:]):
        print(f"{id}/{l} ({id/l*100}%):",g_p)
        ar, ext, t = ps.transph_grd(grids_path.Fichier[g_p])
        ultrabig_data[:,:,id+1]  = ar[:,:]
        ar, ext, t = ps.transph_grd(grids_path2.Fichier[g_p])
        ultrabig_data2[:,:,id+1]  = ar[:,:]
        #
    np.savez(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/GRILLE/All_NO2.npz", data = ultrabig_data, extent = ext)
    np.savez(f"{mod_p.pathToSirane}/{mod_p.inputs.FICH_DIR_RESUL}/GRILLE/All_NO2.npz", data = ultrabig_data2,extent = ext)
else :
    with np.load(f'{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/GRILLE/All_NO2.npz') as data:
        ext = data['extent']
        ultrabig_data = data["data"]
    with np.load(f'{mod_p.pathToSirane}/{mod_p.inputs.FICH_DIR_RESUL}/GRILLE/All_NO2.npz') as data:
        ultrabig_data2 = data["data"]


data = pd.DataFrame({"Avant":ultrabig_data[y,x,:], "Apres" : ultrabig_data2[y,x,:]}, index=grids_path.index)
data["Heure"] = data.index.hour
data["Jour"] = data.index.dayofweek + data["Heure"]/24
data["Evol"] = (data.Apres-data.Avant) / data.Avant * 100
data.groupby(["Jour"]).mean()[["Avant","Apres","Evol"]].plot()
plt.show(block = False)

#Moyenne Annuelle
NO2_moy_after =ultrabig_data2.mean(axis=2)
NO2_moy_before = ultrabig_data.mean(axis=2)
diff = (NO2_moy_after-NO2_moy_before)
max_diff = np.abs(diff).max()
rel_dif = diff/NO2_moy_before * 100
max_rel_dif = np.abs(rel_dif).max()

with rio.open(grids_path.Fichier.iloc[0]) as src:
    kwargs = src.meta.copy()
    kwargs["crs"] = "EPSG:31370"
    kwargs["driver"] = "GTiff"
    os.makedirs(resultfiles+f"/maps", exist_ok=True)
    with rio.open(resultfiles+f"/maps/NO2_moy_after.tif", "w", **kwargs) as dst:
        dst.write(NO2_moy_after.copy(), 1)
    with rio.open(resultfiles+f"/maps/NO2_moy_before.tif", "w", **kwargs) as dst:
        dst.write(NO2_moy_before.copy(), 1)
    with rio.open(resultfiles+f"/maps/NO2_moy_diff.tif", "w", **kwargs) as dst:
        dst.write(diff.copy(), 1)
    with rio.open(resultfiles+f"/maps/NO2_moy_rel_dif.tif", "w", **kwargs) as dst:
        dst.write(rel_dif.copy(), 1)


network = gpd.read_file(f"{mod_p.pathToSirane}/{mod_p.inputs.FICH_DIR_RESUL}/reseau-rues.shp")

fig,axs = map_subplots(2,2)
ax = axs.flatten()
map0 = ax[0].imshow(NO2_moy_before, extent=ext, vmax=45, cmap="RdYlGn_r")
plt.colorbar(map0,ax=ax[0], extend="max").set_label( label="Concentration moyenne annuelle\nen NO₂ avant [µg/m³]", size=15)
network.plot(ax=ax[0], edgecolor="k", linewidth=0.5)
map1 = ax[1].imshow(NO2_moy_after, extent=ext, vmax=45, cmap="RdYlGn_r")
plt.colorbar(map1,ax=ax[1], extend="max").set_label( label="Concentration moyenne annuelle\nen NO₂ après [µg/m³]", size=15)
network.plot(ax=ax[1], edgecolor="k", linewidth=0.5)
map2 = ax[2].imshow(diff, extent=ext, cmap="coolwarm", vmin = -max_diff, vmax=max_diff)
plt.colorbar(map2,ax=ax[2], ).set_label( label="Différence des concentrations\n moyennes annuelles [µg/m³]", size=15)
network.plot(ax=ax[2], edgecolor="k", linewidth=0.5)
map3 = ax[3].imshow(rel_dif, extent=ext, cmap="coolwarm", vmin = -max_rel_dif, vmax=max_rel_dif)
plt.colorbar(map3,ax=ax[3]).set_label(label="Différence relative des concentrations\n moyennes annuelles [%]", size=15)
network.plot(ax=ax[3], edgecolor="k", linewidth=0.5)
plt.savefig(resultfiles+"MoyennesAnnuelles.png")
plt.close("all")


norm = colors.BoundaryNorm(boundaries=[0,10,40,60], ncolors=256)
fig,axs = map_subplots(1,2)
m1 = axs[0].imshow(NO2_moy_before, extent = ext, norm=norm, cmap='RdBu_r')
plt.colorbar(m1, ax=axs[0], shrink=0.45,label="Concentration moyenne annuelle\nen NO₂ avant [µg/m³]")
m2 = axs[1].imshow(NO2_moy_after, extent = ext, norm=norm, cmap='RdBu_r')
plt.colorbar(m2, ax=axs[1],shrink=0.4
5, label="Concentration moyenne annuelle\nen NO₂ après [µg/m³]")
plt.savefig(resultfiles+"MoyennesAnnuelles2.png")
plt.close("all")


#Moyenne Annuelle en heure d'ouvertures des écoles
grids_path["Heure"] = grids_path.index.hour
grids_path["weekday"] = grids_path.index.dayofweek
OpenHours = grids_path.reset_index()
OpenHours.set_index("index").Heure.plot()
OpenHoursID = list(OpenHours[(OpenHours.Heure >= 7)&(OpenHours.Heure < 18)&(OpenHours.weekday < 5)].index)

NO2_moy_after_OH =ultrabig_data2[:,:,OpenHoursID].mean(axis=2)
NO2_moy_before_OH = ultrabig_data[:,:,OpenHoursID].mean(axis=2)
diff_OH = (NO2_moy_after_OH-NO2_moy_before_OH)
max_diff_OH = np.abs(diff_OH).max()
rel_dif_OH = diff_OH/NO2_moy_before_OH * 100
max_rel_dif_OH = np.abs(rel_dif_OH).max()

fig,axs = map_subplots(2,2)
ax = axs.flatten()
map0 = ax[0].imshow(NO2_moy_before_OH, extent=ext, vmax=50, cmap="RdYlGn_r")
plt.colorbar(map0,ax=ax[0], extend="max").set_label( label="Concentration moyenne annuelle\nen NO₂ avant [µg/m³]", size=15)
network.plot(ax=ax[0], edgecolor="k", linewidth=0.5)
map1 = ax[1].imshow(NO2_moy_after_OH, extent=ext, vmax=50, cmap="RdYlGn_r")
plt.colorbar(map1,ax=ax[1], extend="max").set_label( label="Concentration moyenne annuelle\nen NO₂ après [µg/m³]", size=15)
network.plot(ax=ax[1], edgecolor="k", linewidth=0.5)
map2 = ax[2].imshow(diff_OH, extent=ext, cmap="coolwarm", vmin = -max_diff_OH, vmax=max_diff_OH)
plt.colorbar(map2,ax=ax[2], ).set_label( label="Différence des concentrations\n moyennes annuelles [µg/m³]", size=15)
network.plot(ax=ax[2], edgecolor="k", linewidth=0.5)
map3 = ax[3].imshow(rel_dif_OH, extent=ext, cmap="coolwarm", vmin = -max_rel_dif_OH, vmax=max_rel_dif_OH)
plt.colorbar(map3,ax=ax[3]).set_label(label="Différence relative des concentrations\n moyennes annuelles [%]", size=15)
network.plot(ax=ax[3], edgecolor="k", linewidth=0.5)
plt.savefig(resultfiles+"MoyennesAnnuelles_7h-18h.png")
plt.close("all")


####Effet long terme sur la santé
#Import population
pop = gpd.read_file(mainfile+"densite-de-population_2018_quartier.geojson")
pop["surface"] = pop.area/10**6
pop["Popu"] = pop.value * pop["surface"]
pop[pop.name=="Josaphat"]

(left, right, bottom, top) = [int(v) for v in ext]
x = np.linspace(left, right,diff.shape[1], endpoint=False)
x += (x[1]-x[0])/2
y = np.linspace(bottom, top,diff.shape[0], endpoint=False)
y += (y[1]-y[0])/2

from shapely.geometry import Point
a = gpd.GeoDataFrame(geometry=[Point(xv,yv) for xv in x for yv in y], crs="EPSG:31370")
a_sj = a.sjoin(pop)
a_sj["X"] = ((a_sj.geometry.x - (x[1]-x[0])/2 - left)/(x[1]-x[0])).astype(int)
a_sj["Y"] = ((top-a_sj.geometry.y - (y[1]-y[0])/2)/(y[1]-y[0])).astype(int)
data_avant = {}
data_après = {}
for shp in set(pop["name"]):
    data_avant[shp] = ultrabig_data[a_sj[a_sj["name"] == shp].Y,a_sj[a_sj["name"] == shp].X,:].mean(axis=0)
    data_après[shp] = ultrabig_data2[a_sj[a_sj["name"] == shp].Y,a_sj[a_sj["name"] == shp].X,:].mean(axis=0)

data = pd.concat([pd.DataFrame(data_avant, index = grids_path.index).stack(),pd.DataFrame(data_après, index = grids_path2.index).stack()],
                  axis = 1, keys=["Avant", "Après"])
data = data.reset_index().rename(columns={"level_0":"Date","level_1":"Quartier"})
data["Day"] = data.Date.dt.date
data.groupby(["Quartier", "Day"]).mean()

#Nbr dépassement horaire de 200µg/m3 :
pic_C_av = (ultrabig_data>= 200).sum(axis=2)
pic_C_ap = (ultrabig_data2>= 200).sum(axis=2)
(pic_C_ap-pic_C_av).max()
(pic_C_ap-pic_C_av).min()


fig,axs = map_subplots(1,3)
m1 = axs[0].imshow(pic_C_av, extent = ext)
rues.plot(ax=axs[0], color="grey")
plt.colorbar(m1, ax=axs[0])
m2 = axs[1].imshow(pic_C_ap, extent = ext)
plt.colorbar(m2, ax=axs[1])
m3 = axs[2].imshow(pic_C_ap-pic_C_av, extent = ext)
plt.savefig(resultfiles+"nbrdepassement200ugm3.png")
plt.close("all")

#Limite journaliaire de 25µg/m3 :
grids_path["id"] = range(len(grids_path))
grids_path["day"] =  grids_path.index.date
idperdays = grids_path.groupby("day")["id"].apply(lambda df: df.values)
if not os.path.isfile(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/../meanday.npz"):
    meanday_av = np.ndarray(shape=(ultrabig_data.shape[0],ultrabig_data.shape[1],len(idperdays)), dtype=np.float16)
    meanday_ap = np.ndarray(shape=(ultrabig_data.shape[0],ultrabig_data.shape[1],len(idperdays)), dtype=np.float16)
    for i,ids in enumerate(idperdays):
        print(i)
        meanday_av[:,:,i] = ultrabig_data[:,:,ids].mean(axis=2)
        meanday_ap[:,:,i] = ultrabig_data2[:,:,ids].mean(axis=2)
        np.savez(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/../meanday.npz", meanday_av = meanday_av, meanday_ap=meanday_ap, dates= idperdays.index)
    else:
        with np.load(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/../meanday.npz",allow_pickle=True) as data:
            meanday_av = data["meanday_av"]
            meanday_ap = data["meanday_ap"]
            dates = data["dates"]

            meanday_C_av = (meanday_av> 25).sum(axis=2)
            meanday_C_ap = (meanday_ap> 25).sum(axis=2)
            (meanday_C_ap-meanday_C_av).max()
            (meanday_C_ap-meanday_C_av).min()


rues.plot(ax=axs[0], color="grey")
fig,axs = map_subplots(1,3)
m1 = axs[0].imshow(meanday_C_av, extent = ext)
plt.colorbar(m1, ax=axs[0], shrink=0.45,label="Nombre de jour dépasant une concentration\nmoyenne de 25µg/m³ (Avant)")
m2 = axs[1].imshow(meanday_C_ap, extent = ext)
plt.colorbar(m2, ax=axs[1],shrink=0.45, label="Nombre de jour dépasant une concentration\nmoyenne de 25µg/m³ (Après)")
m3 = axs[2].imshow(meanday_C_ap-meanday_C_av, extent = ext)
plt.colorbar(m3, ax=axs[2],shrink=0.45, label="Nombre de jour dépasant une concentration\nmoyenne de 25µg/m³ (Après-Avant)")
plt.savefig(resultfiles+"nbrdepassement25ugm3.png")
plt.close("all")

#Limite Annuelle de 40µg/m3 : (EU)
meanYear_C_av = (ultrabig_data.mean(axis=2)> 40)
meanYear_C_ap = (ultrabig_data2.mean(axis=2)> 40)
meanYear_C_av.min()
meanYear_C_ap.min()
fig,axs = map_subplots(1,3)
m1 = axs[0].imshow(meanYear_C_av, extent = ext)
plt.colorbar(m1, ax=axs[0], shrink=0.45,label="Lieux dépasant une concentration\nmoyenne annuelle de 40µg/m³ (Avant)")
m2 = axs[1].imshow(meanYear_C_ap, extent = ext)
plt.colorbar(m2, ax=axs[1],shrink=0.45, label="Lieux dépasant une concentration\nmoyenne annuelle de 40µg/m³ (Après)")
m3 = axs[2].imshow(meanYear_C_ap-meanYear_C_av, extent = ext)
plt.colorbar(m3, ax=axs[2],shrink=0.45, label="Lieux dépasant une concentration\nmoyenne annuelle de 40µg/m³ (Après-Avant)")
plt.savefig(resultfiles+"nbrdepassement40ugm3.png")
plt.close("all")

#Limite Annuelle de 10µg/m3 : (OMS)
meanYearOMS_C_av = (ultrabig_data.mean(axis=2)> 10)
meanYearOMS_C_ap = (ultrabig_data2.mean(axis=2)> 10)
meanYearOMS_C_av.min()
meanYearOMS_C_ap.min()


####Effet courts terme sur la santé
if not os.path.isfile(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/../quant95.npz"):
    q_ap = np.quantile(ultrabig_data2,0.95, axis=2)
    q_av = np.quantile(ultrabig_data ,0.95, axis=2)
    np.savez(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/../quant95.npz", q_ap=q_ap, q_av=q_av)
else:
    with np.load(f"{mod_v.pathToSirane}/{mod_v.inputs.FICH_DIR_RESUL}/../quant95.npz",allow_pickle=True) as data:
        q_ap = data["q_ap"]
        q_av = data["q_av"]

diff=q_ap-q_av
max_diff = np.abs(diff).max()
rel_dif = diff/q_av * 100
max_rel_dif = np.abs(rel_dif).max()
fig,axs = map_subplots(2,2)
axs = axs.flatten()
m1 = axs[0].imshow(q_av, extent = ext)
plt.colorbar(m1, ax=axs[0], shrink=1,label="Quantile 95% des concentrations de \nNO₂ [µg/m³] (Avant)")
m2 = axs[1].imshow(q_ap, extent = ext)
plt.colorbar(m2, ax=axs[1],shrink=1, label="Quantile 95% des concentrations de \nNO₂ [µg/m³] (Après)")
m3 = axs[2].imshow(diff, extent = ext, cmap="coolwarm", vmin = -max_diff, vmax=max_diff)
plt.colorbar(m3, ax=axs[2],shrink=1, label="Difference des quantiles 95% des concentrations de \nNO₂ [µg/m³] (Après-Avant)")
m4 = axs[3].imshow(rel_dif, extent = ext, cmap="coolwarm", vmin = -max_rel_dif, vmax=max_rel_dif)
plt.colorbar(m4, ax=axs[3],shrink=1, label="Difference relative des quantiles 95% des concentrations de \nNO₂ [%] (Après-Avant)/Avant")
plt.savefig(resultfiles+"Quantils0.95.png")
plt.close("all")
