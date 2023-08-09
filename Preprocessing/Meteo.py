"""
Ce script montre comment j'ai agrégée par heure les données météorologiques de l'IRM.

Les données météo issue de l'IRM sont a mettre dans le fichier "Preprocessing".
#"""

import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pytz
SIRANE_tz = pytz.FixedOffset(+1*60) #GMT+1 : le fuseau horaire de sirane
Data_tz = "CET"
mainfile = "Preprocessing/"


ls = []
for file in [f"{mainfile}Meteo_Uccle_2010-2019.csv",f"{mainfile}Meteo_Uccle_2020-2021.csv"]: #
    dataIRM = pd.read_csv(f"{mainfile}/{file}", sep=";", skiprows=9)
    dataIRM = dataIRM[dataIRM.PARAMETER_ELEVATION!=1.5] #pour enlever la T° à 1.5m du sol
    dataIRM = dataIRM[~((dataIRM.PARAMETER_ELEVATION==2)&(dataIRM.PARAMETER=="WIND_SPEED"))] #pour enlever la vitesse du vent à 2 du sol
    dataIRM["Date"] = dataIRM.DATE + " " + dataIRM.HOUR
    dataIRM["d"] = dataIRM.groupby(["Date","STEP","PARAMETER"])["VALUE"].transform(lambda x: range(len(x))) #est égale à 0 sauf pour les deuxiemmes données lors du changement d'heure d'autoomne. Je suppose que les données de la 2e 02:00, sont après celles de la 1ere
    dataIRM_table = dataIRM[["PARAMETER","VALUE","Date","d"]].reset_index().pivot_table(index=["Date","d"], columns="PARAMETER", values="VALUE", aggfunc=np.mean) #Moyenne données des demis-heures
    dataIRM_table["PRECIPITATION"] = dataIRM_table["PRECIPITATION"]*2 #car la somme des précipitation est le double de la moyenne (2 données/h)
    dataIRM_table = dataIRM_table.reset_index()
    dataIRM_table.index = pd.to_datetime(dataIRM_table.Date, format = "%d/%m/%Y %H:%M")
    dataIRM_table.index = dataIRM_table.index.tz_localize(Data_tz, ambiguous="infer") #ambiguous="infer" will attempt to infer fall dst-transition hours based on order
    dataIRM_table = dataIRM_table[['GLOBAL_RADIATION', 'PRECIPITATION', 'PRESSURE','TEMPERATURE', 'WIND_DIRECTION', 'WIND_SPEED']]
    dataIRM_table= dataIRM_table.rename(columns={"GLOBAL_RADIATION":"SolRad","PRECIPITATION":"Precip","TEMPERATURE":"Temp","WIND_DIRECTION":"Dir", "WIND_SPEED":"U"})
    ls.append(dataIRM_table)


file=f"{mainfile}Meteo_Uccle_2021-2023.csv"
dataIRM = pd.read_csv(f"{mainfile}/{file}", sep=";", header=[0,1], index_col=[0])
dataIRM = dataIRM.iloc[:-4,:-2].astype(float)#drop les 4 dernières lignes étranges et les deux dernières colones
dataIRM[dataIRM == -9999] = np.nan #enlever les nan
dataIRM["SolRad"] = dataIRM["DSRIRM"]["J/cm2"].clip(0,100000000000) * 10**4/3600 #convertire les J/cm2 en W/m2
dataIRM.index = pd.to_datetime(dataIRM.index, format = "%d/%m/%Y %H:%M")
dataIRM.index = dataIRM.index.tz_localize(SIRANE_tz) + timedelta(hours=0.5) #+ timedelta(hours=-0.5) #J'ai l'impression que lui est déjà en GMT+1 et avec un décalage de 30 min ?

dataIRM = dataIRM.droplevel(1,axis=1) #elever les unitées
dataIRM["Date"] = pd.Series(dataIRM.index, index=dataIRM.index).astype(str).apply(lambda x : x[:-11]+"00") #pour ne plus avoir les demis-heures
dataIRM[["SolRad","RAINRM"]] = dataIRM[["SolRad","RAINRM"]]*2 #on fait fois 2 comme ca quand on fait la moyenne c'est comme la somme. La somme des irradiance car elle étaient en J recue par demie heure
dataIRM_table = dataIRM.groupby(["Date"]).mean()
dataIRM_table = dataIRM_table.rename(columns={"SolRad":"SolRad","RAINRM":"Precip","TT0IRM":"Temp","WD0IRM":"Dir", "WSHIRM":"U"}) #pour la vitesse du vent, j'ai pris celle qui sont les plus grandes, je suppose qu'il y a encore deux hauteur de mesure et que la plus grande est celle à 30m
dataIRM_table[["U","Dir","Temp","SolRad","Precip"]]
dataIRM_table.index = pd.to_datetime(dataIRM_table.index, format = "%Y-%m-%d %H:%M").tz_localize(SIRANE_tz) #tjs en GMT+1 et je fait ça pour avoir tout en format datetime
ls.append(dataIRM_table)

for df in ls:
    df.index = df.index.tz_convert(SIRANE_tz) #une fois les jours remplis, on reconvertit toutes les timeszones en GMT +1

full_data = pd.concat(ls)
###Pour avoir toutes les heures des données même dans celles ou il manque 100%
date_range = pd.date_range(start='2010-01-01', end='2023-12-31 23:00:00', freq='H', tz=SIRANE_tz) # Création de la plage de dates avec une fréquence horaire: fuseau horaire GMT+1
values = pd.Series(np.nan, index=date_range) #avoir une série avec touts les indexs mais des données absantes
values.index.name = "Date"
full_data = pd.concat([values,full_data], axis=1)[["U","Dir","Temp","SolRad","Precip"]] #lier les deux pour avoir toutes les heurs dans les données
#full_data.index= full_data.index + timedelta(hours=1)

full_data = full_data.sort_index() #pour être sur qu'elles soient dans l'ordre alphabétique
full_data.plot(subplots=True)
plt.show()

full_data.to_csv(f"{mainfile}/Meteo_New.dat",sep="\t",index=True, date_format='%d/%m/%Y %H:%M', na_rep=-9999)
