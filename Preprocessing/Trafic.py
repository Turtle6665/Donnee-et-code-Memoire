"""
Ce script montre comment j'ai modifié la dynamique de trafic pour la modulation hebdomadaire.

Les données de comptage du trafic sont à mettre dans le fichier préprocessing.

Pour les différentes espèces de polluants, relancer le script en changant la variable "Espece".
"""

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pytz
SIRANE_tz = pytz.FixedOffset(+1*60) #GMT+1 : le fuseau horaire de sirane
mainfile = "Preprocessing/"
Espece = "NO2"

# Création de la plage de dates pour l'année 2018 avec une fréquence horaire
date_range = pd.date_range(start='2018-01-01', end='2018-12-31 23:00:00', freq='H', tz="CET") #fuseau horaire de la belgique AVEC changement d'heure (GMT+1 hiver et GMT +2 été)

# Initialisation de la série avec des valeurs nulles
values = pd.Series(np.nan, index=date_range)
values.index.name = "Date"
values[24*83:24*83+4]
values[24*300:24*300+4]
# Remplacement des valeurs pour chaque jour de la semaine
values[values.index.weekday < 5] = 1  # Jours de semaine
values[values.index.weekday == 5] = 2  # Samedis
values[values.index.weekday == 6] = 3  # DImanche
# Remplacement des valeurs pour les jours février belges
belgian_holidays = [
    '2018-01-01',
    '2018-04-02',
    '2018-05-01',
    '2018-05-08',
    '2018-05-10',
    '2018-05-21',
    '2018-07-14',
    '2018-08-15',
    '2018-11-01',
    '2018-11-11',
    '2018-12-25'
]
belgian_holidays = [datetime.strptime(date_str, '%Y-%m-%d') for date_str in belgian_holidays]
for holiday in belgian_holidays:
    values.loc[values.index.date == holiday.date()] = 3  # Jours février belges

values.head()

df = pd.DataFrame(values, columns=["TypeJour"])
df["hourID"] = (df["TypeJour"]*24 + df.index.hour).astype(int)

# Application du nbr de voiture par type de jours
ls = []
for i in ["semaine","samedi","ferie_dimanche"]:
    ls.append(pd.read_csv(f"Complet2018_{i}.csv", sep=";"))

Sem=pd.concat(ls)
len(set(Sem.Poste))
count = Sem.groupby(["Poste","TypeJour"]).count().Heure
postes = count[count==24].unstack().dropna().index

data = Sem.loc[Sem.Poste.apply(lambda x : x in postes),:].groupby(["TypeJour","Heure"])[["Volume"]].sum().reset_index()
data["TypeJour"] = data.TypeJour.map({"JOS":1,"SAM":2,"DIM":3})
data["hourID"] = data["TypeJour"]*24 + data["Heure"]
data_map = data.sort_values("hourID").set_index("hourID")["Volume"]


df["Coef"] = df.hourID.map(data_map)/df.hourID.map(data_map).sum()
df.index = df.index.tz_convert(SIRANE_tz) #une fois les jours remplis, on reconvertit en GMT +1 pour bien avoir la dynamique décalée d'une heure entre été/hivers

##import de la série temporelle de base
Emmissions = pd.read_csv(f"{mainfile}/../SIRANE/INPUT/EMISSIONS/EMIS_LIN/Mod_Temp_Trafic{Espece}.dat", sep="\t")
Emmissions.index = pd.to_datetime(Emmissions["Date"], format="%d/%m/%Y %H:%M")
Emmissions.index=Emmissions.index.tz_localize(SIRANE_tz) #Données en GMT+1
Emmissions = Emmissions["2018":"2018"]
Emmissions["mois"] = Emmissions.index.month
Emmissions["jour_sem"]= Emmissions.index.day_of_week
Emmissions["TypeJour"]=Emmissions["jour_sem"].map({6:"DIM",5:"SAM"})
Emmissions.loc[Emmissions.TypeJour.isna(),"TypeJour"]="JOS"
Emmissions["moyenne_mois"] = Emmissions.mois.map(Emmissions.groupby(["mois","TypeJour"])["Coeff_Modul"].mean().reset_index("mois").loc["JOS",:].set_index("mois")["Coeff_Modul"]/Emmissions.groupby("TypeJour")["Coeff_Modul"].mean().loc["JOS"])

#Modulation des données de comptage en fonction des données initiales (pour avoir la dynamique des mois)
df["Coeff_Modul"] = df["Coef"]/df["Coef"].sum()*Emmissions["Coeff_Modul"].sum() * Emmissions["moyenne_mois"]/Emmissions["moyenne_mois"].mean()
df["Coeff_Modul_axel"] = Emmissions["Coeff_Modul"]

df[["Coeff_Modul","Coeff_Modul_axel"]].sum() #la "petite" différence est due aux jours ferier il me semble.

df[["Coeff_Modul","Coeff_Modul_axel"]].plot(marker="o")
plt.show()
df["Coeff_Modul"].to_csv(f"{mainfile}/Mod_Temp_Trafic{Espece}.dat",sep="\t",index=True, date_format='%d/%m/%Y %H:%M', na_rep=-9999)
