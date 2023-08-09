"""
Fichier pour calculer les émissions de surface.

Il se base sur la météo calculé avec le fichier meteo.py et sur la répartition du batit du travail préparatoire de A.Briffaut.txt
Il faut ces deux fichiers (météo_new et EMIS_ring.txt) dans le dossier Préprocessing
"""
#import des données météos
mainfile = "Preprocessing/"
data_meteo = pd.read_csv(f"{mainfile}/Meteo_New.dat", sep="\t", index_col=0, na_values=-9999)
data_meteo.index = pd.to_datetime(data_meteo.index, format="%d/%m/%Y %H:%M")
data_meteo["Date-day"] = data_meteo.index.date

#Calcule des degré-jours équivalent
data_DJ = data_meteo.groupby("Date-day")[["Temp"]].mean()
data_DJ["Dj"] = (15 - data_DJ["Temp"]).clip(0,10000) #eq 12
data_DJ["Dj_eq"] = data_DJ["Dj"].rolling(3).apply(lambda x: 0.6*x[0] + 0.3*x[1] + 0.1*x[2]) #eq 11
data_meteo["Dj_eq"] = data_meteo["Date-day"].map(data_DJ["Dj_eq"])
data_meteo["Year"] = data_meteo.index.year

#Répartition des émissions totales selon les degré-jours équivalent
ESurf_attendu = pd.Series({2010:2.2979169405472, 2011: 1.80005239053116, 2012: 1.95466673174922,
                                     2013:2.01859403082654, 2014: 1.59518184243388, 2015: 1.71492074367155,
                                     2016:1.69382160652324, 2017: 1.62871586439113, 2018: 1.58872590747962,
                                     2019:1.57634946546155}) #kt/ans
data_meteo["Coeff_Modul_kt_h"] = (data_meteo["Dj_eq"] / data_meteo.groupby("Year")["Dj_eq"].transform("sum") * data_meteo.Year.map(ESurf_attendu)) #kt/h
#transforme en g/s
data_meteo["Coeff_Modul"]  = data_meteo.Coeff_Modul_kt_h / 3600 * 10**9 #g/s = kt/h / (s/h) * (g/kt)
data_meteo[["Coeff_Modul"]].to_csv(f"{mainfile}/Mod_Temp_Cadastre.dat", na_rep=-9999, sep="\t",date_format="%d/%m/%Y %H:%M")

#Normalisation de la répartition spatiale
E_surface = pd.read_csv(f"{mainfile}/EMIS_ring.txt", sep="\t")
E_surface["NO"] = E_surface["NO"]  / E_surface["NO"].sum() * 0.95 #-
E_surface["NO2"]= E_surface["NO2"] / E_surface["NO2"].sum()* 0.05 #-
E_surface["NOandNO2"] = E_surface[["NO","NO2"]].sum(axis=1)
E_surface.sum()
E_surface[["X","Y","NO","NO2","O3"]].to_csv(f"{mainfile}/EMIS_ring_New.txt", sep="\t", index=False, na_rep=-9999)
