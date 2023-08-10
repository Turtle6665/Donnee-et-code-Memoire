"""
Ce script permet de lancer et parametrer SIRANE depuis un script python.
/!\ Une version plus à jour est disponible dans le github dédiée.
"""

import subprocess
from os import makedirs
from os.path import exists,isfile

import io,glob
import pandas as pd
import numpy as np
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
try:
    from osgeo import gdal
except:
    print("Gdal not found, please install it to use transph_grd() !")

pathToSirane = "SIRANE/"
PySiWrapPath = "PySiWrap"

def get_input_keys(input_key_file):
    with io.open(f"{input_key_file}",'r',encoding='utf8') as f:
        data = f.read()
    data = pd.Series(data.split(sep="\n"))
    def sepdf(str_val):
        str_val = str_val.split(";")
        if len(str_val) != 7:
            None #return ["","","","","","",""]
        else :
            return str_val
        None #return 'oups should not be there'
    strip = lambda x : x.strip()
    sprip_v = np.vectorize(strip) #to remove spaces in front and at the back
    data = sprip_v(np.array(data.apply(sepdf).dropna().to_list()))
    data = pd.DataFrame(data[1:,1:], index = data[1:,0], columns = data[0,1:])
    data["Key"] = data.index
    return data

def read_dat(file_path, Nindex = 0, Ncol = 0, mod_input = None, mod_result = None):
    """
    fonction pour lire les séries temporelles en .dat
    INPUT :
        file_path[str]: the path to the .dat file
        Nindex[int]: the index column number
        Ncol[int]: the colomns names number
        mod_input = None [model]: the model so that the file_path is relative to the input file
        mod_result = None [Model]:the model so that the file_path is relative to the result file
    OUTPUT :
        the data table as a df
    """
    if mod_input is not None:
        file_path = f"{mod_input.pathToSirane}/{mod_input.inputs.FICH_DIR_INPUT}/{file_path}"
    if mod_result is not None:
        file_path = f"{mod_result.pathToSirane}/{mod_result.inputs.FICH_DIR_RESUL}/{file_path}"

    data = pd.read_csv(file_path, sep='\t',index_col=Nindex,header=Ncol,na_values=["ND", "-9999", "-999"])
    if type(data.index[0]) is str:
        if "/" in data.index[0] and ":" in data.index[0]:
            data.index = pd.to_datetime(data.index, format="%d/%m/%Y %H:%M")
    return data

class Model(object):
    """
        INPUT :
            name [str]:
            pathToSirane = {pathToSirane} [str]: the relative path to the folder containing SIRAN.exe
            InputFile [str] : the donnees.dat file relative to the input folder
    """

    def __init__(self, name, pathToSirane=pathToSirane, InputFile = "Donnees.dat"):
        super(Model, self).__init__()
        self._lst_modifiable_inputs = []
        self._lst_changes_inputs=[]
        self.name = name
        self.pathToSirane = pathToSirane
        self.InputFile = InputFile
        self._inputs_full = self._read_input(self.InputFile)
        self.inputs = self._inputs_full.Value

        self._site_meteo_full = self._read_site(self.inputs.FICH_SITE_METEO)
        self.site_meteo = self._site_meteo_full.Value
        self._site_disp_full = self._read_site(self.inputs.FICH_SITE_DISP)
        self.site_disp = self._site_disp_full.Value

    def __repr__(self):
        return f"PySiWrap.Model named '{self.name}' at '{self.pathToSirane}'"

    def run(self, silent= False):
        """
        Fonction pour lancer Sirane depuis un fichier python
        INPUT :
            silent [bool = ]
        OUTPUT :
            Bool : True = the model has been lunched; false = an error occured due to inputs
        """

        #Safeguard if the inputs have not been well set
        try:
            self._inputs_full.Value == self.inputs
            self.site_meteo == self._site_meteo_full.Value
            self.site_disp == self._site_disp_full.Value
        except :
            reflist = pd.concat([self._inputs_full, self._site_meteo_full, self._site_disp_full])
            modlist = pd.concat([self.inputs, self.site_meteo, self.site_disp])
            isnotinlist = np.vectorize(lambda x : not x in reflist.index)
            print(f"You cannot add an imput that as not been set in the default parameter ! Those are not present: {modlist.index[isnotinlist(modlist.index)].to_list()}")
            return False

        self._write_all_dat()
        self._write_input(self._inputs_full, self.InputFile)
        self._write_site(self._site_meteo_full, self.inputs.FICH_SITE_METEO)
        self._write_site(self._site_disp_full, self.inputs.FICH_SITE_DISP)

        self.log_file_path = f"{self.inputs.FICH_DIR_RESUL}/log.txt"
        self.PathToInput = f"{self.inputs.FICH_DIR_INPUT}/{self.InputFile}"
        makedirs(f"{self.pathToSirane}/{self.inputs.FICH_DIR_RESUL}", exist_ok=True)
        #make a copy of all the input file that have changed in the result file
        self._export_Mod_input_files(f"{self.pathToSirane}/{self.inputs.FICH_DIR_RESUL}/INPUT/")
        if silent==False: #run in an other terminal
            a = subprocess.run(["start", f"SIRANE - {self.name}","/wait", "cmd", f"/k cd {self.pathToSirane} & sirane {self.PathToInput} {self.log_file_path} & exit /b"], shell = True, capture_output=True)
        elif silent=="Colab":
            os.chdir(self.pathToSirane)
            a = subprocess.run(["wine", "sirane.exe", self.PathToInput, self.log_file_path], capture_output=True)
            os.chdir("..")
        else :
            a = subprocess.run(["cd",self.pathToSirane, "&", "sirane", self.PathToInput, self.log_file_path], shell = True, capture_output=True)

        #back to default parameter on all the input files
        self._Del_Mod_input_file()
        return not "ERREUR FATALE" in self.get_log(printlogs=False)

    def get_log(self, printlogs=True):
        with io.open(f"{self.pathToSirane}/{self.log_file_path}",'r',encoding='utf8') as f:
            logs = f.read()
        if printlogs==True:
            print(logs)
            return None
        else:
            return logs

    """
        Partie liée a la lecture et l'ecriture du fichier d'inputs de Sirane
    """
    def _read_input(self,input_file):
        """
        Fonction pour lire les données en python
        INPUT :
            input_file [str] : the input file relative to the input folder
        OUTPUT :
            The input file data as a np array.
        """
        with io.open(f"{self.pathToSirane}/INPUT/{input_file}",'r',encoding='utf8') as f:
            data = f.read()
        data = pd.Series(data.split(sep="\n"))
        def sepdf(str_val):
            if "=" in str_val:
                return str_val.split("=")
            else :
                return [str_val,""]
            return 'oups should not be there'
        strip = lambda x : x.strip()
        sprip_v = np.vectorize(strip) #to remove spaces in front and at the back
        data = sprip_v(np.array(data.apply(sepdf).to_list()))
        data = pd.DataFrame(data, columns=["Name","Value"])
        # add the keyyysss
        keys = get_input_keys(f"{PySiWrapPath}/DefaultParam/Don_Defaut_FR.dat")
        data = data.merge(keys, how="right", right_on="Description", left_on="Name").reindex(columns=["Description","Key","Value"])
        data.index = data.Key
        del data["Key"]
        return data

    def _write_input(self, data, input_dat, keepDef=True):
        """
        Fonction pour écrire le fichier de données
        INPUT :
            data [np.array of str] : the data that should be writen in the file
            input_dat [str] : the input file relative to the input folder
        OUTPUT :
            None
        """
        if keepDef==True:
            if input_dat in self._get_Mod_input_files():
                None
            else :
                self._Add_Mod_input_files(input_dat)
        data.dropna().to_csv(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{input_dat}",sep="=", header = False, index=False)
        return None

    """
        Partie liée a la lecture et l'ecriture des fichiers de site (disp et meteo)
    """
    def _read_site(self,site_file):
        """
        Fonction pour lire le fichier de site
        INPUT :
            site_file [str] : the site file relative to the input folder
        OUTPUT :
            The site file data as a np array.
        """
        with io.open(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{site_file}",'r',encoding='utf8') as f:
            data = f.read()
        data = pd.Series(data.split(sep="\n"))
        def sepdf(str_val):
            if "=" in str_val:
                return str_val.split("=")
            else :
                return [str_val,""]
            return 'oups should not be there'
        strip = lambda x : x.strip()
        sprip_v = np.vectorize(strip) #to remove spaces in front and at the back
        data = sprip_v(np.array(data.apply(sepdf).to_list()))
        data = pd.DataFrame(data, columns=["Name","Value"])
        # add the keyyysss
        keys = get_input_keys(f"{PySiWrapPath}/DefaultParam/Site_Defaut_FR.dat")
        data = keys.merge(data, how="left", left_on="Description", right_on="Name").reindex(columns=["Description","Key","Value"])
        data.index = data.Key
        del data["Key"]
        #data.Value = data.Value.astype(float)
        return data

    def _write_site(self, data, site_file, keepDef=True):
        """
        Fonction pour écrire le fichier de site
        INPUT :
            data [df] : the data that should be writen in the file
            site_file [str] : the site file relative to the input folder
        OUTPUT :
            None
        """
        if keepDef==True:
            if site_file in self._get_Mod_input_files():
                None
            else :
                self._Add_Mod_input_files(site_file)
        data.dropna().to_csv(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{site_file}",sep="=", header = False, index=False)
        return None

    """
        This part alow to copy a data file to a temporary file and therfore keep a default value of the data
    """
    def _Add_Mod_input_files(self, PathToFile):
        def add(PathToFile):
            """
                Make a safly modfifiable files from the PathToFile. It can be a string of the path or a list of stings
            """
            PathToFile = str(PathToFile)
            if PathToFile in self._lst_modifiable_inputs:
                raise Exception(f"{PathToFile} is already in a modfifiable way")
            elif exists(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{PathToFile}-def")  :
                raise Exception(f"{PathToFile}-def already exist but {PathToFile} is not in a modfifiable way. Use 'PySiWrap.force_all_files_to_def(self)' to change to default file")
            else :
                shutil.copy(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{PathToFile}", f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{PathToFile}-def")
                self._lst_modifiable_inputs.append(PathToFile)
        if type(PathToFile) is list:
            for path in PathToFile:
                add(path)
        else :
            add(PathToFile)

    def _Del_Mod_input_file(self, PathToFile=None, force_to_def = False):
        """
            put back to default the safly modfifiable files and make it not safly modfifiable
            if PathToFile = None -> change all the files
            if force_to_def = True -> change the file to default even if not registerd as a safly motifiable
        """
        def rem(path):
            path = str(path)
            if path in self._lst_modifiable_inputs:
                #replace file
                shutil.move(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{path}-def",f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{path}")
                self._lst_modifiable_inputs.remove(path)
                return None
            elif force_to_def :
                shutil.move(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{path}-def",f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{path}")
                return None
            else:
                raise Exception(f"{path} is not in a modfifiable way")

        if PathToFile is None:
            PathToFile = self._get_Mod_input_files()

        if type(PathToFile) is list:
            for path in PathToFile:
                rem(path)
        else :
            rem(PathToFile)

    def _get_Mod_input_files(self):
        """To get all path to the safly modfifiable files already changed"""
        return self._lst_modifiable_inputs.copy()

    def get_Mod_changes_inputs(self, data=False):
        """
        To get all path to the safly modfifiable files that will be changed
        INPUT:
            data[bool]: if true, return all the safly modfifiable files and not only the paths
        OUTPUT:
            data = True:
                List of [data, file_path, sep, index]
            else :
                list of file_path
        """
        if data:
            return self._lst_changes_inputs.copy()
        else:
            return [a[1] for a in self._lst_changes_inputs].copy()

    def _export_Mod_input_files(self, pathToExport, files = None):
        """to export the safly modifiable files to the file pathToExport (if it doesn't exist, it will create a new and if it exist, it removes everything from it!)"""
        if files is None:
            files = self._get_Mod_input_files()
        #rmv all the files in the path to export
        shutil.rmtree(pathToExport, ignore_errors=True)
        def exp(PathToFile):
            dir = "".join([str(i)+"/" for i in f"{pathToExport}/{PathToFile}".split("/")[:-1]])
            if dir != "/" :
                makedirs(dir, exist_ok=True)
            shutil.copy(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{PathToFile}", f"{pathToExport}/{PathToFile}")

        if type(files) is list:
            for path in files:
                exp(path)
        else :
            exp(files)

    def save_dat(self, data, file_path, sep='\t', index=True, na_rep=-9999):
        """
        Fonction pour sauver les inputs modifiés sans les écrire directement
        INPUT:
            data [df] : the data to save
            file_path[str] : the path to the .dat file relative to inpute file
            index[bool] : if True, the index will be in the file
            na_rep[str,int,...] : The value by witch the NaN should be replaced in the file (Some are "NULL", -9999,...)
        OUTPUT:
            None
        """
        if file_path in [a[1] for a in self._lst_changes_inputs]:
            self.delete_dat(file_path)
        self._lst_changes_inputs.append([data, file_path, sep, index, na_rep])
        return None

    def delete_dat(self, file_path):
        """
        Fonction pour supprimer les inputs modifiés et revenir a la normal
        INPUT:
            file_path[str] : the path to the .dat file relative to inpute file
        OUTPUT:
            None
        """
        self._lst_changes_inputs = [i for i in self._lst_changes_inputs if i[1] != file_path]
        return None

    def _write_all_dat(self):
        """
        Fonction pour ecrire tous les inputs modifiées
        INPUT:
            None
        OUTPUT:
            None
        """
        for element in self._lst_changes_inputs:
            self._write_dat(data = element[0],file_path= element[1], sep= element[2], index= element[3],na_rep= element[4])

    def _write_dat(self,data,file_path, sep='\t', index=True,na_rep="",keepDef=True):
        """
        Fonction pour écrire les séries temporelles en .dat
        INPUT :
            data [df] : the data to save
            file_path[str] : the path to the .dat file relative to inpute file
            index[bool] : if True, the index will be in the file
            keepDef[bool]: if True, the data will be written in a safe way (keep a default value as {file_path}-def)
        OUTPUT :
            None
        """
        if keepDef==True:
            if file_path in self._get_Mod_input_files():
                None
            else :
                self._Add_Mod_input_files(file_path)
        data.to_csv(f"{self.pathToSirane}/{self.inputs.FICH_DIR_INPUT}/{file_path}",sep=sep,index=index, date_format='%d/%m/%Y %H:%M', na_rep=na_rep)


class Output(object):
    """
    Une classe qui traite les sorties de modèle de qualité de l'air pour plusieurs stations de mesure.
    INPUT:
        Model : instance de classe de modèle PySiWrap
        AddHourOfDay [bool = False] : ajout de l'heure du jours, le jour du mois et le mois comme variable.
        Relative_Residual [bool = True] : Calcule les résidus en relatif : (Obs-Mod)/(Obs+Mod)*2
        //Filter []: un filtre pour elever des données des calcules

    ATTRIBUTS:
        mod : instance d'une classe de modèle de qualité de l'air
        alldata : un DataFrame qui contient les données de réception pour toutes les stations de mesure
        list_of_Recept : une liste de noms de stations de mesure qui ont des données de réception
        list_of_species : une liste des noms d'espèces pour lesquelles les données sont disponibles dans alldata
        list_of_VarMeteo : une liste des noms de variables météorologiques

    MÉTHODES:
        scatterplots : Affiche des graphiques de dispersion pour les espèces et les récepteurs spécifiés.
        qQplots : Affiche les graphiques quantile-quantile pour chaque combinaison d'espèces et de récepteurs.
        residualplots : Affiche les graphiques des résidus en boîtes selon différentes variables météorologiques pour chaque espèce.
        indicateurs : Calcule les indcateurs 'FB', 'MG', 'NMSE', 'VG', 'R' et 'FAC2' des données pour differents recepteur et especes.
    """
    def __init__(self, Model, AddHourOfDay = False, Relative_Residual=True, Filter = None, extended_meteo=False):
        """
        Initialise une instance de la classe Output.
        """
        super(Output, self).__init__()
        self.mod = Model
        if self.mod.inputs.FICH_RECEPT in self.mod.get_Mod_changes_inputs():
            ls = list([i[0] for i in self.mod.get_Mod_changes_inputs(data=True) if i[1]==self.mod.inputs.FICH_RECEPT][0].index)
        else:
            ls = list(read_dat(self.mod.inputs.FICH_RECEPT, mod_input=self.mod).index)
        self.list_of_Recept = [i for i in ls if isfile(f"{self.mod.pathToSirane}/{self.mod.inputs.FICH_DIR_RESUL}/RECEPT/Recept_{i}.dat")]

        data = [read_dat(f"RECEPT/Recept_{i}.dat", mod_result=self.mod) for i in self.list_of_Recept]
        self.alldata = pd.concat(data, keys=self.list_of_Recept).reset_index(level=0)
        self.alldata = self.alldata.rename(columns = {'level_0':'Recept'})
        self.list_of_species=[sp for sp in set([i.split("_")[0] for i in self.alldata.columns if (i not in ["Recept", "Date"])]) if len(self.alldata.loc[:,[f"{sp}_Mod", f"{sp}_Mes"]].dropna()) != 0]
        self.list_of_species= np.sort(self.list_of_species)
        for species in self.list_of_species:
            if Relative_Residual:
                self.alldata[f"{species}_Res"] = (self.alldata[f"{species}_Mod"] / self.alldata[f"{species}_Mes"]) #/ (self.alldata[f"{species}_Mes"])
            else:
                self.alldata[f"{species}_Res"] = self.alldata[f"{species}_Mod"] - self.alldata[f"{species}_Mes"]
        self.list_of_VarMeteo = None

        self.AddHourOfDay = AddHourOfDay
        if AddHourOfDay:
            self.alldata["Heure"] = self.alldata.index.hour
            self.alldata["Jour_semaine"] = self.alldata.index.dayofweek
            self.alldata["Jour_mois"] = self.alldata.index.day
            self.alldata["Mois"] = self.alldata.index.month

        meteo = read_dat(self.mod.inputs.FICH_METEO,mod_input=self.mod)
        meteo = meteo.loc[self.alldata.index,:]
        if extended_meteo == True:
            meteo2 = read_dat("METEO/Resul_Meteo.dat",mod_result=self.mod).loc[self.alldata.index,["Hcla","Ustar","SigmaTheta","Cld","H0","Lmo","Thetastar","k1","k3"]]
            meteo = pd.concat([meteo, meteo2], axis=1)
        self.list_of_VarMeteo = list(meteo.columns)
        if self.AddHourOfDay: #ajout de l'heure comme variable "météo"
            self.list_of_VarMeteo += ["Heure","Jour_semaine","Jour_mois","Mois"]
        self.alldata = pd.concat([self.alldata, meteo], axis=1)

    def __repr__(self):
        return f"PySiWrap.Output for model named '{self.mod.name}'"

    def scatterplots(self, list_of_Recept=None, list_of_species=None):
        """
        Affiche des graphiques de dispersion pour les espèces et les récepteurs spécifiés.
        INPUT:
            list_of_Recept [list]: Une liste des identifiants de récepteur pour lesquels les graphiques de dispersion doivent être affichés. Si cette liste n'est pas fournie, tous les récepteurs seront inclus.
            list_of_species [list]: Une liste des espèces pour lesquelles les graphiques de dispersion doivent être affichés. Si cette liste n'est pas fournie, toutes les espèces seront incluses.
        OUTPUT:
            fig [matplotlib.figure.Figure]: Figure contenant les graphiques de dispersion.
            axs [numpy.ndarray]: Tableau des axes des sous-graphiques pour chaque espèce.
        """
        if list_of_Recept is None:
            list_of_Recept = self.list_of_Recept
        if list_of_species is None:
            list_of_species = self.list_of_species
        if len(list_of_Recept) <= 10:
            sns.set_palette(sns.color_palette("tab10"), len(list_of_Recept))
        else:
            sns.set_palette(sns.color_palette("tab20"), len(list_of_Recept))
        if len(list_of_species)<=0 :
            return (None,None)
        fig, axs = plt.subplots(1,len(list_of_species))
        if len(list_of_Recept) > 1:
            legend = True
        else:
            legend = False

        filter_recept = self.alldata.Recept.apply(lambda x : x in list_of_Recept)
        for id,species in enumerate(list_of_species):
            min = self.alldata.loc[filter_recept,[f"{species}_Mod", f"{species}_Mes"]].min().min()
            max = self.alldata.loc[filter_recept,[f"{species}_Mod", f"{species}_Mes"]].max().max()
            min = min - (max-min)/20
            max = max + (max-min)/20
            axs[id].plot([min,max],[min,max], ls="--", c=".5")
            axs[id] = sns.scatterplot(data= self.alldata[filter_recept], y=f"{species}_Mod", x=f"{species}_Mes", hue="Recept", ax=axs[id], legend=legend, s=10)
            axs[id].set_aspect("equal")
            axs[id].set(xlim=(min, max), ylim=(min, max))
            axs[id].set_ylabel('Concentrations modélisées',fontsize=15)
            axs[id].set_xlabel('Concentrations mesurées',fontsize=15)
            axs[id].set_title(f"{species}",fontsize=15)
            if (id == 0) and (len(list_of_Recept) > 1):
                legend = False
                axs[id].legend(title="Recepteurs")
        return fig, axs

    def qQplots(self, list_of_Recept=None, list_of_species=None):
        """
        Affiche les graphiques quantile-quantile pour chaque combinaison d'espèces et de récepteurs.

        INPUT :
            list_of_Recept (list[str]) : Liste des récepteurs à inclure dans le graphique. Si la valeur est None, toutes les récepteurs seront affichées.
            list_of_species (list[str]) : Liste des noms des espèces à afficher. Si la valeur est None, toutes les espèces seront affichées.
        OUTPUT :
            fig (matplotlib.figure.Figure) : L'objet Figure contenant le tracé.
            axs (numpy.ndarray) : Un tableau d'objets AxesSubplot représentant les sous-graphes de la figure.
        """
        if list_of_Recept is None:
            list_of_Recept = self.list_of_Recept
        if list_of_species is None:
            list_of_species = self.list_of_species
        quantile = lambda series: np.quantile(series.dropna(), np.arange(0,1,1/len(series)))# - series.min())/(series.max()-series.min())
        if len(list_of_Recept)*len(list_of_species) == 0:
            return (None,None)
        fig, axs = plt.subplots(len(list_of_Recept),len(list_of_species))
        if type(axs) != type(np.ndarray([])):
            d = np.ndarray(shape=(2,2),dtype=object )
            d[0,0] = axs
            axs = d
        elif len(axs.shape) == 1:
            if len(list_of_Recept) == 1 :
                d = np.ndarray(shape=(2,axs.shape[0]),dtype=object)
                d[0,:] = axs
                axs = d
            else:
                d = np.ndarray(shape=(axs.shape[0],2),dtype=object)
                d[:,0] = axs
                axs = d

        for n_el, species in enumerate(list_of_species):
            for n_rec, recept in enumerate(list_of_Recept):
                if len(self.alldata.loc[self.alldata["Recept"] == recept,[f"{species}_Mes",f"{species}_Mod"]].dropna()) == 0:
                    fig.delaxes(axs[n_rec, n_el])
                else:
                    observed_data = quantile(self.alldata.loc[self.alldata["Recept"] == recept,f"{species}_Mes"])
                    predicted_data = quantile(self.alldata.loc[self.alldata["Recept"] == recept,f"{species}_Mod"])
                    min = np.array((observed_data.min(),predicted_data.min())).min()
                    max = np.array((observed_data.max(),predicted_data.max())).max()
                    min = min - (max-min)/20
                    max = max + (max-min)/20
                    axs[n_rec, n_el].scatter(observed_data, predicted_data, color='blue')
                    axs[n_rec, n_el].plot([min,max],[min,max], ls="--", c=".5")
                    axs[n_rec, n_el].set(xlim=(min, max), ylim=(min, max))
                    axs[n_rec, n_el].set_ylabel('Concentrations modélisées',fontsize=15)
                    axs[n_rec, n_el].set_xlabel('Concentrations mesurées',fontsize=15)
                    axs[n_rec, n_el].set_title(f"{species}",fontsize=15)
                    axs[n_rec, n_el].set_aspect("equal")

        #fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        #fig.set_size_inches(10, 10)
        #fig.subplots_adjust(hspace=0)
        #print(axs)
        return fig, axs[:len(list_of_Recept),:len(list_of_species)]

    def residualplots(self, list_of_VarMeteo = None, list_of_species=None, list_of_Recept=None, max_bins=35):
        """
        Affiche les graphiques des résidus en boîtes selon différentes variables météorologiques pour chaque espèce.
        INPUT:
        - list_of_VarMeteo (list[str]) : Liste des noms des variables météorologiques à afficher. Si la valeur est None, toutes les variables météorologiques seront affichées.
        - list_of_species (list[str]) : Liste des noms des espèces à afficher. Si la valeur est None, toutes les espèces seront affichées.
        - list_of_Recept (list[str]) : Liste des noms des recpeteurs à prendre en compte. Si la valeur est None, toutes les recepteurs seront pris en compte.
        - max_bins (int) : Nombre de bins maximum par plot
        OUTPUT:
        - fig (matplotlib.figure.Figure) : Figure matplotlib contenant les graphiques de résidus.
        - axs (numpy.ndarray) : Tableau numpy contenant les axes de chaque graphique de résidus.
        """
        if list_of_VarMeteo is None:
            list_of_VarMeteo = self.list_of_VarMeteo
        if list_of_species is None:
            list_of_species = self.list_of_species
        if list_of_Recept is None:
            list_of_Recept = self.list_of_Recept

        filter_recept = self.alldata.Recept.apply(lambda x : x in list_of_Recept)
        fig, axs = plt.subplots(len(list_of_VarMeteo),len(list_of_species))
        fig.set_size_inches(12.4 * len(list_of_species) /2 ,17.5 * len(list_of_VarMeteo) / 5)
        if type(axs) != type(np.ndarray([])):
            d = np.ndarray(shape=(2,2),dtype=object)
            d[0,0] = axs
            axs = d
        elif len(axs.shape) == 1:
            if len(list_of_VarMeteo) == 1 :
                d = np.ndarray(shape=(2,axs.shape[0]),dtype=object)
                d[0,:] = axs
                axs = d
            else:
                d = np.ndarray(shape=(axs.shape[0],2),dtype=object)
                d[:,0] = axs
                axs = d
        def nbr_bins(data):
            nbr_diff = len(set(data))
            if (nbr_diff>0)&(nbr_diff < max_bins):
                return nbr_diff
            else:
                return int(max_bins)
        #sns.set_palette(sns.color_palette("tab20"))
        for n_var, variable in enumerate(list_of_VarMeteo):
            nbr_xvalue = len(set(self.alldata[filter_recept][variable]))
            if (nbr_xvalue>0)&(nbr_xvalue <= max_bins):
                self.alldata[f"{variable}_Group"] = self.alldata[filter_recept][variable]
                retbins = np.sort(np.array(list(set(self.alldata[filter_recept][variable]))))
                xlables = mticker.AutoLocator().tick_values(retbins.min(),retbins.max())
                xticks  = (xlables-retbins[1])/(retbins[-1]-retbins[1])*len(retbins[2:])+1
            else:
                self.alldata[f"{variable}_Group"],retbins = pd.cut(self.alldata[filter_recept][variable], nbr_bins(self.alldata[filter_recept][variable]), retbins=True)
                xlables = mticker.AutoLocator().tick_values(retbins.min(),retbins.max())
                xticks = (xlables-retbins[1])/(retbins[-1]-retbins[1])*len(retbins[2:])+0.5
            for n_el, species in enumerate(list_of_species):
                if len(self.alldata.loc[filter_recept,[f"{species}_Res"]].dropna()) == 0:
                    fig.delaxes(axs[n_rec, n_el])
                else:
                    axs[n_var, n_el] = sns.boxplot(data = self.alldata[filter_recept], x=f"{variable}_Group", y=f"{species}_Res", orient="v", ax = axs[n_var, n_el], color = "palegreen")
                    axs[n_var, n_el].set_xticks(xticks)
                    axs[n_var, n_el].set_xticklabels(xlables)
                    axs[n_var, n_el].set_xlabel(f"{variable}", fontsize=15)
                    axs[n_var, n_el].set_ylabel(f"Résidus de {species}", fontsize=15)

        return fig,axs[:len(list_of_VarMeteo),:len(list_of_species)]

    def indicateurs(self, list_of_Recept=None, list_of_species=None, yearlyMQI=False):
        """
        Calcule les indcateurs 'FB', 'MG', 'NMSE', 'VG', 'R' et 'FAC2' des données pour differents recepteur et especes.
        INPUT :
            list_of_Recept [list] : la liste des recepteur ou il faut calculer les indicateurs.
            list_of_species [list]: la liste des especes pour lesquelles il faut calculer les indicateurs.
            yearlyMQI [bool] : si vrais, les MQI annuelles de chaque espece et indicateur sera calculé
        OUTPUT :
            [dict] un dictionnaire contenant par especes, un df avec les indicateurs par recepteurs
        """
        if list_of_Recept is None:
            list_of_Recept = self.list_of_Recept
        if list_of_species is None:
            list_of_species = self.list_of_species
        output = {}
        for esp in list_of_species:
            output[esp] = pd.DataFrame([calculate_indicators(self.alldata[self.alldata.Recept == rec], f"{esp}_Mod", f"{esp}_Mes", esp, yearlyMQI=yearlyMQI) for rec in list_of_Recept], index = list_of_Recept)
        return output

def transph_grd(file_path_grd, file_path_tif = None, mod_input = None, mod_result=None):
    """
    Fonction pour lire et transformer les grd file en .tif
    INPUT :
        file_path_grd [str] : The path to the grd file (can't be rotated)
        file_path_tif [str] : the path where the tif file should be writen. If set to None -> no file conversion.
    OUTPUT :
        a tuple of values [np.array] and extend [lst] of the grid file.
    """
    if mod_input is not None:
        file_path_grd = f"{mod_input.pathToSirane}/{mod_input.inputs.FICH_DIR_INPUT}/{file_path_grd}"
    if mod_result is not None:
        file_path_grd = f"{mod_result.pathToSirane}/{mod_result.inputs.FICH_DIR_RESUL}/{file_path_grd}"

    t = gdal.Open(file_path_grd)
    value = t.ReadAsArray()
    GT = t.GetGeoTransform()
    extent = [GT[0], GT[0] + value.shape[1] * GT[1], GT[3] + value.shape[0] * GT[5],  GT[3]]
    if file_path_tif != None:
        t.SetProjection("EPSG:31370") #To set the Projection to Lambert 72
        t = gdal.Translate(file_path_tif, t)
    return (value, extent, t)

def force_all_files_to_def(model):
    """
    Restore tout les fichiers modifié dans une autre instance (finisant par "-def") par les fichier par défaut
    INPUT:
        model [PySiWrap.Model] : le modèle aillant des "-def" dans les fichiers INPUTS
    OUTPUT:
        liste de tous les fichiers remis à leurs valeur par défaut.
    """
    Paths = [i[:-4].split("INPUT\\")[-1] for i in glob.glob(f"{model.pathToSirane}\\{model.inputs.FICH_DIR_INPUT}/**/*-def", recursive=True)]
    [model._Del_Mod_input_file(PathToFile=i, force_to_def = True) for i in Paths]
    return Paths

#indcateurs stats:
def calculate_indicators(df, C_p, C_o, esp, LOQ = None, yearlyMQI = False):
    """
    Calcule les 'FB', 'MG', 'NMSE', 'VG', 'R', 'FAC2', 'NAD_t', 'NAD' des paires de données.
    INPUT:
        df[pd.DataFrame] : Un DataFrame avec les données observée et prédite. Les paires de données sont sur la même ligne
        C_p[str] : le nom de la colone aillant les données prédite par le models
        C_o[str] : le nom de la colone aillant les données observées
        esp[str] : Le nom de l'espece concerné (peut être [NO2, O3, PM10, PM2.5], si autre, alors MQI est np.nan)
        LOQ[float] : La limite of quantification des detecteurs (si None, NAD_t est np.nan)
        yearlyMQI[bool] : Si vrais, le MQI annuel est aussi donnée
    OUTPUT:
        [dict] : un dictionaire aillant les valeurs des indicateurs. Pour les MQI, si plus de 25% des paires de données sont manquantes, les valeurs sont rendue négative.
    """
    df = df[[C_p,C_o]].copy()
    df.loc[df.isna().sum(axis=1)>0,:] = np.nan
    #df = df.dropna()
    # calculer les moyennes et écart-types
    C_p_mean = df[C_p].mean()
    C_o_mean = df[C_o].mean()
    C_p_std = df[C_p].std()
    C_o_std = df[C_o].std()
    # calculer FB
    FB = (C_o_mean - C_p_mean) / (0.5 * (C_o_mean + C_p_mean))
    # calculer MG
    MG = np.exp(np.log(df[C_o]).mean() - np.log(df[C_p]).mean())
    # calculer NMSE
    NMSE = ((df[C_o] - df[C_p])**2).mean() / (C_o_mean * C_p_mean)
    # calculer VG
    VG = np.exp(((np.log(df[C_o]) - np.log(df[C_p]))**2).mean())
    # calculer R
    R = (((df[C_o] - C_o_mean) * (df[C_p] - C_p_mean)).mean()) / (C_o_std * C_p_std)
    # calculer FAC2
    FAC2 = ((df[C_p] / df[C_o] >= 0.5) & (df[C_p] / df[C_o] <= 2.0)).mean()
    # Calculer NAD_t (treshold)
    if LOQ is not None:
        CT = 3 * LOQ
        false_positives = ((df[C_p] > CT) & (df[C_o] < CT)).sum()
        false_negatives = ((df[C_p] < CT) & (df[C_o] > CT)).sum()
        AF = (false_positives + false_negatives) / len(df)
        overlap = ((df[C_p] > CT) & (df[C_o] > CT)).sum()
        AOV = overlap / len(df)
        NAD_t = AF / (AF + AOV)
    else :
        NAD_t = np.nan
    #Calculer NAD sans treshold
    NAD = (df[C_o] - df[C_p]).abs().mean()/(C_p_mean + C_o_mean)

    #Calculer le MQI
    def U(Oi, esp):
        """
        Calcule l'incertitude de mesure
        """
        params = pd.DataFrame({"NO2" : {"Ur" : 0.24, "RV": 200, "a" : 0.2, "Np":5.2, "Nnp":5.5},
                               "O3"  : {"Ur" : 0.18, "RV": 120, "a" : 0.79, "Np":11, "Nnp":3},
                               "PM10": {"Ur" : 0.28, "RV": 50, "a" : 0.25, "Np":20, "Nnp":1.5},
                               "PM2.5":{"Ur" : 0.36, "RV": 25, "a" : 0.5, "Np":20, "Nnp":1.5}})
        return params[esp]["Ur"] * ((1-params[esp]["a"]**2) * Oi**2 + params[esp]["a"]**2 * params[esp]["RV"]**2)**0.5

    if esp in ["NO2", "PM10", "PM2.5"]:
        MQI = ((df[C_o] - df[C_p])**2).sum()**0.5 / 2 / (U(df[C_o], esp)**2).sum()**0.5
        if (df[[C_o,C_p]].isna().sum(axis=1)>=1).sum()/len(df) > 0.25:
            MQI = - MQI
    elif esp == "O3":
        #maximum journalier des moyennes glissantes sur huit heures
        dfO3 = df[[C_o,C_p]].rolling(8, min_periods=6).mean()
        dfO3["date"] = dfO3.index.date
        dfO3_day = dfO3.groupby("date").max()
        dfO3_day.loc[dfO3.groupby("date").count().sum(axis=1)>=18,:] = np.nan
        MQI = (np.sum(dfO3_day[C_o] - dfO3_day[C_p])**2)**0.5 / 2 / (U(dfO3_day[C_o], esp)**2).sum()**0.5
        if dfO3[C_o].isna().sum()/len(dfO3) > 0.25:
            MQI = - MQI
    else :
        MQI = np.nan

    #Calculer le MQI Annuel
    def U_year(Oi,esp):
        """
        Calcule l'incertitude de mesure
        """
        params = pd.DataFrame({"NO2" : {"Ur" : 0.24, "RV": 200, "a" : 0.2, "Np":5.2, "Nnp":5.5},
                               "O3"  : {"Ur" : 0.18, "RV": 120, "a" : 0.79, "Np":11, "Nnp":3},
                               "PM10": {"Ur" : 0.28, "RV": 50, "a" : 0.25, "Np":20, "Nnp":1.5},
                               "PM2.5":{"Ur" : 0.36, "RV": 25, "a" : 0.5, "Np":20, "Nnp":1.5}})
        return params[esp]["Ur"] * ((1-params[esp]["a"]**2) * Oi**2 / params[esp]["Np"] + params[esp]["a"]**2 * params[esp]["RV"]**2 /params[esp]["Nnp"])**0.5
    if esp in ["NO2", "O3", "PM10", "PM2.5"] and yearlyMQI:
        MQI_y = np.abs(C_p_mean - C_o_mean)/2/U_year(C_o_mean,esp)
        if (df[[C_o,C_p]].isna().sum(axis=1)>=1).sum()/len(df) > 0.25:
            MQI_y = - MQI_y
    else :
        MQI_y = np.nan

    return {'FB': FB, 'MG': MG, 'NMSE': NMSE, 'VG': VG, 'R': R, 'FAC2': FAC2, 'NAD_t': NAD_t, "NAD" : NAD, "MQI": MQI, "MQI_y":MQI_y}
