"""
Created on Wed May 20 08:11:42 2020

@author: fcuevas
"""
import numpy as np
import pandas as pd

# librerías de Bokeh
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource,TableColumn,DataTable,TextInput,Select,NumberFormatter,Range1d,HoverTool,MultiChoice,CustomJS,TapTool,StringFormatter,LegendItem, CustomJS, Legend
from bokeh.events import Tap
from bokeh.io import curdoc, show
from bokeh.layouts import column, row, Spacer
from bokeh.models.widgets import Button

from bokeh.palettes import Category20

import geopandas as gpd
import fiona

from shapely.geometry import Point, LineString


fiona.drvsupport.supported_drivers['KML'] = 'rw'

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#--------------------------------------------------- PATH -------------------------------------------------------------
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
import sys 
import os
from os.path import basename, normpath, join, dirname

sys.path
p = lambda x: x if basename(normpath(x))=='sun4heat' else p(os.path.dirname(x))
sun4heat_path = p(os.getcwd())
script_path = os.path.join(sun4heat_path,'scripts')
sys.path.append(script_path)

Path_e=lambda x: x if basename(normpath(x))=='visualizaciones' else Path_e(os.path.dirname(x))
path=Path_e(os.getcwd())
path = os.path.join(path,'')

datos_path= os.path.join(sun4heat_path,'datos')
path_to_folder= path+'2022_mapa_emisiones/'


tiles = ["STAMEN_TERRAIN","ESRI_IMAGERY",'OSM']

archivo = 'Indus_II_2022.csv'

from funciones_bokeh_plots import create_table

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#--------------------------------------------------- FUNCIONES -------------------------------------------------------------
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

# leer archivo indus_II_2021 o "emisiones-fuentes-puntuales-2021" y crear indus_II_2021
def ReadIndus():
    

    """
    Función que lee la base de emisiones entregada por el ministerio (variable por trimestre) 
    y entrega un DF para manipular.
    Parameters
    ----------
        None

    Returns
    -------
    indus: Data Frame

    Headers DataFrame 
    -----------------
    'AÑO', 'CODIGO VU', 'NOMBRE ESTABLECIMIENTO', 'RAZON SOCIAL',
           'RUT RAZON SOCIAL', 'CODIGO CIIU6', 'CIIU6', 'RUBRO', 'REGION',
           'PROVINCIA', 'COMUNA', 'COORDENADA NORTE', 'COORDENADA ESTE',
           'CODIGO DE FUENTE', 'TIPO DE FUENTE', 'COMBUSTIBLE PRIMARIO',
           'CCF8 PRIMARIO', 'EMISION PRIMARIO', 'COMBUSTIBLE SECUNDARIO',
           'CCF8 SECUNDARIO', 'EMISION SECUNDARIO', 'CCF8 MATERIA PRIMA',
           'EMISION MATERIA PRIMA', 'CONTAMINANTE', 'TOTAL EMISION', 'ORIGEN'
           
    Es posible que sea necesario cambiar el nombre de las columnas
    
    cols_2022(['Unnamed: 0', 'año', 'razon_social', 'rut_razon_social',
           'nombre_establecimiento', 'id_vu', 'ciiu6', 'rubro_vu', 'region',
           'provincia', 'comuna', 'id_comuna', 'latitud', 'longitud',
           'cantidad_toneladas', 'unidad', 'contaminantes', 'id_contaminantes',
           'fuente_emisora_general', 'tipo_fuente', 'id_fuente_emisora',
           'combustible_primario', 'emision_primario', 'combustible_secundario',
           'emision_secundario', 'emision_materia_prima', 'sistema'],
          dtype='object')

    
    cols_2021 (['AÑO', 'CODIGO VU', 'NOMBRE ESTABLECIMIENTO', 'RAZON SOCIAL',
           'RUT RAZON SOCIAL', 'CODIGO CIIU6', 'CIIU6', 'RUBRO', 'REGION',
           'PROVINCIA', 'COMUNA', 'COORDENADA NORTE', 'COORDENADA ESTE',
           'CODIGO DE FUENTE', 'TIPO DE FUENTE', 'COMBUSTIBLE PRIMARIO',
           'CCF8 PRIMARIO', 'EMISION PRIMARIO', 'COMBUSTIBLE SECUNDARIO',
           'CCF8 SECUNDARIO', 'EMISION SECUNDARIO', 'CCF8 MATERIA PRIMA',
           'EMISION MATERIA PRIMA', 'CONTAMINANTE', ' TOTAL EMISION', 'ORIGEN'],
          dtype='object')
    
    cols={"RUBRO_VU":"RUBRO", "ID_VU":"CODIGO VU", "TIPO_FUENTE":"TIPO DE FUENTE",
          "CONTAMINANTES":"CONTAMINANTE", "SISTEMA":"ORIGEN"}

    """
    try:
        indus = pd.read_csv(
            datos_path + "/RETC/"+archivo, encoding="utf-8-sig")
        del indus[indus.columns[0]]
        
        
        
    except:
        ####################################################################################################
        ######################################## LECTURA EMISIONES #########################################
        
        
        indus = pd.read_csv(datos_path + "/RETC/emisiones-fuentes-puntuales-2022.csv",
                            decimal=",", encoding="utf-8-sig")

        indus.columns = indus.columns.str.strip() ##Elimina los espacios laterales de las columnas
        
        indus.rename(str.upper, axis='columns', inplace=True)
        
        indus['TOTAL EMISION'] = pd.to_numeric(indus['CANTIDAD_TONELADAS'], errors="coerce")
        
        # Se eliminan del dataset aquellas fuentes que no reporten emisiones
        indus = indus[indus['TOTAL EMISION'] >= 0.1]

# =============================================================================
#         indus.rename(columns={'COORDENADA NORTE': 'LATITUD',
#                      'COORDENADA ESTE': 'LONGITUD'}, inplace=True)
# =============================================================================
        indus.rename(columns={"RUBRO_VU":"RUBRO", "ID_VU":"CODIGO VU", "TIPO_FUENTE":"TIPO DE FUENTE",
              "CONTAMINANTES":"CONTAMINANTE", "SISTEMA":"ORIGEN"}, inplace=True)
        
        indus['LATITUD'] = pd.to_numeric(indus['LATITUD'], errors="coerce")
        indus['LONGITUD'] = pd.to_numeric(indus['LONGITUD'], errors="coerce")
        
        
        indus.drop(["AÑO", "EMISION_SECUNDARIO", "ORIGEN",
                    "EMISION_PRIMARIO","EMISION_MATERIA_PRIMA"],
                   axis=1, inplace=True)
        
        indus.dropna(
            subset=(["TOTAL EMISION", "LONGITUD", "LATITUD"]), inplace=True)

        indus['EQUIPO'] = indus['ID_FUENTE_EMISORA'].str[:2]
        indus['CALNUM'] = indus['EQUIPO'].apply(lambda x: 1 if x in ['CF', 'CA', 'CG', 'CR', 'IN'] else 0)
        
        indus = indus[(indus.LATITUD>=-1000)&(indus.LONGITUD>-1000)&(indus.LATITUD<1000)&(indus.LONGITUD<1000)]  ### si hay localidades con latitud y longitud del orden de los -300000


        ##############################################################################################################
        ######################################## LECTURA EMISIONES CCF8 ##############################################

        base = pd.read_excel(datos_path + "/RETC/fe_ccf8.xlsx").iloc[:,1:]
        
        ####################################################################################################
        ################################### Emission to energy #############################################

        XX,YY=[], []
       
        def fe_tmp(row):
            x = row['EQUIPO']
            y = row['COMBUSTIBLE_PRIMARIO']
            
            peso_comb=list(base.fe_cff8[(base.equipo == x)&(base.combustible == y)])

            if peso_comb==[]:
                if x in ['CF', 'CA', 'CG', 'CR', 'IN']:
                    XX.append(x)
                    YY.append(y)
                return -1
            else:
                return peso_comb[0]
        
        if len(XX)!=0:
            df_fe = pd.DataFrame()
            df_fe['EQUIPO']=XX
            df_fe['COMBUSTIBLE_PRIMARIO']=YY
            df_fe.drop_duplicates(inplace=True)
            
            print('')
            print('WARNING!')
            print('Modificar archivo {}'.format(datos_path + "/RETC/fe_ccf8.xlsx"))
            print('Guardando archivo con datos faltantes:...')
            df_fe.to_excel(path_to_folder+'Warnings_info_faltante.xlsx')
            print('')


        indus['fe_tmp'] = indus.apply(fe_tmp, axis=1)
        indus['peso_combustible'] = (indus['TOTAL EMISION']*10**3)/indus.fe_tmp

        temp = []
        for i, j in zip(indus.peso_combustible, indus['COMBUSTIBLE_PRIMARIO']):
            if i < 0:
                temp.append(0)
            else:
                PCI=list(base['PCI(KWh/Kg)'][base.combustible==j])[0]
                temp.append((PCI*i)*10**(-3))

        indus['ener_cons_CO2'] = temp

        indus.columns = indus.columns.str.replace(' ', '_')
        indus.to_csv(datos_path + "/RETC/"+archivo)

    return indus

def Read_GeoData():
    # empty GeoDataFrame
    df = gpd.GeoDataFrame()
    
    # iterate over layers
    string=''
    L_tipo=[]
    
    for layer in fiona.listlayers(path_to_folder+"SISTEMA ELÉCTRICO NACIONAL.kml"):
        try:
            string_aux=string+'-'+layer if layer!='Subestaciones' else ' -Subestacion'
            s = gpd.read_file(path_to_folder+"SISTEMA ELÉCTRICO NACIONAL.kml", driver='KML', layer=layer)
            df = pd.concat([df, s], ignore_index=True)
            
            ltipo=[string_aux]*len(s.index)
            L_tipo+=ltipo
            
        except:
            string=layer
    
    df['Tipo']=L_tipo
    df['Categoria']=df['Tipo'].apply(lambda x: x.split('-')[0])
    df['Clase']=df['Tipo'].apply(lambda x: x.split('-')[1])
    
    df['Nombre']=df['Name'].str.extract(r"([^_][A-Z]*[/]?[A-Z]?(\s[A-Z]*).*)").iloc[:,0]
    df['Categoria']=df['Categoria'].replace(['Centrales', 'Compensadores Activos'], ['Central','Compensador Activo'])
    
    df['Clase']=df['Clase'].replace(['Eólicos','Solares Fotovoltaicos', 'Geotermia','Subestaciones'], ['Eólica','Solar Fotovoltaico', 'Geotermal','Subestacion'])
    df['Clase']=df['Clase'].replace(r'Hidroeléctricas', 'Hidroeléctrica', regex=True)
    df['Clase']=df['Clase'].replace(r'Termoeléctricas', 'Termoeléctrica', regex=True)
    df['Clase']=df['Clase'].replace(r'Almacenamientos', 'Almacenamiento', regex=True)
    
    df['Clase']=df['Clase'].apply(lambda x: 'Termoeléctrica' if 'Termoeléctrica' in x else x)
    
    gama_colores={
        'Eólica': 'blue',
        'Geotermal':'brown',
        'Hidroeléctrica de Embalse':'purple',
        'Hidroeléctrica de Pasada':'purple', 
        'Solar Fotovoltaico':'yellow',
        'Termoeléctrica':'orange', 
        'Termosolar de Concentración':'gray',
        'Almacenamiento de Energía':'green', 
        '500 kV':'#084594', 
        '345 kV':'#2171b5', 
        '220 kV':'#4292c6',
        '154 kV':'#6baed6', 
        '110 kV':'#9ecae1', 
        '66 kV':'#c6dbef', 
        '44 ': '#deebf7', 
        'Subestacion':'green',
        'Conexión en Derivación':'#f7fbff',
        }
    
    df['Color']=df.Clase.apply(lambda x: gama_colores[x])
    
    
    def getLineCoords(row, geom, coord_type):
        """Returns a list of coordinates ('x' or 'y') of a LineString geometry"""
        if coord_type == 'x':
            return list( row[geom].coords.xy[0] )
        elif coord_type == 'y':
            return list( row[geom].coords.xy[1] )
        
    def getPointCoords(row, geom, coord_type):
        """Calculates coordinates ('x' or 'y') of a Point geometry"""
        if coord_type == 'x':
            return row[geom].x
        elif coord_type == 'y':
            return row[geom].y   
    
    
    trazados= df[df['geometry'].apply(lambda x: isinstance(x, LineString))]
    puntos= df[df['geometry'].apply(lambda x: isinstance(x, Point))]
    
    trazados=trazados.copy()
    trazados['x'] = trazados.apply(getLineCoords, geom='geometry', coord_type='x', axis=1)
    trazados['y'] = trazados.apply(getLineCoords, geom='geometry', coord_type='y', axis=1)
    
    puntos=puntos.copy()
    puntos['x'] = puntos.apply(getPointCoords, geom='geometry', coord_type='x', axis=1)
    puntos['y'] = puntos.apply(getPointCoords, geom='geometry', coord_type='y', axis=1)
    
    trazados = trazados.drop('geometry', axis=1).copy()
    puntos = puntos.drop('geometry', axis=1).copy()
    
    k = 6378137
    
    trazados['x']=trazados['x'].apply(lambda x: np.array(x)*(k * np.pi / 180.0))
    trazados['y'] = trazados['y'].apply(lambda y: np.log(np.tan((90 + np.array(y))* np.pi/360.0))* k)
    
    puntos['x']=puntos['x'].apply(lambda x: x * (k * np.pi / 180.0))
    puntos['y']=puntos['y'].apply(lambda y: np.log(np.tan((90 + y) * np.pi / 360.0)) * k)
    
    return puntos, trazados

# filtrar según el equipo o mercado a analizar
def FiltEquip(df, mkt):
    """
    Función que filtra los equipos del mercado según la DF entregada (emisiones_aire_año_cart.csv).

    Equipos / Mercado
    -----------------
    ACTUALIZAR EQUIPOS!!!
        Caldera Calefacción ('CA')
        Caldera Industrial ('IN')
        Mercado Solar ('IN', 'CF', 'CA')
        Mercado H2 ('IN','CF','CA','PC','PS')
        Generación Eléctrica ('GE')
        Todo ('CA', 'IN', 'PC', 'CF', 'PS', 'GE')



    Parameters
    ----------
    df : DataFrame
        Corresponde al DF de emisiones.
    mkt : String
        Corresponde al equipo o mercado a analizar.

    Returns
    -------
    df : Data Frame
        DF con los filtros aplicados.

    """

    nl = []

    if "Todo" in mkt:
        pass
    else:

        if "Antorcha" in mkt:
            nl.append("AN")

        if "Caldera de Fluido Térmico" in mkt:
            nl.append("CF")

        if "Caldera de Generación Eléctrica" in mkt:
            nl.append("CG")

        if "Calentador" in mkt:
            nl.append("CL")

        if "Caldera Recuperadora" in mkt:
            nl.append("CR")

        if "Convertidor Teniente (CT)" in mkt:
            nl.append("CV")

        if "Convertidor Pierce Smith (CPS)" in mkt:
            nl.append("CV")

        if "Caldera Calefacción (CA)" in mkt:
            nl.append("CA")

        if "Grupo Electrógeno" in mkt:
            nl.append("EL")

        if "Horno" in mkt:
            nl.append("HR")

        if "Incinerador" in mkt:
            nl.append("IC", "MO")

        if "Molino de Rodillo" in mkt:
            nl.append("MC")

        if "Marmita de Calcinación" in mkt:
            nl.append("MC", "MO")

        if "Caldera Industrial (IN)" in mkt:
            nl.append("IN")

        if "Motor Generación Eléctrica" in mkt:
            nl.append("MG")

        if "Turbina de Gas" in mkt:
            nl.append("TG")

        if "Regenerador Cracking Catalítico (FCCU)" in mkt:
            nl.append("RG")

        if "Secadores" in mkt:
            nl.append("SC")

        if "Mercado Solar" in mkt:
            nl.append("IN", "CF", "CA")

        if "Mercado H2" in mkt:
            nl.append("IN", "CF", "CA", "PC", "PS")

        df = df[df.EQUIPO.isin(nl)]

    return df

def IndusFilt(df, min_ton, max_ton):
    """
    Función que filtra el DF entregado (emisiones de industrias) según un rango de toneladas (min_ton, max_ton). 
    También agrupa según el ID, sumando las columnas de ton_emision y n_equip, agregando columna de max_emision (conjunto de empresas).

    Parameters
    ----------
    df : DataFrame
        DF de industrias.
    min_ton : Float
        Cantidad mínima de toneladas de emisiones a filtrar.
    max_ton : Float
        Cantidad máxima de toneladas de emisiones a filtrar.

    Returns
    -------
    df : DataFrame
        Nueva agrupación de datos por ID y sumando ton_emision y n_equip


    """
    df=df.copy()
    df["n_equip"] = 1
    #df = df.sort_values("TOTAL_EMISION", ascending=False).drop_duplicates("TIPO_DE_FUENTE")
    df = df.sort_values("TOTAL_EMISION", ascending=False)

    indus_gr = df.groupby(["CODIGO_VU"]).agg(  # group by ID?
        {
            "TOTAL_EMISION": "sum",
            'ener_cons_CO2': "sum",
            "n_equip": "sum",
            "CALNUM": "sum",
            "RAZON_SOCIAL": "first",
            "RUT_RAZON_SOCIAL": "first",
            "NOMBRE_ESTABLECIMIENTO": "first",
            "RUBRO": "first",
            "CIIU6": "first",  # no hay ciiu4, hay ciiu6
            "REGION": "first",
            "PROVINCIA": "first",
            "COMUNA": "first",
            "COMBUSTIBLE_PRIMARIO": "first",
            "COMBUSTIBLE_SECUNDARIO": "first",
            "TIPO_DE_FUENTE": "first",
            "CONTAMINANTE": "first",
            "LATITUD": "first",
            "LONGITUD": "first"
        }
    )

    indus_max = df.groupby(["CODIGO_VU"])["TOTAL_EMISION"].max()
    indus_gr["max_emision"] = indus_max

    df = indus_gr
    df["orden"] = 1

    if max_ton > min_ton:
        df = df[(df["TOTAL_EMISION"] > min_ton) &
                (df["TOTAL_EMISION"] <= max_ton)]

    else:
        df = df[df["TOTAL_EMISION"] > min_ton]

    return df

# filtrar por categorias
def Filtrbr(df, rbr, max_empr):
    """
    Función que filtra según las categorias presentes en 'emisiones_aire_año_cart.csv'.

    Categorias
    ----------
            'Otras actividades',
             'Comercio minorista',
             'Captación, tratamiento y distribución de agua',
             'Otras industrias manufactureras',
             'Pesca y acuicultura',
             'Plantas de tratamiento de aguas servidas',
             'Comercio mayorista',
             'Producción agropecuaria',
             'Ventas y mantención de vehículos automotores',
             'Construcción',
             'Minería',
             'Termoeléctricas',
             'Otras centrales de generación eléctrica',
             'Industria del papel y celulosa',
             'Industria química, de plástico y caucho',
             'Industria de la madera y silvicultura',
             'Refinería de petróleo',
             'Gestores de residuos',
             'Producción de cemento, cal y yeso'

    Parameters
    ----------
    df : DataFrame
        DF en donde se filtra la categoría correspondiente.
    rbr : List
        Lista de strings que contiene las categorías a filtrar.
    max_empr : int
        Por ahora ninguna función.

    Returns
    -------
    df : DataFrame
        DF con la categoría filtrada.

    """

    if rbr == ["Todo"]:
        rbr = list(indus.RUBRO.unique())

    else:
        pass

    df = df[df.RUBRO.isin(rbr)]

    return df

# filtrar por region
def FiltRegion(df, rgn, latNor, latSur):
    """
    Función que filtra por región o según rango de latitud.

    Parameters
    ----------
    df : DataFrame
        Base de datos 'emisiones_aire_año_cart.csv'.
    rgn : string
        Region a filtrar.
    latNor : float
        Latitud norte.
    latSur : float
        Latitud sur.

    Returns
    -------
    df_filt : DataFrame.
        DF con el filtro de región/latitud realizados.       
    """

    if "Rango latitud" in rgn:
        df_filt = df[(df.LATITUD < latNor) & (df.LATITUD > latSur)]
    elif "Todas" in rgn:
        df_filt = df
    else:
        df_filt = df[df.REGION.isin(rgn)]

    return df_filt

def TableResumen(df):
    tot_emis = df['TOTAL_EMISION'].sum()/1000
    tot_ener_cons = df.ener_cons_CO2.sum()
    tot_empresas = len(df)

    Data = [[tot_emis, tot_ener_cons, tot_empresas, df.CALNUM.sum()]]

    table_res = pd.DataFrame(Data, columns=['Miles de Toneladas de emisión anual (miles Ton/año)',
                                            'Energía consumida anual (MWh/año)', 'N° Fuentes puntuales', 'N° de calderas'])

    return table_res

def wgs84_to_web_mercator(df, lon="LONGITUD", lat="LATITUD"):
    """
    Función que convierte las escalas de longitud en una variable 'x' y latitud
    en una variable 'y', generando nuevas columnas 'x' 'y' en el DF entregado.

    Parameters
    ----------
    df : DataFrame
        DF en el que se generará las variables 'x', 'y' ('emisiones_aire_año_cart.csv').
    lon : Column
        Pertenece a la columna "Longitud" del DF.
    lat : Column
        Pertenece a la columna "Latitud" del DF.

    Returns
    -------
    df : DataFrame
        DF con las nuevas columnas 'x', 'y' generadas.

    """
    k = 6378137
    df["x"] = df[lon] * (k * np.pi / 180.0)
    df["y"] = np.log(np.tan((90 + df[lat]) * np.pi / 360.0)) * k

    return df

def Inside(x,y,lista,radio):
    for punto in lista:
        if x<punto[0]+radio:
            if x>punto[0]-radio:
                if y<punto[1]+radio:
                    if y>punto[1]-radio:
                        return punto[2]
                else:
                    continue
            else:
                continue
        else:
            continue
    return 'Nan'

def Closeness(df_1,lista, radio):
    df=df_1.copy()
    df['close']=df[['x','y']].apply(lambda z: Inside(z['x'],z['y'],lista,radio), axis=1)
    return df
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#--------------------------------------------------- BOTONES Y MENUS -------------------------------------------------------------
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
wdt = 345

# crear dataframe (df) indus
indus = ReadIndus()

puntos, trazados=Read_GeoData()

# paleta de colores para los gráficos
clr = dict(zip(indus.RUBRO.unique(), Category20[20]))

# crear lista de contaminantes
ctms_opt = ["Todo"] + list(indus.CONTAMINANTE.unique())
dropDownCtms = Select(value="Dióxido de carbono (CO2)", title="Contaminante", options=ctms_opt, width=wdt)

# crear lista de combustibles
comb_list = ['Todo']+list(indus['COMBUSTIBLE_PRIMARIO'].unique())
dropDownComb = MultiChoice(value=["Todo"], title="Combustible Primario", options=comb_list, width=700)

minTon = TextInput(value=str(100), title="Mínimo emisiones anuales [ton/año]", width=wdt)
maxTon = TextInput(value=str(1e6), title="Máximo emisiones anuales [ton/año]", width=wdt)

mrc = ["Todo",
    "Antorcha",
    "Caldera de Fluido Térmico",
    "Caldera de Generación Eléctrica",
    "Calentador",
    "Caldera Recuperadora",
    "Convertidor Teniente (CT)",
    "Convertidor Pierce Smith (CPS)",
    "Caldera Calefacción (CA)",
    "Grupo Electrógeno",
    "Horno de Panadería",
    "Incinerador",
    "Molino de Rodillo",
    "Marmita de Calcinación",
    "Caldera Industrial (IN)",
    "Motor Generación Eléctrica",
    "Turbina de Gas",
    "Regenerador Cracking Catalítico (FCCU)",
    "Secadores",
    "Mercado Solar",
    "Mercado H2",
    "Generación eléctrica"]

Equipos_por_defecto=["Caldera Calefacción (CA)",
                    "Caldera Industrial (IN)",
                    "Caldera de Fluido Térmico",
                    "Calentador"]



dropdownEquip = MultiChoice(value=Equipos_por_defecto, title="Equipo térmico", options=mrc, width=700)


rubro = ["Todo"] + list(indus.RUBRO.unique())

Rubros_por_defecto = ['Captación, tratamiento y distribución de agua',
                      'Industria del papel y celulosa',
                      'Industria química, de plástico y caucho',
                      'Industria de la madera y silvicultura',
                      'Fundiciones de cobre',
                      'Pesca y acuicultura',
                      'Gestores de residuos',
                      'Otras industrias manufactureras',
                      'Minería',
                      'Plantas de tratamiento de aguas servidas',
                      'Producción agropecuaria',
                      'Producción de cemento, cal y yeso']

rbr_multi_choice = MultiChoice(title="Rubro", value=Rubros_por_defecto, options=rubro, width=700)

region = ["Todas","Rango latitud"]+list(indus.REGION.unique())

dropdownRegion = MultiChoice(value=["Todas"], title="Region o Rango de latitud", options=region, width=wdt)

latNorte = TextInput(value=str(-18.4), title="Latitud norte (Opción rango latitud)", width=wdt)
latSur = TextInput(value=str(-35), title="Latitud sur (Opción rango latitud)", width=wdt)
maxEmpr = TextInput(value=str(1000), title="Total empresas", width=wdt)

# boton para filtrar y resetear mapa.
buttCalcUpdate = Button(label="Filtrar", button_type="success", width=wdt)
radio_closeness = TextInput(value=str(1000), title="Radio cercania a S/E (Km)", width=wdt)
dropDownTiles = Select(value="ESRI_IMAGERY", title="Tipo mapa", options=tiles)

indus_ft = FiltEquip(indus, dropdownEquip.value)    
indus_ft = IndusFilt(indus_ft, 100, 1e6) #Cruza la base agrupada con la categoría de actividad
indus_ft = wgs84_to_web_mercator(indus_ft)
indus_ft = Filtrbr(indus_ft, rbr_multi_choice.value, int(maxEmpr.value))
indus_ft = FiltRegion(indus_ft, dropdownRegion.value,float(latNorte.value),float(latSur.value))
indus_ft["pt_size"] = np.log(indus_ft['TOTAL_EMISION']+10)
indus_ft["clr"] = indus_ft.RUBRO.map(clr)
indus["f_ind"] = indus['TIPO_DE_FUENTE']
indus = indus.set_index("f_ind") # Definir nuevo ID por fuente de emisión



lista_xy_sub=np.array(puntos[['x','y','Name']][puntos.Clase=='Subestacion'])
radius=float(radio_closeness.value)*1000

df_aux=Closeness(indus_ft,lista_xy_sub, radius)
df_aux=df_aux[df_aux.close!='Nan']

df_by_region=df_aux.copy().groupby(['REGION']).sum(numeric_only=True)
df_by_region=df_by_region[['TOTAL_EMISION','ener_cons_CO2','CALNUM']]
df_by_region['ener_cons_CO2']=df_by_region['ener_cons_CO2']/1000
df_by_region['TOTAL_EMISION']=df_by_region['TOTAL_EMISION']/1000
                                            
df_by_region=df_by_region.rename(columns={'TOTAL_EMISION':'Emisión anual (miles Ton/año)',
                             'ener_cons_CO2':'Energía consumida anual (GWh/año)',
                             'CALNUM':'N° de calderas'})

df_by_region['N° Fuentes puntuales']=list(df_aux.groupby(['REGION']).count().iloc[:,0])

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#--------------------------------------------------- PATH -------------------------------------------------------------
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

# definir titulo de columnas de una tabla

columns_dict={"NOMBRE_ESTABLECIMIENTO":'Nombre',
              "RAZON_SOCIAL":'Razon social',
              "RUT_RAZON_SOCIAL":"Rut Razon social",
              "TOTAL_EMISION":"Emisiones (ton/año)",
              "ener_cons_CO2":"Energía consumida anual (MWh/año)",
              "REGION":"Región",
              "COMBUSTIBLE_PRIMARIO":"Combustible Primario",
              "RUBRO":"Rubro",
              "CIIU6":"CIIU6",
              "close":"S/E Cerca",
              'LATITUD': 'Latitud'}

# iniciar tabla con columnas y fuente de datos ds source_indus
data_table = create_table(df_aux[columns_dict.keys()],columns_dict=columns_dict,form='nro',height_table=300,width_1c=200, width_c=130)
data_tableres = create_table(TableResumen(df_aux),height_table=80,form='nro',width_1c=220, width_c=140)

DataTable_by_Region=create_table(df_by_region.reset_index(), width_1c=80, width_c=150)
DataTable_by_Region.columns[-2].width=100
DataTable_by_Region.columns[-3].width=100
DataTable_by_Region.columns[-4].width=190


# ########################################################################################

#######################################################################################
# iniciar mapa
#tile_provider = get_provider('ESRI_IMAGERY')
p1 = figure(
    width=1200,
    height=1000,
    tools=["pan,wheel_zoom,box_zoom,reset,save"],
    x_axis_type="mercator",
    y_axis_type="mercator",
    x_range=(-9000000, -6000000),
    y_range=(-6000000, -1200000),
)
p1.add_tile('ESRI_IMAGERY')


rsource = ColumnDataSource(puntos[puntos.Name.isin(df_aux.close)])
Radio=p1.circle(x='x',y='y',source=rsource,radius=radius+50, line_color='Color', legend_label='Radio {} km'.format(radius),radius_units='data', line_width=5, fill_color=None)
Radio.visible=False

for r in df_aux.RUBRO.unique():
    source_indus_aux = ColumnDataSource(df_aux[df_aux['RUBRO']==r])
    sct=p1.scatter(x='x',y='y',fill_color='clr',fill_alpha=0.8, legend_label=r,source=source_indus_aux, size='pt_size')
    p1.add_tools(HoverTool(
        renderers=[sct],
        tooltips=[
            ("Nombre: ", "@NOMBRE_ESTABLECIMIENTO"),
            ("Region: ", "@REGION"),
            ("Emisiones (ton/año): ", "@TOTAL_EMISION{0.00}"),
            ("Energía consumida MWh/año): ", "@ener_cons_CO2{0.00}"),
            ("Rubro: ", "@RUBRO"),
            ("Combustible Primario: ", "@COMBUSTIBLE_PRIMARIO"),
        ],
    ))
    
for est in puntos.Clase.unique():
    if est=='Subestacion':
        psource = ColumnDataSource(puntos[puntos.Name.isin(df_aux.close)])
    else:
        psource = ColumnDataSource(puntos[puntos['Clase']==est])
        
    Estacion=p1.scatter(x='x',y='y',fill_color='Color', legend_label=est,source=psource, size=10)
    p1.add_tools(HoverTool(
        tooltips="""
        <div>
            <div>
                <span style="font-size: 18px; font-weight: bold;">@Nombre</span>
                
            </div>
            <div>
                <span style="font-size: 16px; color: green;">@Categoria @Clase</span>
            </div>
        </div>
        """, renderers=[Estacion]))

for traz in trazados.Clase.unique():
    msource = ColumnDataSource(trazados[trazados['Clase']==traz])
    Trazado=p1.multi_line(xs="x", ys="y",line_color='Color',legend_label=traz,source=msource, line_width=5)
    p1.add_tools(HoverTool(
        tooltips="""
        <div>
            <div>
                <span style="font-size: 18px; font-weight: bold;">@Name</span>
        """, renderers=[Trazado]))
        

p1.legend.click_policy = "hide"
taptool = p1.select(type=TapTool)
p1.on_event(Tap)

########################################################################

# iniciar tabla específica de empresa
empr1 = indus_ft["NOMBRE_ESTABLECIMIENTO"].iloc[0]
df_empr = indus_ft[(indus_ft['NOMBRE_ESTABLECIMIENTO'] == empr1)]

source_empr = ColumnDataSource(data=df_empr)

columns2_dict={"NOMBRE_ESTABLECIMIENTO":'Nombre',
              "TIPO_DE_FUENTE":'Fuente emisión',
              "TOTAL_EMISION":"Emisiones (ton/año)",
              "ener_cons_CO2":"Energía consumida anual (MWh/año)",
              "CONTAMINANTE":"Contaminante",
              "COMBUSTIBLE_PRIMARIO":"Combustible Primario",
              "COMBUSTIBLE_SECUNDARIO":"Combustible Secundario"}

data_tableEmpr = create_table(df_empr[columns2_dict.keys()],columns_dict=columns2_dict,form='nro',height_table=200,width_1c=200, width_c=150)


#############################################################################################
# definir coordenadas del mapa específico de una empresa
lat = df_empr.x
lon = df_empr.y

offSet = 5000
ymin = lon.iloc[0] - offSet
ymax = lon.iloc[0] + offSet

yrng = Range1d()  # ?
yrng.start = ymin
yrng.end = ymax

xmin = lat.iloc[0] - offSet
xmax = lat.iloc[0] + offSet
xrng = Range1d()  # ?
xrng.start = xmin
xrng.end = xmax

p = figure(
    width=700,
    height=700,
    tools=["pan,wheel_zoom,box_zoom,reset,save"],
    x_axis_type="mercator",
    y_axis_type="mercator",
    x_range=xrng,
    y_range=yrng,
)
p.add_tile('ESRI_IMAGERY')

source = ColumnDataSource(data=dict(lat=lat, lon=lon))

p.circle(x="lat", y="lon", size=10, fill_color="blue",
         fill_alpha=0.8, source=source)
###############################################################################################

###################
# crear funcion para cambiar mapa de empresa específica (cambio al clickear empresa en tabla superior)


def function_source(attr, old, new):
    """
    Función que permite cambiar el mapa de la empresa específica (la cual cambia al clickear empresa en la tabla superior)

    Parameters
    ----------
    attr : TYPE
        DESCRIPTION.
    old : TYPE
        DESCRIPTION.
    new : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    try:
        selected_index = data_table.source.selected.indices[0]
        name_selected = data_table.source.data["NOMBRE_ESTABLECIMIENTO"][selected_index]
        
        if dropDownCtms.value== 'Todo':
            df_empr =indus[indus['NOMBRE_ESTABLECIMIENTO'] == name_selected]
        else:
            df_empr =indus[(indus['NOMBRE_ESTABLECIMIENTO'] == name_selected)&
                           (indus['CONTAMINANTE'] == dropDownCtms.value)]
            
        df_empr = wgs84_to_web_mercator(df_empr)
        data_tableEmpr.source.data = df_empr
        table_res = TableResumen(df_empr)
        data_tableres.source.data = table_res

        lat = df_empr.x
        lon = df_empr.y
        new_data = dict(lat=lat, lon=lon)
        source.data = new_data

        ymin = lon.iloc[0] - offSet
        ymax = lon.iloc[0] + offSet

        xmin = lat.iloc[0] - offSet
        xmax = lat.iloc[0] + offSet
        xrng.update(start=xmin, end=xmax)
        yrng.update(start=ymin, end=ymax)
        ##############

    except IndexError:
        pass


# crear funcion para cambiar mapa y tabla general
data_table.source.selected.on_change("indices", function_source)


def UpdateTable():
    """
    Función creada para un boton. Permite leer y procesar el archivo de 'emisiones_aire_año_cart.csv', 
    aplicando todos los filtros colocados.

    Returns
    -------
    None.

    """
    
    ind_aux=indus.copy()
    
    ctm = dropDownCtms.value        #Contaminante + Todo
    comb = dropDownComb.value       #Combustible primario + Todo
    mkt = dropdownEquip.value       #Equipo + Todo
    min_ton = float(minTon.value)   #Minimo de emisiones
    max_ton = float(maxTon.value)   #Maximo de emisiones
    rbr = rbr_multi_choice.value    #Rubro + Todo
    rn = dropdownRegion.value       #Region + Todo
    latN = float(latNorte.value)    #Latitud sur
    latS = float(latSur.value)      #Latitud norte
    max_empr = int(maxEmpr.value)   #Maximo de empresas
    
    
    if ctm == "Todo":
        ind_aux=indus
        
    else:
        ind_aux = ind_aux[ind_aux.CONTAMINANTE == ctm]

    if comb == ["Todo"]:
        pass
    else:
        ind_aux = ind_aux[ind_aux.COMBUSTIBLE_PRIMARIO.isin(comb)]


    indus_tmp_aux = FiltEquip(ind_aux, mkt)    
    indus_ft_aux = IndusFilt(indus_tmp_aux, min_ton, max_ton)
    indus_ft_aux = wgs84_to_web_mercator(indus_ft_aux)
    indus_ft_aux = Filtrbr(indus_ft_aux, rbr, max_empr)
    indus_ft_aux = FiltRegion(indus_ft_aux, rn, latN, latS)
    
    data_tableEmpr.source.data = indus_ft_aux

    indus_ft_aux = wgs84_to_web_mercator(indus_ft_aux)
    pt_size = np.log(indus_ft_aux['TOTAL_EMISION']+10)
    indus_ft_aux["pt_size"] = pt_size
    indus_ft_aux["clr"] = indus_ft_aux.RUBRO.map(clr)
    
    lista_xy_sub=np.array(puntos[['x','y','Name']][puntos.Clase=='Subestacion'])
    radius=float(radio_closeness.value)*1000
    
    df_aux=Closeness(indus_ft_aux,lista_xy_sub, radius)
    df_aux=df_aux[df_aux.close!='Nan']
    

    table_res = TableResumen(indus_ft_aux)

    data_tableres.source.data = table_res
    data_tableEmpr.source.data = indus_ft_aux
    data_table.source.data = df_aux
        
    p1.renderers=[]
    p1.legend.items = []
    p1.tools=p1.tools[:5]
    
    
    rsource = ColumnDataSource(puntos[puntos.Name.isin(df_aux.close)])
    Radio=p1.circle(x='x',y='y',source=rsource,radius=radius+50, line_color='Color', radius_units='data',legend_label='Radio {} Km'.format(radio_closeness.value), line_width=5, fill_color=None)
    Radio.visible=False
    
    for r in df_aux.RUBRO.unique():
        source_indus_aux = ColumnDataSource(df_aux[df_aux['RUBRO']==r])
        sct=p1.scatter(x='x',y='y',fill_color='clr',fill_alpha=0.8, legend_label=r,source=source_indus_aux, size='pt_size')
        p1.add_tools(HoverTool(
            renderers=[sct],
            tooltips=[
                ("Nombre: ", "@NOMBRE_ESTABLECIMIENTO"),
                ("Region: ", "@REGION"),
                ("Emisiones (ton/año): ", "@TOTAL_EMISION{0.00}"),
                ("Energía consumida MWh/año): ", "@ener_cons_CO2{0.00}"),
                ("Rubro: ", "@RUBRO"),
                ("Combustible Primario: ", "@COMBUSTIBLE_PRIMARIO"),
            ],
        ))
        
    for est in puntos.Clase.unique():
        if est=='Subestacion':
            psource = ColumnDataSource(puntos[puntos.Name.isin(df_aux.close)])
        else:
            psource = ColumnDataSource(puntos[puntos['Clase']==est])
            
        Estacion=p1.scatter(x='x',y='y',fill_color='Color', legend_label=est,source=psource, size=10)
        p1.add_tools(HoverTool(
            tooltips="""
            <div>
                <div>
                    <span style="font-size: 18px; font-weight: bold;">@Nombre</span>
                    
                </div>
                <div>
                    <span style="font-size: 16px; color: green;">@Categoria @Clase</span>
                </div>
            </div>
            """, renderers=[Estacion]))

    for traz in trazados.Clase.unique():
        msource = ColumnDataSource(trazados[trazados['Clase']==traz])
        Trazado=p1.multi_line(xs="x", ys="y",line_color='Color',legend_label=traz,source=msource, line_width=3)
        p1.add_tools(HoverTool(
            tooltips="""
            <div>
                <div>
                    <span style="font-size: 18px; font-weight: bold;">@Name</span>
            """, renderers=[Trazado]))
    

    p1.add_tile(dropDownTiles.value)
    
    df_by_region=df_aux.copy().groupby(['REGION']).sum(numeric_only=True)
    df_by_region=df_by_region[['TOTAL_EMISION','ener_cons_CO2','CALNUM']]
    df_by_region['ener_cons_CO2']=df_by_region['ener_cons_CO2']/1000
    df_by_region['TOTAL_EMISION']=df_by_region['TOTAL_EMISION']/1000
                                                
    df_by_region=df_by_region.rename(columns={'TOTAL_EMISION':'Emisión anual (miles Ton/año)',
                                 'ener_cons_CO2':'Energía consumida anual (GWh/año)',
                                 'CALNUM':'N° de calderas'})

    df_by_region['N° Fuentes puntuales']=list(df_aux.groupby(['REGION']).count().iloc[:,0])

    DataTable_by_Region.source.data=df_by_region
    
    
    
    
    

def DownloadButton():
    """
    Función que filtra y elimina las columnas no deseables para la descarga del df

    Returns
    -------
    nw : Data Source
        DS similar al source_indus pero con menos columnas (de no interés)
    """
    nw=data_table.source.data[[
        "NOMBRE_ESTABLECIMIENTO",
        "RAZON_SOCIAL",
        "RUBRO",
        "CIIU6",
        "REGION",
        "PROVINCIA",
        "TIPO_DE_FUENTE",
        "COMBUSTIBLE_PRIMARIO",
        "COMBUSTIBLE_SECUNDARIO",
        "CONTAMINANTE",
        "TOTAL_EMISION",
    ]]
    return nw

#############################################
buttCalcUpdate.on_click(UpdateTable)
button = Button(label="Descargar", button_type="success")
button.js_on_click(CustomJS(args=dict(source=data_table.source),
                            code=open(join(dirname(__file__), "download.js")).read()))
#############################################


spc = 50
layout = column(
    row(dropDownCtms, minTon, maxTon, maxEmpr),
    Spacer(height=spc),
    row(dropdownEquip, dropDownComb),
    row(rbr_multi_choice),
    Spacer(height=spc),
    row(dropdownRegion),
    row(latNorte, latSur),
    row(dropDownTiles),
    row(radio_closeness),
    row(buttCalcUpdate, button),
    # row(b),
    Spacer(height=spc - 20),
    row(p1),
    row(data_table),
    row(p,data_tableEmpr),
    row(data_tableres),
    row(DataTable_by_Region),
    # Spacer(height= 10),
)

#show(layout)
curdoc().add_root(layout)

#%%

output_file(filename="mapa_emisiones.html", title="Zeus Tech")
save(layout)

#%%

# =============================================================================
# 
# guia_met = pd.read_csv(datos_path + "/RETC/datos_guia_met_fe1.csv", sep=",").iloc[:,1:]
# 
# def agregar(guia_met, combustible, fe):
#     l=len(guia_met)
#     guia_met.loc[l]=['CA',combustible,fe]
#     guia_met.loc[l+1]=['CF',combustible,fe]
#     guia_met.loc[l+2]=['IN',combustible,fe]
#     guia_met.loc[l+3]=['CG',combustible,fe]
#     
#     return guia_met
# 
# def agregar(guia_met, combustible, fe):
#     l=len(guia_met)
#     guia_met.loc[l]=['CA',combustible,fe]
#     guia_met.loc[l+1]=['CF',combustible,fe]
#     guia_met.loc[l+2]=['IN',combustible,fe]
#     guia_met.loc[l+3]=['CG',combustible,fe]
#     
#     return guia_met
# 
# 
# 
# 
# 
# =============================================================================

