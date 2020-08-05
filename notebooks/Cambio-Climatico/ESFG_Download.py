# _Autor:_    __Salavador Navas__
# _Revisión:_ __05/08/2020__


import os
import xarray as xr
from pyesgf.logon import LogonManager
from pyesgf.search import SearchConnection
import tqdm
import numpy as np
from netCDF4 import Dataset
from math import *

def download_ESGF_data(Open_ID, password, server, project, experiment,time_frequency, variable, domain, path_output):
    """Esta función nos permite descargar masivamente mediante WGET los diferentes ficheros netcdf que contienen los servidrores de ESGF sobre cambio climático.
    
    El servidor por defecto que se va a utilizar es: https://esgf-data.dkrz.de/projects/esgf-dkrz/.
    
    Parámetros:
    ---------------------
    Open_ID         : string. ID de tu usuario para acceder a la base de datos correspondiente del servidor
    password        : string. Contraseña correspondiente a la ID
    server          : string. Servidor del que se desea descargar la información. Ejemplo: https://esgf-data.dkrz.de/esg-search
    project         : string. Proyecto dentro del servidor del que se quiere descargar los datos. Ejemplo: CORDEX, CMIP5, CMIP6
    experiment      : string. Escenarios de cambio climático. Ejemplo: historical, rcp26, rcp45, rcp85
    time_frequency  : string. Frecuencia de la base de datos que se quiere. Ejemplo: 1hr, 6hr, day, mon
    variable        : string. Variable que se desea descargar: tasmax, tasmin, pr 
    domain          : string. En el caso de que se desee descargar CORDEX, se debe de incluir el nombre de la malla. Ejemplo: EUR-11
    path_output     : string. Directorio donde se desean guardar los ficheros
    
    Salidas:
    ----------------------
    Ficheros netcdf para cada uno de los escenarios y modelos solicitados
    
    """
    conn = SearchConnection('https://esgf-data.dkrz.de/esg-search', distrib=True)
    lm = LogonManager()
    lm.logoff()
    lm.is_logged_on()
    lm.logon_with_openid(Open_ID, password)
    lm.is_logged_on()
    if project=='CORDEX':
        ctx = conn.new_context(
        project=project,
        experiment=experiment,
        time_frequency=time_frequency,
        variable=variable,
        domain=domain,)
    else:
        ctx = conn.new_context(
        project=project,
        experiment=experiment,
        time_frequency=time_frequency,
        variable=variable,)
        
    with open('wget-plantilla-ESGF.sh', "r+") as out_file:
        lines = out_file.readlines()
        
    lines_first=lines[:27]
    lines_end=lines[28:]
    
    for ct in tqdm.tqdm(range(ctx.hit_count)):
        files_list=list()
        result = ctx.search()[ct]
        lines[23]="earch_url=https://esgf-data.dkrz.de/esg-search/wget/?distrib=false&dataset_id="+result.dataset_id+'\n'

        files = result.file_context().search()
        ntcf_ds=list()
        ntcf_name=list()
        for file in files:
            if variable in file.opendap_url:
                files_list.append("'"+file.filename+"'"+' '+"'"+file.download_url+"'"+' '+"'"+file.checksum_type+"'"+' '+"'"+file.checksum+"'"+'\n')
        with open(path_output+"Download.sh", "w") as fh:
            for line in (lines_first+files_list+lines_end):
                fh.write(line)
                
        conn = SearchConnection('https://esgf-data.dkrz.de/esg-search', distrib=True)
        lm = LogonManager()
        lm.logoff()
        lm.is_logged_on()

        lm.logon_with_openid(Open_ID, password)
        lm.is_logged_on()
        os.chdir(path_output)
        os.system('bash '+'Download.sh'+' H '+Open_ID+' '+ password)   
        
        


def rotated_grid_transform(grid_in, option, SP_coor):
    """Permite tranformar coordenadas rotadas al sistema de corrdenadas WGS84 o viceversa.
    
    Parámetros:
    ---------------------
    grid_in : array. Coordendas que se desean tranformar
    option  : int. 1: # WGS84 -> Rotadas  2: # Rotadas -> WGS84 
    SP_coor : list. grid_north_pole #SP_lon = NP_lon - 180, SP_lat = -NP_lat. Los parámetros NP_ se obtienen del netcdf o fichero de datos.
    
    Salidas:
    ----------------------
    lon_new , lat_new : Coordendas tansformadas

    """
    
    lon = grid_in[0]
    lat = grid_in[1];

    lon = (lon*pi)/180; # Convert degrees to radians
    lat = (lat*pi)/180;

    SP_lon = SP_coor[0];
    SP_lat = SP_coor[1];
    
    #SP_lon = NP_lon - 180, SP_lat = -NP_lat.

    theta = 90+SP_lat; # Rotation around y-axis
    phi = SP_lon; # Rotation around z-axis

    theta = (theta*pi)/180;
    phi = (phi*pi)/180; # Convert degrees to radians

    x = cos(lon)*cos(lat); # Convert from spherical to cartesian coordinates
    y = sin(lon)*cos(lat);
    z = sin(lat);

    if option == 1: # Regular -> Rotated

        x_new = cos(theta)*cos(phi)*x + cos(theta)*sin(phi)*y + sin(theta)*z;
        y_new = -sin(phi)*x + cos(phi)*y;
        z_new = -sin(theta)*cos(phi)*x - sin(theta)*sin(phi)*y + cos(theta)*z;

    else:  # Rotated -> Regular

        phi = -phi;
        theta = -theta;

        x_new = cos(theta)*cos(phi)*x + sin(phi)*y + sin(theta)*cos(phi)*z;
        y_new = -cos(theta)*sin(phi)*x + cos(phi)*y - sin(theta)*sin(phi)*z;
        z_new = -sin(theta)*x + cos(theta)*z;



    lon_new = atan2(y_new,x_new); # Convert cartesian back to spherical coordinates
    lat_new = asin(z_new);

    lon_new = (lon_new*180)/pi; # Convert radians back to degrees
    lat_new = (lat_new*180)/pi;

    return lon_new , lat_new

def extract_CORDEX_EUR11(path_output,area=False,lon_min_area=None,lat_min_area=None,lon_max_area=None,lat_max_area=None,point=True,
                         lon_point=None,lat_point=None,name_point=None):
    """ Esta función permite descargar de G los datos de cambio climático CORDEX de la malla EUR-11 en el punto o area que introduzcamos. La propia función va a elegir el puntos más cercano.
        Además se puede extraer un area seleccionada.
        
        Parámetros:
        --------------------- 
        path_output  : string. Path donde se desean guardar los resultados
        area         : True o False. Opción para saber que tipo de resultados se quiere resultados en un area
        lon_min_area : float. Longitud mínima del area si area está identificada como True
        lat_min_area : float. Latitud mínima del area si area está identificada como True
        lon_max_area : float. Longitud máxima del area si area está identificada como True
        lat_max_area : float. Latitud máxima del area si area está identificada como True
        point        : True o False. Opción para saber que tipo de resultados se quiere resultados en un punto
        lon_point    : float. Longitud del punto 
        lat_point    : float. Latitud del punto
        name_point   : string. Nombre que se quiere dar al punto
        
        Salidas:
        ----------------------
        Coordenadas en sistema WGS84: Fichero CSV para cada uno de los modelos
        Variables para cada modelo: Fichero CSV con los datos extraidos el area o punto introducidos
        
        """
    
    variables=['pr','tasmax','tasmin']
    Experiment=['historical','rcp45','rcp85']
    modelos=[['CNRM-CERFACS-CNRM-CM5','SMHI-RCA4','r1i1p1'],
         ['CNRM-CERFACS-CNRM-CM5','KNMI-RACMO22E','r1i1p1'],
         ['ICHEC-EC-EARTH','SMHI-RCA4','r12i1p1'],
         ['ICHEC-EC-EARTH','KNMI-RACMO22E','r1i1p1'],
         ['ICHEC-EC-EARTH','KNMI-RACMO22E','r12i1p1'],
         ['ICHEC-EC-EARTH','DMI-HIRHAM5','r3i1p1'],
         ['ICHEC-EC-EARTH','CLMcom-CCLM4-8-17','r12i1p1'],
         ['IPSL-IPSL-CM5A-MR','SMHI-RCA4','r1i1p1'],
         ['IPSL-IPSL-CM5A-MR','IPSL-WRF381P','r1i1p1'],
         ['MOHC-HadGEM2-ES','SMHI-RCA4','r1i1p1'],
         ['MOHC-HadGEM2-ES','DMI-HIRHAM5','r1i1p1'],
         ['MOHC-HadGEM2-ES','CLMcom-CCLM4-8-17','r1i1p1'],
         ['MPI-M-MPI-ESM-LR','SMHI-RCA4','r1i1p1'],
         ['MPI-M-MPI-ESM-LR','CLMcom-CCLM4-8-17','r1i1p1'],
         ['NCC-NorESM1-M','DMI-HIRHAM5','r1i1p1']]
    itera=0
    for i in tqdm.tqdm(modelos):
        for ex in Experiment:
            for v in variables:
                direc, root =find('G:/CLIMA/02_HYDRO-CLIMATE/CORDEX/',i[0],i[1],i[2],ex,v)
                if area==True:
                    [lon_min_,lat_min_]= rotated_grid_transform((lon_min_area,lat_min_area), 1, (18,-39.25))
                    [lon_max_,lat_max_]= rotated_grid_transform((lon_max_area,lat_max_area), 1, (18,-39.25))
                    
                    if os.path.exists(path_output+'/'+v+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+
                           direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'.csv')==True:
                        continue
                    else:
                        print('### Descargando: '+v+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'.csv')
                        da=xr.open_mfdataset(root[0]+'/'+'*.nc')
                        da = da.sel(rlat=slice(lat_min_-1.5, lat_max_+1.5), rlon=slice(lon_min_-0.5, lon_max_+0.5))
                        [XX,YY]=np.meshgrid(da['rlon'].values,da['rlat'].values)
                        [lon_new,lat_new]=rotated_grid_transform((XX,YY), 2, (18,-39.25))
                        try:
                            time=da.indexes['time'].to_datetimeindex()
                        except:
                            time=da.indexes['time']
                        if v=='pr':
                            csv=pd.DataFrame((da[v].values*86400).flatten().reshape(len(time),np.size(XX)),
                                             index=time)
                        else:
                            csv=pd.DataFrame((da[v].values-273).flatten().reshape(len(time),np.size(XX)),
                                             index=time)
                        Coordenadas=pd.DataFrame(index=np.arange(0,len(lon_new.flatten())),columns=['Lon','Lat'])
                        Coordenadas.iloc[:,0]=lon_new.flatten()
                        Coordenadas.iloc[:,1]=lat_new.flatten()

                        csv.to_csv(path_output+v+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+
                                   direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'.csv')
                        Coordenadas.to_csv(path_output+'Coordenadas'+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+
                                   direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'.csv')
                        if itera==0:
                            fig, ax = plt.subplots(figsize=(14, 10),subplot_kw=dict(projection=ccrs.PlateCarree()))
                            ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
                            ax.add_feature(cfeature.BORDERS.with_scale('10m'))
                            cs=ax.plot(Coordenadas.Lon,Coordenadas.Lat, '.r',alpha=0.5)
                            plt.show()
                        itera=itera+1

                        da.close()
                        
                elif point==True:
                    [lon_min_,lat_min_]= rotated_grid_transform((lon_point,lat_point), 1, (18,-39.25))
                    if os.path.exists(path_output+v+'/'+ex+'/'+v+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+
                               direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'_'+name_point+'.csv')==True:
                        continue
                    else:
                        print('### Descargando: '+v+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'_'+name_point+'.csv')
                        da=xr.open_mfdataset(root[0]+'/'+'*.nc')
                        da=da.sel(rlat=lat_min_, rlon=lon_min_, method='nearest')
            
                        try:
                            time=da.indexes['time'].to_datetimeindex()
                        except:
                            time=da.indexes['time']
                        if v=='pr':
                            csv=pd.DataFrame((da[v].values*86400).flatten(),index=time)
                        else:
                            csv=pd.DataFrame((da[v].values-273).flatten(),index=time)
                        
                        csv.to_csv(path_output+v+'_EUR-11_'+i[0]+'_'+ex+'_'+i[2]+'_'+i[1]+'_day_'+
                                   direc[0][-20:-12]+'_'+direc[-1][-11:-3]+'_'+name_point+'.csv')

                        da.close()
