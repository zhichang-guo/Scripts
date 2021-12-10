#!/usr/bin/env python3
###################################################################################################
## Create a plot on a map with irregularly spaced data stored in netcdf files, examples:
##   1. python plot_snapshot_nc.py -i file_name.nc -v vname 
##      1D variable is assumed to be stationary vector data. Also assume lon/lat and the variable
##      can be found in the same file. For 2D(vector,time), 2D(lon,lat), 3D(lon,lat,time), 
##      3D(vector,level,time), 4D(lon,lat,lev,time) variable, the code creates a geographical plot
##      with data averaged over the time domain.
##   2. python plot_snapshot_nc.py -i data_name.nc -v vname1-vname2 -t 2,-3 -factor 86400
##      Make a graphical plot of difference, multiplied by 86400, between vname1 and vname2 and
##      averaged over the time segment from the second record to the third last record.
##   3. python plot_snapshot_nc.py -i data.nc -g geo.nc -v vname -lon lon1,lon2 -lat lat1,lat2
##      The code will make a geographical plot over a specific spatial domain. The variable vname
##      is assumed to be stored in data.nc while the geographic information is stored in geo.nc.
##      If the files are located in the same directory, the directory path just needs to be
##      specified in one file.
##   4. python plot_snapshot_nc.py -i data.nc -g geo.nc -v vname -s no -t 128
##      The code makes a geographical plot for the time step 128 and output it to a png file.
##   5. python plot_snapshot_nc.py -i data.nc -g geo.nc -v vname -llvn lon,lat
##      The code makes a geographical plot using "lon" and "lat" as the longitude and latitude
##      variable names.
##   6. python plot_snapshot_nc.py -i file -v vname
##      The code makes a geographical plot with data stored in multiple netcdf files (file*.nc*)
##   7. python plot_snapshot_nc.py -i file -v vname -f samefield -r radian
##      The code makes a geographical plot assuming lon and lat having the same format with data 
##      stored in netcdf files (x-axis,y-axis,tiles) and lon/lat is in radian
##   8. python plot_snapshot_nc.py -i file -v vname -f samefield -r radian -z 2
##      The code makes a geographical plot at the vertical level 2 (starting from 0)
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
###################################################################################################
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
import argparse
import glob
import os
import sys
import ntpath

def lonF2S(flon):
    flon = int(flon*100)/100.
    if flon > 180:
        flon = 360. - flon
        flon = int(flon*100)/100.
        strLon = str(flon) + 'W'
    elif flon < 0.0:
        strLon = str(-flon) + 'W'
    else:
        strLon = str(flon) + 'E'
    return strLon

def lonS2F(strLon):
    strLon = strLon.upper()
    strLon.replace('E','')
    if 'W' in strLon:
        flon = -1.0*float(strLon.replace('W',''))
    else:
        flon = float(strLon)
    return flon

def latF2S(flat):
    flat = int(flat*100)/100.
    if flat < 0.0:
        strLat = str(-flat) + 'S'
    else:
        strLat = str(flat) + 'N'
    return strLat

def latS2F(strLat):
    strLat = strLat.upper()
    strLat.replace('N','')
    if 'S' in strLat:
        flat = -1.0*float(strLat.replace('S',''))
    else:
        flat = float(strLat)
    return flat

def tstepS2I(strTStep, tds):
    if ',' in strTStep:
        tsteps = strTStep.split(',')
        tstepBeg = int(tsteps[0])
        tstepEnd = int(tsteps[1])
        if tstepBeg < 0:
            tstepBeg += tds
        if tstepEnd < 0:
            tstepEnd += tds
        tstepBeg = max(0,min(tstepBeg,tds-1))
        tstepEnd = max(0,min(tstepEnd,tds-1))
    else:
        tstepBeg = 0
        tstepEnd = tds-1
    if tstepBeg > tstepEnd:
        return tstepEnd, tstepBeg
    else:
        return tstepBeg, tstepEnd

def last_6chars(x):
    return(x[-6:])

def plot_world_map(lons, lats, data, metadata, plotpath, screen, lonr, latr, extreme, comment):
    # plot generic world map
    if screen.upper() == "NO":
        matplotlib.use('agg')
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    listWE = lonr.split(',')
    listSN = latr.split(',')
    lonBeg = lonS2F(listWE[0])
    lonEnd = lonS2F(listWE[1])
    latBeg = latS2F(listSN[0])
    latEnd = latS2F(listSN[1])
    scaleWE = 360/(lonEnd-lonBeg)*10
    scaleSN = 180/(latEnd-latBeg)*10
    scaleMax = max(scaleWE,scaleSN)
    ax.set_extent([lonBeg, lonEnd, latBeg, latEnd])
    if not extreme == '':
        strMinMax = extreme.split(',')
        vmin_c = float(strMinMax[0])
        vmax_c = float(strMinMax[1])
        invalid = np.logical_or(data > vmax_c, data < vmin_c)
        data[invalid] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)
    vmax = np.nanmean(data)+np.nanstd(data)*2
    vmin = np.nanmean(data)-np.nanstd(data)*2
    print("vmax, vmin: ",vmax, ' ', vmin)
    cmap = 'viridis'
    cbarlabel = '%s' % (metadata['var'])
    if comment == '':
        plttitle = 'Variable: %s' % (metadata['var'])
    else:
        plttitle = 'Variable: %s; %s' % (metadata['var'], comment)
    cs = plt.scatter(lons, lats, c=data, s=scaleMax,
                     cmap=cmap, transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax)
    cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
    cb.set_label(cbarlabel, fontsize=12)
    plt.title(plttitle)
    if screen.upper() == "NO" or screen.upper() == "ALL":
        plt.savefig(plotpath)
        if screen.upper() == "ALL":
            plt.show()
        plt.close('all')
    else:
        plt.show()

def read_var(datapath, geopath, varname, tstep, fov, llvn, fact, radian, level):
    data_name, data_extension = os.path.splitext(datapath)
    if geopath == "":
        geo_name, geo_extension = os.path.splitext(datapath)
    else:
        geo_name, geo_extension = os.path.splitext(geopath)
    if not 'nc' in data_extension:
        datapath += '*'
    if not 'nc' in geo_extension:
        geopath += '*'
    opath, obsfname = ntpath.split(datapath)
    if geopath == "":
        gpath = ''
        geofname = ''
    else:
        gpath, geofname = ntpath.split(geopath)
    if opath == '' and gpath == '':
        cpath = os.getcwd()
        obsfiles = glob.glob(os.path.join(cpath,obsfname))
        if geofname == '':
            geofiles = glob.glob(os.path.join(cpath,obsfname))
        else:
            geofiles = glob.glob(os.path.join(cpath,geofname))
    elif opath != '' and gpath == '':
        obsfiles = glob.glob(datapath)
        if geofname == '':
            geofiles = obsfiles
        else:
            geofiles = glob.glob(os.path.join(opath,geofname))
    elif opath == '' and gpath != '':
        obsfiles = glob.glob(os.path.join(gpath,datapath))
        geofiles = glob.glob(geopath)
    else:
        obsfiles = glob.glob(datapath)
        if geofname == '':
            geofiles = obsfiles
        else:
            geofiles = glob.glob(geopath)
    geofiles = sorted(geofiles, key=last_6chars)
    obsfiles = sorted(obsfiles, key=last_6chars)
    comment = ''
    llvns = llvn.split(',')
    lats = np.array([])
    lons = np.array([])
    data = np.array([])
    zid = 0
    if not level == '':
        zid = int(level)
    if not geopath == "":
        if 'SAME' in fov.upper():
            for g in geofiles:
                geonc = nc.Dataset(g)
                lontmp = geonc.variables[llvns[0]][:]
                lattmp = geonc.variables[llvns[1]][:]
                dim_ll = len(lontmp.shape)
                if dim_ll == 3:
                    for tid in range(len(lontmp)):
                        lats = np.concatenate((lats,lattmp[tid].ravel()))
                        lons = np.concatenate((lons,lontmp[tid].ravel()))
                elif dim_ll == 2:
                    lats = np.concatenate((lats,lattmp.ravel()))
                    lons = np.concatenate((lons,lontmp.ravel()))
                else:
                    lats = np.concatenate((lats,lattmp))
                    lons = np.concatenate((lons,lontmp))
                geonc.close()
        else:
            for g in geofiles:
                geonc = nc.Dataset(g)
                lontmp = geonc.variables[llvns[0]][:]
                lattmp = geonc.variables[llvns[1]][:]
                lats = np.concatenate((lats,lattmp))
                lons = np.concatenate((lons,lontmp))
                geonc.close()
    for f in obsfiles:
        datanc = nc.Dataset(f)
        if geopath == "":
            lontmp = datanc.variables[llvns[0]][:]
            lattmp = datanc.variables[llvns[1]][:]
            if len(lontmp.shape) == 2:
                lons = np.concatenate((lons,lontmp.ravel()))
                lats = np.concatenate((lats,lattmp.ravel()))
            else:
                if 'FIELD' in fov.upper():
                    xds = len(lontmp)
                    yds = len(lattmp)
                    lds = xds*yds
                    lons = np.zeros(lds)
                    lats = np.zeros(lds)
                    for j in range(yds):
                        for i in range(xds):
                            k = j*xds + i
                            lons[k] = lontmp[i]
                            lats[k] = lattmp[j]
                else:
                    lats = np.concatenate((lats,lattmp))
                    lons = np.concatenate((lons,lontmp))
        if '-' in varname:
            varnames = varname.split('-')
            datatmp = datanc.variables[varnames[0]][:]
            datatmp2 = datanc.variables[varnames[1]][:]
            datatmp -= datatmp2
        else:
            datatmp = datanc.variables[varname][:]
        datanc.close()
        dims = len(datatmp.shape)
        if dims == 1:
            data = np.concatenate((data,datatmp))
            comment = ''
        elif dims == 2:
            if 'VECTOR' in fov.upper():
                tds = len(datatmp)
                lds = len(datatmp[0])
                if ',' in tstep:
                    tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                    tsteps = tstepEnd - tstepBeg + 1
                    data_new = np.zeros(lds)
                    for lid in range(lds):
                        cnt = 0.0
                        for i in range(tsteps):
                            tid = i + tstepBeg
                            data_new[lid] += datatmp[tid][lid]
                            cnt += 1.0
                        if cnt > 0.5:
                            data_new[lid] /= cnt
                        else:
                            data_new[lid] = np.nan
                    data = np.concatenate((data,data_new))
                    comment = 'Time average: %s - %s' % (str(tstepBeg), str(tstepEnd))
                else:
                    timestep = int(tstep)
                    if timestep < 0:
                        timestep += len(datatmp)
                    timestep = max(timestep,0)
                    timestep = min(timestep,tds-1)
                    data = np.concatenate((data,datatmp[timestep]))
                    comment = 'Time step: %s of 0 - %s' % (str(timestep), str(tds-1))
            elif 'FIELD' in fov.upper():
                data = np.concatenate((data,datatmp.ravel()))
                comment = 'Stationary Field'
            else:
                sys.exit("Error: invalid fov option")
        elif dims == 3:
            if 'VECTOR' in fov.upper():
                tds = len(datatmp)
                zds = len(datatmp[0])
                lds = len(datatmp[0][0])
                zid = max(0,min(zid,zds-1))
                if ',' in tstep:
                    tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                    tsteps = tstepEnd - tstepBeg + 1
                    data_new = np.zeros(lds)
                    for lid in range(lds):
                        cnt = 0.0
                        for i in range(tsteps):
                            tid = i + tstepBeg
                            data_new[lid] += datatmp[tid][zid][lid]
                            cnt += 1.0
                        if cnt > 0.5:
                            data_new[lid] /= cnt
                        else:
                            data_new[lid] = np.nan
                    data = np.concatenate((data,data_new))
                    comment = 'Time average: %s - %s; Level: %s' % (str(tstepBeg), str(tstepEnd), str(zid))
                else:
                    timestep = int(tstep)
                    if timestep < 0:
                        timestep += len(datatmp)
                    timestep = max(timestep,0)
                    timestep = min(timestep,tds-1)
                    data = np.concatenate((data,datatmp[timestep][zid]))
                    comment = 'Time step: %s of 0 - %s; Level: %s' % (str(timestep), str(tds-1), str(zid))
            elif 'FIELD' in fov.upper():
                tds = len(datatmp)
                yds = len(datatmp[0])
                xds = len(datatmp[0][0])
                data_new = np.zeros((yds*xds))
                if ',' in tstep:
                    tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                    tsteps = tstepEnd - tstepBeg + 1
                    for yid in range(yds):
                        for xid in range(xds):
                            lid = yid*xds + xid
                            cnt = 0.0
                            for i in range(tsteps):
                                tid = i + tstepBeg
                                data_new[lid] += datatmp[tid][yid][xid]
                                cnt += 1.0
                            if cnt > 0.5:
                                data_new[lid] /= cnt
                            else:
                                data_new[lid] = np.nan
                    comment = 'Time average: %s - %s' % (str(tstepBeg), str(tstepEnd))
                else:
                    timestep = int(tstep)
                    if timestep < 0:
                        timestep += len(datatmp)
                    timestep = max(timestep,0)
                    timestep = min(timestep,tds-1)
                    data_new = datatmp[timestep].ravel()
                    comment = 'Time step: %s of 0 - %s' % (str(timestep), str(tds-1))
                data = np.concatenate((data,data_new))
            else:
                sys.exit("Error: invalid fov option")
        elif dims == 4:
            tds = len(datatmp)
            zds = len(datatmp[0])
            yds = len(datatmp[0][0])
            xds = len(datatmp[0][0][0])
            zid = max(0,min(zid,zds-1))
            data_new = np.zeros((yds*xds))
            if ',' in tstep:
                tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                tsteps = tstepEnd - tstepBeg + 1
                for yid in range(yds):
                    for xid in range(xds):
                        lid = yid*xds + xid
                        cnt = 0.0
                        for i in range(tsteps):
                            tid = i + tstepBeg
                            data_new[lid] += datatmp[tid][zid][yid][xid]
                            cnt += 1.0
                        if cnt > 0.5:
                            data_new[lid] /= cnt
                        else:
                            data_new[lid] = np.nan
                comment = 'Time average: %s - %s; Level: ' % (str(tstepBeg), str(tstepEnd), str(zid))
            else:
                timestep = int(tstep)
                if timestep < 0:
                    timestep += len(datatmp)
                timestep = max(timestep,0)
                timestep = min(timestep,tds-1)
                data_new = datatmp[timestep][zid].ravel()
                comment = 'Time step: %s of 0 - %s; Level: %s' % (str(timestep), str(tds-1), str(zid))
            data = np.concatenate((data,data_new))
        else:
            sys.exit("cannot handle variables with dimensions more than 4")
    if not fact == '1':
        factor = float(fact)
        for lid in range(len(data)):
            data[lid] *= factor
        comment += '; Factor: %s'%(fact)
    if 'RADIAN' in radian.upper():
        rad2deg = 180.0/3.1415926535
        for lid in range(len(lons)):
            lons[lid] *= rad2deg
        for lid in range(len(lats)):
            lats[lid] *= rad2deg
    return data, lons, lats, comment

def gen_figure(inpath, geopath, outpath, varname, screen, tstep, fov, llvn, lonr, latr, fact, radian, level, extreme):
    # read the files to get the 2D array to plot
    data, lons, lats, comment = read_var(inpath, geopath, varname, tstep, fov, llvn, fact, radian, level)
    plotpath = outpath+'/%s.png' % (varname)
    metadata = {
                'var': varname
                }
    plot_world_map(lons, lats, data, metadata, plotpath, screen, lonr, latr, extreme, comment)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-g', '--geo', help="path to the geographic info file", default="")
    ap.add_argument('-v', '--variable', help="variable name to plot", required=True)
    ap.add_argument('-s', '--screen', help="no if plot to file", default="yes")
    ap.add_argument('-t', '--tstep', help="time step for plotting", default="0,-1")
    ap.add_argument('-f', '--fov', help="field or vector data", default="vector")
    ap.add_argument('-lon', '--longitude', help="longitude range for plotting", default="-180,180")
    ap.add_argument('-lat', '--latitude', help="latitude range for plotting", default="-90,90")
    ap.add_argument('-llvn', '--llvname', help="lonitude/latitude variable name", default="longitude,latitude")
    ap.add_argument('-factor', '--fact', help="factor for units conversion", default="1")
    ap.add_argument('-e', '--extreme', help="minimum and maximum limits", default="")
    ap.add_argument('-r', '--radian', help="radian or degree for lon/lat", default="degree")
    ap.add_argument('-z', '--level', help="vertial level index", default="0")
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.variable, MyArgs.screen, MyArgs.tstep, MyArgs.fov, MyArgs.llvname, MyArgs.longitude, MyArgs.latitude, MyArgs.fact, MyArgs.radian, MyArgs.level, MyArgs.extreme)
