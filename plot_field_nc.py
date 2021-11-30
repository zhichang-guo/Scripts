#!/usr/bin/env python3
###################################################################################################
## Create a plot on a map with gridded data stored in netcdf files, examples:
##   1. python plot_field_nc.py -i file_name.nc -v vname 
##      2D variable is assumed to be stationary field data. Assumng longitude/latitude and the
##      variable can be found in the same file. For 3D variables, the dimensions are assumed to 
##      be (lon,lat,time)
##   2. python plot_field_nc.py -i data.nc -g geo.nc -v vname -lon lon1,lon2 -lat lat1,lat2
##      The code will make a geographical plot over a specific spatial domain. The variable vname
##      is assumed to be stored in data.nc while the geographic information is stored in geo.nc.
##      If the files are located in the same directory, the directory path just needs to be
##      specified in one file.
##   3. python plot_field_nc.py -i data.nc -g geo.nc -v vname -s no -t 128
##      The code makes a geographical plot for the time step 128 and output it to a png file.
##   4. python plot_field_nc.py -i data.nc -g geo.nc -v vname -llvn lon,lat
##      The code makes a geographical plot using "lon" and "lat" as the longitude and latitude
##      variable names.
##   5. python plot_field_nc.py -i file -v vname
##      The code makes a geographical plot with data stored in multiple netcdf files (file*.nc*). 
##      It's useful for fv3 data plotting
##   6. python plot_field_nc.py -i file -v vname -e 230,330
##      The code makes a geographical plot for the values between extremes (230-330). It's useful
##      for excluding missing values
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

def plot_world_map(tiles, lons, lats, data, metadata, plotpath, screen, lonr, latr, extreme, comment):
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
    ax.set_extent([lonBeg, lonEnd, latBeg, latEnd])
    # dealing with missing values (assuming absolute values > 1.E8)
    data_new = data.copy()
    data_new[abs(data) > 1.E8] = np.nan
    vmin = np.nanmin(data_new)
    vmax = np.nanmax(data_new)
    #vmax = np.nanmean(data_new)+np.nanstd(data_new)*2
    #vmin = np.nanmean(data_new)-np.nanstd(data_new)*2
    if not extreme == '':
        strMinMax = extreme.split(',')
        vmin = float(strMinMax[0])
        vmax = float(strMinMax[1])
        invalid = np.logical_or(data > vmax, data_new < vmin)
        data_new[invalid] = np.nan
        data_new = np.ma.masked_where(np.isnan(data_new), data_new)
    cmap = 'viridis'
    cbarlabel = '%s' % (metadata['var'])
    if comment == '':
        plttitle = 'Variable: %s' % (metadata['var'])
    else:
        plttitle = 'Variable: %s; %s' % (metadata['var'], comment)
    for t in range(0,tiles):
        cs = ax.pcolormesh(lons[...,t], lats[...,t], data_new[...,t],vmin=vmin,vmax=vmax,cmap=cmap)
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

def read_var(datapath, geopath, varname, tstep, llvn):
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
    tile = 0
    llvns = llvn.split(',')
    for f in geofiles:
        tile += 1
        if tile == 1:
            tmpdata = nc.Dataset(f,'r')
            tmplon = tmpdata.variables[llvns[0]][:]
            tmplat = tmpdata.variables[llvns[1]][:]
            tmpdata.close()
    tiles = tile
    dims_lon = len(tmplon.shape)
    dims_lat = len(tmplat.shape)
    if dims_lon == 1 and dims_lat == 1:
        yds = len(tmplat)
        xds = len(tmplon)
        arrayshape = (yds,) + (xds,) + (tiles,)
    elif dims_lon == 2 and dims_lat == 2:
        arrayshape = tmplat.shape + (tiles,)
    elif dims_lon == 3 and dims_lat == 3:
        arrayshape = tmplat[0,...].shape + (tiles,)
    else:
        sys.exit("Error: cannot handle variables with dimensions more than 3 or less than 1 for lon/lat")

    dataout = np.empty(arrayshape)
    lonout = np.empty(arrayshape)
    latout = np.empty(arrayshape)
    for t in range(0,tiles):
        datafile = obsfiles[t]
        geofile = geofiles[t]
        geonc = nc.Dataset(geofile)
        datanc = nc.Dataset(datafile)
        lats = geonc.variables[llvns[1]][:]
        lons = geonc.variables[llvns[0]][:]
        data = datanc.variables[varname][:]
        dims = len(data.shape)
        tds = len(data)
        if tiles == 6 and t in [2,3,4]:
            lats = np.rot90(lats)
            lons = np.rot90(lons)
            if dims == 2:
                data = np.rot90(data)
            elif dims == 3:
                for i in range(tds):
                    data[i,...] = np.rot90(data[i,...])
        if dims_lon == 1 and dims_lat == 1:
            yds = len(lats)
            xds = len(lons)
            for yid in range(yds):
                for xid in range(xds):
                    latout[yid,xid,t-1] = lats[yid]
                    lonout[yid,xid,t-1] = lons[xid]
        elif dims_lon == 2 and dims_lat == 2:
            latout[:,:,t-1] = lats
            lonout[:,:,t-1] = lons
        elif dims_lon == 3 and dims_lat == 3:
            latout[:,:,t-1] = lats[0,...]
            lonout[:,:,t-1] = lons[0,...]
        if dims == 2:
            dataout[:,:,t-1] = data
        elif dims == 3:
            yds = len(data[0])
            xds = len(data[0][0])
            data_new = np.zeros((yds,xds))
            if ',' in tstep:
                tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                tsteps = tstepEnd - tstepBeg + 1
                for yid in range(yds):
                    for xid in range(xds):
                        cnt = 0.0
                        for i in range(tsteps):
                            tid = i + tstepBeg
                            data_new[yid,xid] += data[tid][yid][xid]
                            cnt += 1.0
                        if cnt > 0.5:
                            data_new[yid,xid] /= cnt
                        else:
                            data_new[yid,xid] = np.nan
                comment = 'Time average: %s - %s' % (str(tstepBeg), str(tstepEnd))
            else:
                timestep = int(tstep)
                if timestep < 0:
                    timestep += len(datatmp)
                timestep = max(timestep,0)
                timestep = min(timestep,tds-1)
                data_new = data[timestep,:,:]
                comment = 'Time step: %s of 0 - %s' % (str(timestep), str(tds-1))
            dataout[:,:,t-1] = data_new
        else:
            sys.exit("Error: cannot handle variables with dimensions more than 3 or less than 2")
        geonc.close()
        datanc.close()
    return tiles, dataout, lonout, latout, comment

def gen_figure(inpath, geopath, outpath, varname, screen, tstep, llvn, lonr, latr, extreme):
    # read the files to get the 2D array to plot
    tiles, data, lons, lats, comment = read_var(inpath, geopath, varname, tstep, llvn)
    plotpath = outpath+'/%s.png' % (varname)
    metadata = {
                'var': varname
               }
    plot_world_map(tiles, lons, lats, data, metadata, plotpath, screen, lonr, latr, extreme, comment)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-i', '--input', help="path to prefix of input files (no tileN.nc)", required=True)
    ap.add_argument('-g', '--geo', help="path to prefix of input files for geographic info", default="")
    ap.add_argument('-v', '--variable', help="variable name to plot", required=True)
    ap.add_argument('-s', '--screen', help="no if plot to file", default="yes")
    ap.add_argument('-t', '--tstep', help="time step for plotting", default="0")
    ap.add_argument('-lon', '--lon', help="longitude range for plotting", default="-180,180")
    ap.add_argument('-lat', '--lat', help="latitude range for plotting", default="-90,90")
    ap.add_argument('-llvn', '--llvname', help="lonitude/latitude variable name", default="longitude,latitude")
    ap.add_argument('-e', '--extreme', help="minimum and maximum limits", default="")
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.variable, MyArgs.screen, MyArgs.tstep, MyArgs.llvname, MyArgs.lon, MyArgs.lat, MyArgs.extreme)
