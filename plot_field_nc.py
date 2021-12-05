#!/usr/bin/env python3
###################################################################################################
## Create a plot on a map with gridded data (regularly gridded, Gaussian, fv3 tiles, vectorized 
## global, vectorized regional) stored in netcdf files, examples:
##   1. python plot_field_nc.py -i file_name.nc -v vname 
##      The code will make a geographical plot for a variable "vname". 2D variable is assumed to 
##      be stationary field data and longitude/latitude can be found in the same file. For 3D 
##      variables, the dimensions are assumed to be (lon,lat,time)
##   2. python plot_field_nc.py -i data.nc -g geo.nc -v vname -lon lon1,lon2 -lat lat1,lat2
##      The code will make a geographical plot over a specific spatial domain. The variable vname
##      is assumed to be stored in data.nc while the geographic information is stored in geo.nc.
##      If the files are located in the same directory, the directory path just needs to be
##      specified in one file.
##   3. python plot_field_nc.py -i data.nc -g geo.nc -v vname -ll lon,lat
##      The code makes a geographical plot using "lon" and "lat" as the longitude and latitude
##      variable names.
##   4. python plot_field_nc.py -i file.nc -v vname
##      The code makes a geographical plot with data stored in multiple netcdf files 
##      (file.tile*.nc). It's useful for fv3 data plotting
##   5. python plot_field_nc.py -i file -v vname -e 230,330
##      The code makes a geographical plot for the values between extremes (230-330). It's useful
##      for excluding missing values
##   6. python plot_field_nc.py -i file.nc -f same -g geo.nc -v vname
##      The code makes a geographical plot for vectorized fv3 data stored in file.nc. The information
##      cube_i/cube_j/cube_tile are stored in the same file (or specified in -f file). The global
##      lon/lat information is required and specified in geo.tile*.nc
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
import sys
import ntpath
import os
from os.path import exists

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
    if not extreme == '':
        strMinMax = extreme.split(',')
        vmin = float(strMinMax[0])
        vmax = float(strMinMax[1])
        invalid = np.logical_or(data_new == np.nan, np.logical_or(data_new > vmax, data_new < vmin))
        data_new[invalid] = np.nan
        data_new = np.ma.masked_where(np.isnan(data_new), data_new)
    vmin = np.nanmin(data_new)
    vmax = np.nanmax(data_new)
    #vmax = np.nanmean(data_new)+np.nanstd(data_new)*2
    #vmin = np.nanmean(data_new)-np.nanstd(data_new)*2
    cmap = 'viridis'
    cbarlabel = '%s' % (metadata['var'])
    if comment == '':
        plttitle = 'Variable: %s' % (metadata['var'])
    else:
        plttitle = 'Variable: %s; %s' % (metadata['var'], comment)
    for t in range(tiles):
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

def read_var(datapath, geopath, varname, tstep, llvn, vecpath):
    ntile = 6
    comment = ''
    llvns = llvn.split(',')
    if not vecpath == '':
        tiles = ntile
        nx = 96
        ny = 96
        if geopath == '':
            sys.exit("The global geographic info is required!")
        else:
            for tile in range(0, ntile):
                geofile = geopath.replace(".nc", ".tile" + str(tile+1) + ".nc")
                if not exists(geofile):
                    geofile = geopath.replace(".nc4", ".tile" + str(tile+1) + ".nc4")
                tmpgeo  = nc.Dataset(geofile, "r", format="NETCDF4")
                tmplon  = tmpgeo.variables[llvns[0]][...]
                tmplat  = tmpgeo.variables[llvns[1]][...]
                dims_geo = len(tmplon.shape)
                if dims_geo == 3:
                    lons = tmplon[0,:,:]
                    lats = tmplat[0,:,:]
                elif dims_geo == 2:
                    lons = tmplon[:,:]
                    lats = tmplat[:,:]
                elif dims_geo == 1:
                    yds_ll = len(tmplat)
                    xds_ll = len(tmplon)
                    lons = np.zeros((yds_ll, xds_ll))
                    lats = np.zeros((yds_ll, xds_ll))
                    for yid in range(yds_ll):
                        for xid in range(xds_ll):
                            lats[yid,xid] = tmplat[yid]
                            lons[yid,xid] = tmplon[xid]
                else:
                    sys.exit("The dimension "+str(dims_geo)+" is not expected!")
                if tile == 0:
                    shp_geo  = lons.shape
                    yds_geo  = shp_geo[0]
                    xds_geo  = shp_geo[1]
                    lonout   = np.zeros((yds_geo, xds_geo, tiles))
                    latout   = np.zeros((yds_geo, xds_geo, tiles))
                latout[:,:,tile]  = lats
                lonout[:,:,tile]  = lons
                tmpgeo.close()
        dataout = np.zeros((ny, nx, tiles))
        if vecpath == 'same':
            tmpvec = nc.Dataset(datapath, "r")
        else:
            tmpvec = nc.Dataset(vecpath, "r")
        cube_i    = tmpvec['cube_i'][:]
        cube_j    = tmpvec['cube_j'][:]
        cube_tile = tmpvec['cube_tile'][:]
        tmpvec.close()
#
        shp = np.shape(cube_i)
        locations = shp[0]
        dataout = np.zeros((ny, nx, tiles))
        dataout.fill(np.nan)
        ncdata = nc.Dataset(datapath, "r")
        tmpdata = ncdata[varname][...]
        dims = len(tmpdata.shape)
        if dims == 2:
            for ip in range(0, locations):
                i = cube_i[ip] - 1
                j = cube_j[ip] - 1
                tile = cube_tile[ip] - 1
                dataout[j,i,tile] = tmpdate[0,ip]
        elif dims == 1:
            for ip in range(0, locations):
                i = cube_i[ip] - 1
                j = cube_j[ip] - 1
                tile = cube_tile[ip] - 1
                dataout[j,i,tile] = tmpdata[ip]
        else:
            sys.exit("cannot deal with dimensions ("+str(dims)+") other than 1 and 2")
        ncdata.close()
        for tile in range(0, ntile):
           if tile in [2,3,4]:
               lonout[:,:,tile]  = np.rot90(lonout[:,:,tile])
               latout[:,:,tile]  = np.rot90(latout[:,:,tile])
               dataout[:,:,tile] = np.rot90(dataout[:,:,tile])
    else:
        separated_files = -1
        flag_nc4 = 0
        tiles    = 1
        print("datapath: ",datapath)
        if exists(datapath):
            separated_files = 0
        else:
            found = 1
            for tile in range(0, ntile):
                filename = datapath.replace(".nc", ".tile" + str(tile+1) + ".nc")
                if not exists(filename):
                    found = 0
            if found == 0:
                found = 1
                for tile in range(0, ntile):
                    filename = datapath.replace(".nc4", ".tile" + str(tile+1) + ".nc4")
                    if not exists(filename):
                        found = 0
                if found == 1:
                    flag_nc4 = 1
                    separated_files = 1
                    tiles = ntile
                else:
                    sys.exit("files: "+datapath+" not found")
            else:
                separated_files = 1
                tiles = ntile
        for tile in range(0, tiles):
            if tiles == 1:
                if geopath == '':
                    geofile = datapath
                else:
                    geofile = geopath
                datafile = datapath
            else:
                if geopath == '':
                    geofile = datapath.replace(".nc", ".tile" + str(tile+1) + ".nc")
                    if not exists(geofile):
                        geofile = datapath.replace(".nc4", ".tile" + str(tile+1) + ".nc4")
                else:
                    geofile = geopath.replace(".nc", ".tile" + str(tile+1) + ".nc")
                    if not exists(geofile):
                        geofile = geopath.replace(".nc4", ".tile" + str(tile+1) + ".nc4")
                datafile = datapath.replace(".nc", ".tile" + str(tile+1) + ".nc")
                if not exists(datafile):
                    datafile = datapath.replace(".nc4", ".tile" + str(tile+1) + ".nc4")
            tmpgeo  = nc.Dataset(geofile, "r", format="NETCDF4")
            tmplon  = tmpgeo.variables[llvns[0]][...]
            tmplat  = tmpgeo.variables[llvns[1]][...]
            tmpdata = nc.Dataset(datafile, "r", format="NETCDF4")
            tmpvar  = tmpdata.variables[varname][...]
            dims_data = len(tmpvar.shape)
            if dims_data == 4:
                data = tmpvar[0,0,:,:]
            elif dims_data == 3:
                data = tmpvar[0,:,:]
            elif dims_data == 2:
                data = tmpvar[:,:]
            else:
                sys.exit("The dimension "+str(dims_data)+" is not expected!")
            shp_data = data.shape
#
            dims_geo = len(tmplon.shape)
            if dims_geo == 3:
                lons = tmplon[0,:,:]
                lats = tmplat[0,:,:]
            elif dims_geo == 2:
                lons = tmplon[:,:]
                lats = tmplat[:,:]
            elif dims_geo == 1:
                yds_ll = len(tmplat)
                xds_ll = len(tmplon)
                lons = np.zeros((yds_ll, xds_ll))
                lats = np.zeros((yds_ll, xds_ll))
                for yid in range(yds_ll):
                    for xid in range(xds_ll):
                        lats[yid,xid] = tmplat[yid]
                        lons[yid,xid] = tmplon[xid]
            else:
                sys.exit("The dimension "+str(dims_geo)+" is not expected!")
            if tiles == 6 and tile in [2,3,4]:
               lats = np.rot90(lats)
               lons = np.rot90(lons)
               data = np.rot90(data)
            if tile == 0:
                shp_geo  = lons.shape
                yds_geo  = shp_geo[0]
                xds_geo  = shp_geo[1]
                shp_data = data.shape
                yds_data = shp_data[0]
                xds_data = shp_data[1]
                dataout  = np.zeros((yds_data, xds_data, tiles))
                lonout   = np.zeros((yds_geo, xds_geo, tiles))
                latout   = np.zeros((yds_geo, xds_geo, tiles))
            latout[:,:,tile]  = lats
            lonout[:,:,tile]  = lons
            dataout[:,:,tile] = data
            tmpgeo.close()
            tmpdata.close()
    return tiles, dataout, lonout, latout, comment

def gen_figure(inpath, varname, geopath, outpath, screen, tstep, llvn, lonr, latr, extreme, vector):
    # read the files to get the 2D array to plot
    tiles, data, lons, lats, comment = read_var(inpath, geopath, varname, tstep, llvn, vector)
    plotpath = outpath+'/%s.png' % (varname)
    metadata = {
                'var': varname
               }
    plot_world_map(tiles, lons, lats, data, metadata, plotpath, screen, lonr, latr, extreme, comment)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i',    '--datain',   help="path to prefix of input files (no tileN.nc)", type=str, required=True)
    ap.add_argument('-v',    '--variable', help="variable name to plot", type=str, required=True)
    ap.add_argument('-o',    '--output',   help="path to output directory", type=str, default="./")
    ap.add_argument('-g',    '--geoin',    help="path to prefix of input files for geographic info", type=str, default="")
    ap.add_argument('-s',    '--screen',   help="no if plot to file", type=str, default="yes")
    ap.add_argument('-t',    '--tstep',    help="time step for plotting", type=str, default="0")
    ap.add_argument('-lon',  '--lon',      help="longitude range for plotting", type=str, default="-180,180")
    ap.add_argument('-lat',  '--lat',      help="latitude range for plotting", type=str, default="-90,90")
    ap.add_argument('-ll',   '--llvn',     help="lonitude/latitude variable name", type=str, default="longitude,latitude")
    ap.add_argument('-e',    '--extreme',  help="minimum and maximum limits", type=str, default="")
    ap.add_argument('-f',    '--vector',   help="vectorized file or not", default="")
    MyArgs = ap.parse_args()
    print("Input data file: ",MyArgs.datain)
    if not MyArgs.vector == "":
        if MyArgs.vector == "same":
            print("Input v2f file:  ",MyArgs.datain)
        else:
            print("Input v2f file:  ",MyArgs.vector)
    if not MyArgs.geoin == "":
        print("Input geo file:  ",MyArgs.geoin)
    print("Plot variable:   ",MyArgs.variable)
    print("Lon/lat varname: ",MyArgs.llvn)
    gen_figure(MyArgs.datain, MyArgs.variable, MyArgs.geoin, MyArgs.output, MyArgs.screen, MyArgs.tstep, MyArgs.llvn, MyArgs.lon, MyArgs.lat, MyArgs.extreme, MyArgs.vector)
