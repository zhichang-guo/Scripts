#!/usr/bin/env python3
#################################################################################################
## Create a plot on a map with irregularly spaced data stored in netcdf files, examples:
##   1. python plot_snapshot_nc.py -i file_name.nc -v vname 
##      1D variable is assumed to be stationary vector data. Also assume lon/lat and the variable
##      can be found in the same file. For 2D(vector,time) and 3D(lon,lat,time) variable, the 
##      code creates a geographical plot with data averaged over the time domain.
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
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
#################################################################################################
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
            tstepBeg += tds + 1
        if tstepEnd < 0:
            tstepEnd += tds + 1
        tstepBeg = min(tstepBeg,tds)
        tstepEnd = min(tstepEnd,tds)
    else:
        tstepBeg = 0
        tstepEnd = tds
    if tstepBeg > tstepEnd:
        return tstepEnd, tstepBeg
    else:
        return tstepBeg, tstepEnd

def plot_world_map(lons, lats, data, metadata, plotpath, screen, lonr, latr, comment):
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
    vmax = np.nanmean(data)+np.nanstd(data)*2
    vmin = np.nanmean(data)-np.nanstd(data)*2
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
    if screen.upper() == "NO":
        plt.savefig(plotpath)
        plt.close('all')
    else:
        plt.show()

def read_var(datapath, geopath, varname, tstep, fov, llvn, fact):
    file_name, file_extension = os.path.splitext(datapath)
    obsfiles = glob.glob(datapath)
    if not 'nc' in file_extension:
        obsfiles = glob.glob(datapath+'*')
    geofile  = glob.glob(geopath)
    opath, obsfname = ntpath.split(datapath)
    gpath, geofname = ntpath.split(geopath)
    geofile_new  = geofile
    obsfiles_new = obsfiles
    comment = ''
    if opath == '' and gpath == '':
        cpath = os.getcwd()
        obsfiles_new = glob.glob(os.path.join(cpath,datapath))
        geofile_new = glob.glob(os.path.join(cpath,geofname))
    elif opath != '' and gpath == '':
        geofile_new = np.append(geofile_new, [os.path.join(opath,geofname)])
    elif opath == '' and gpath != '':
        obsfiles_new = glob.glob(os.path.join(gpath,datapath))
    llvns = llvn.split(',')
    lats = np.array([])
    lons = np.array([])
    data = np.array([])
    if not geopath == "":
        for g in geofile_new:
            geonc = nc.Dataset(g)
            lontmp = geonc.variables[llvns[0]][:]
            lattmp = geonc.variables[llvns[1]][:]
            lats = np.concatenate((lats,lattmp))
            lons = np.concatenate((lons,lontmp))
            geonc.close()
    for f in obsfiles_new:
        datanc = nc.Dataset(f)
        if geopath == "":
            lontmp = datanc.variables[llvns[0]][:]
            lattmp = datanc.variables[llvns[1]][:]
            lats = np.concatenate((lats,lattmp))
            lons = np.concatenate((lons,lontmp))
        if '-' in varname:
            varnames = varname.split('-')
            datatmp = datanc.variables[varnames[0]][:]
            datatmp2 = datanc.variables[varnames[1]][:]
            if len(datatmp.shape) == 1:
                for lid in range(len(datatmp)):
                    datatmp[lid] -= datatmp2[lid]
            elif len(datatmp.shape) == 2:
                for tid in range(len(datatmp)):
                    for lid in range(len(datatmp[0])):
                        datatmp[tid][lid] -= datatmp2[tid][lid]
            elif len(datatmp.shape) == 3:
                for tid in range(len(datatmp)):
                    for yid in range(len(datatmp[0])):
                        for xid in range(len(datatmp[0][0])):
                            datatmp[tid][yid][xid] -= datatmp2[tid][yid][xid]
            else:
                sys.exit("cannot handle variables with dimensions more than 3 (lon,lat,time or vector,time or time)")
        else:
            datatmp = datanc.variables[varname][:]
        datanc.close()
        dims = len(datatmp.shape)
        if dims == 1:
            data = np.concatenate((data,datatmp))
            comment = ''
        elif dims == 2:
            if 'VECTOR' in fov.upper():
                tds = len(datatmp) - 1
                lds = len(datatmp[0])
                if ',' in tstep:
                    tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                    tsteps = tstepEnd - tstepBeg + 1
                    data = np.zeros(lds)
                    for lid in range(lds):
                        cnt = 0.0
                        for i in range(tsteps):
                            tid = i + tstepBeg
                            data[lid] += datatmp[tid][lid]
                            cnt += 1.0
                        if cnt > 0.5:
                            data[lid] /= cnt
                        else:
                            data[lid] = -999999.99
                    comment = 'Time average: %s - %s' % (str(tstepBeg), str(tstepEnd))
                else:
                    timestep = int(tstep)
                    if timestep < 0:
                        timestep += len(datatmp)
                    timestep = max(timestep,0)
                    timestep = min(timestep,tds)
                    data = np.concatenate((data,datatmp[timestep]))
                    comment = 'Time step: %s of 0 - %s' % (str(timestep), tds)
            elif 'FIELD' in fov.upper():
                yds = len(datatmp)
                xds = len(datatmp[0])
                lons_new = np.zeros((yds*xds))
                lats_new = np.zeros((yds*xds))
                data_new = np.zeros((yds*xds))
                for yid in range(yds):
                    for xid in range(xds):
                        lid = yid*xds + xid
                        lons_new[lid] = lons[xid]
                        lats_new[lid] = lats[yid]
                        data_new[lid] = datatmp[yid][xid]
                lons = lons_new
                lats = lats_new
                data = data_new
                comment = 'Stationary Field'
            else:
                sys.exit("Error: invalid fov option")
        elif dims == 3:
            tds = len(datatmp) - 1
            yds = len(datatmp[0])
            xds = len(datatmp[0][0])
            if ',' in tstep:
                tstepBeg, tstepEnd = tstepS2I(tstep, tds)
                tsteps = tstepEnd - tstepBeg + 1
                data = np.zeros((yds,xds))
                for yid in range(yds):
                    for xid in range(xds):
                        cnt = 0.0
                        for i in range(tsteps):
                            tid = i + tstepBeg
                            data[yid][xid] += datatmp[tid][yid][xid]
                            cnt += 1.0
                        if cnt > 0.5:
                            data[yid][xid] /= cnt
                        else:
                            data[yid][xid] = -999999.99
                comment = 'Time average: %s - %s' % (str(tstepBeg), str(tstepEnd))
            else:
                timestep = int(tstep)
                if timestep < 0:
                    timestep += len(datatmp)
                timestep = max(timestep,0)
                timestep = min(timestep,tds)
                data = np.concatenate((data,datatmp[timestep]))
                comment = 'Time step: %s of 0 - %s' % (str(timestep), tds)
            lons_new = np.zeros((yds*xds))
            lats_new = np.zeros((yds*xds))
            data_new = np.zeros((yds*xds))
            for yid in range(yds):
                for xid in range(xds):
                    lid = yid*xds + xid
                    lons_new[lid] = lons[xid]
                    lats_new[lid] = lats[yid]
                    data_new[lid] = data[yid][xid]
            lons = lons_new
            lats = lats_new
            data = data_new
        else:
            sys.exit("cannot handle variables with dimensions more than 2 (vector, time)")
        if not fact == '1':
            factor = float(fact)
            for lid in range(len(data)):
                data[lid] *= factor
            comment += '; Factor: %s'%(fact)
    return data, lons, lats, comment

def gen_figure(inpath, geopath, outpath, varname, screen, tstep, fov, llvn, lonr, latr, fact):
    # read the files to get the 2D array to plot
    data, lons, lats, comment = read_var(inpath, geopath, varname, tstep, fov, llvn, fact)
    plotpath = outpath+'/%s.png' % (varname)
    metadata = {
                'var': varname
                }
    plot_world_map(lons, lats, data, metadata, plotpath, screen, lonr, latr, comment)

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
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.variable, MyArgs.screen, MyArgs.tstep, MyArgs.fov, MyArgs.llvname, MyArgs.longitude, MyArgs.latitude, MyArgs.fact)
