#!/usr/bin/env python3
###############################################################################################
## plot time series of energy balance terms stored in netcdf files, examples:
##   Note: the code is not generic!
##   1. python plot_fluxnet_energy_nc.py -i file_name.nc
##      1D variable is assumed to be time-dependent. The code draws time series of all terms.
##      2D(vector,time) variable, it draws time series of all terms averaged over all points.
##      3D(lon,lat,time) variable. The code draws domain average time series of all terms.
##   2. python plot_energy_nc.py -i data_name.nc -t 2,-2
##      Draw time series of energy balance terms for the time segment from the second record
##      to the second last record.
##   3. python plot_energy_nc.py -i data_name.nc -g geo_name -lon lon -lat lat
##      The code will find the grid point which is closest to the given location(lon/lat) and 
##      draw the time series. The variables are assumed to be stored in data_name.nc while
##      the geographic information is stored in geo_name. If the files are located in the same
##      directory, the directory path just needs to be specified in one file name.
##   4. python plot_energy_nc.py -i data_name.nc -lon lon1,lon2 -lat lat1,lat2
##      The code will average the variables over the domain (lon1,lon2),(lat1,lat2) and draw
##      the time series
##   5. python plot_energy_nc.py -i file.nc -v Rnet,hflx,evap,gflux,snohf -l dotted,,,,,dashed
##      The code will draw the time series of each budget term, with Rnet line dotted and
##      residual line dashed. Also can use "-v swnet,lwdown,lwup,hflx,evap,gflux,snohf"
##      default setting: "-v Rnet,-hflx,-evap,gflux,snohf"
##   6. python plot_energy_nc.py -i data_name.nc -p point_index -s no
##      The code draws the time series and output it to a png file.
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
###############################################################################################
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
from datetime import timedelta
from datetime import datetime

def plot_time_series(dataX, timeX, dataY, metadata, plotpath, screen, varnames, comment, color, style):
#   defaultColors = np.array(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    defaultColors = np.array(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f','#bcbd22', '#17becf'])
    if screen.upper() == "NO":
        matplotlib.use('agg')
    fig = plt.figure(figsize=(12,8))
    if comment == '':
        plttitle = 'Variable: %s' % (metadata['var'])
    else:
        plttitle = comment
    colors = color.split(',')
    styles = style.split(',')
    r = np.zeros(len(dataY[0]))
    for i in range(len(dataY)):
        factor = 1.
        if varnames[i] in ['lwup', 'hflx', 'evap']:
            factor = -1.
        for lpt in range(len(dataY[0])):
             r[lpt] += factor * dataY[i][lpt]               
#   varnames.append('residual')
#   dataX = np.insert(dataX,len(dataX),dataX[0],0)
#   timeX = np.insert(timeX,len(timeX),timeX[0],0)
#   dataY = np.insert(dataY,len(dataY),r,0)
    for i in range(len(dataY)):
        time = timeX[i]
        x = dataX[i]
        y = dataY[i]
        c = defaultColors[min(i,len(defaultColors)-1)]
        s = 'solid'
        if i == len(dataY) - 1:
            s = 'dotted'
        if len(colors) > i and colors[i] != '':
            c = colors[i]
        if len(styles) > i and styles[i] != '':
            s = styles[i]
        plt.plot(time, y, color=c, linestyle=s, label=varnames[i])
    plt.legend()
    plt.title(plttitle)
    if screen.upper() == "NO":
        plt.savefig(plotpath)
        plt.close('all')
    else:
        plt.show()

def find_nearest_point(strLon, strLat, lons, lats):
    xlon = float(strLon)
    xlat = float(strLat)
    if xlon < 0.0:
        xlon += 360.
    xlona = xlon - 360.
    xlonb = xlon + 360.
    distanceMin = 9999.9
    index = -1
    found = -1
    dist  = -9.9
    for lid in range(len(lons)):
        if lons[lid] < 0.0:
            lons[lid] += 360.
        distLat = (xlat-lats[lid])*(xlat-lats[lid])
        dist  = (xlon-lons[lid])*(xlon-lons[lid]) + distLat
        distA = (xlona-lons[lid])*(xlona-lons[lid]) + distLat
        distB = (xlonb-lons[lid])*(xlonb-lons[lid]) + distLat
        if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
            distanceMin = min(dist,min(distA,distB))
            index = lid
            found = 1
    if found > 0:
        return index
    else:
        print(strLon+" "+strLat+" "+str(dist)+" "+str(len(lons)))
        sys.exit("the nearest point is not found")

def find_nearest_point2D(strLon, strLat, lons, lats):
    xlon = float(strLon)
    xlat = float(strLat)
    if xlon < 0.0:
        xlon += 360.
    xlona = xlon - 360.
    xlonb = xlon + 360.
    distanceMin = 9999.9
    index = -1
    found = -1
    for yid in range(len(lats)):
        distLat = (xlat-lats[yid])*(xlat-lats[yid])
        for xid in range(len(lons)):
            if lons[xid] < 0.0:
                lons[xid] += 360.
            dist  = (xlon-lons[xid])*(xlon-lons[xid]) + distLat
            distA = (xlona-lons[xid])*(xlona-lons[xid]) + distLat
            distB = (xlonb-lons[xid])*(xlonb-lons[xid]) + distLat
            if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
                distanceMin = min(dist,min(distA,distB))
                i_index = xid
                j_index = yid
                found = 1
    if found > 0:
        return i_index, j_index
    else:
        sys.exit("the nearest point is not found")

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
        if tstepEnd < tstepBeg:
            tstepTmp = tstepBeg
            tstepBeg = tstepEnd
            tstepEnd = tstepTmp
    else:
        tstepBeg = 0
        tstepEnd = tds
    return tstepBeg, tstepEnd

def read_single_var(datanc, varname, tstep, lonr, latr, lpt, lons, lats):
    sigma = 5.67E-8
    dataX = np.array([])
    timeX = np.array([])
    dataY = np.array([])
    varname = varname.replace(' ','')
    choice = 0
    if varname.lower() == 'rnet':
        if choice == 0:
            datatmp = datanc.variables['NETRAD'][:]
        else:
            datatmp = datanc.variables['net_radiation'][:]
    elif varname.lower() == 'le':
        if choice == 0:
            datatmp = datanc.variables['LE_F_MDS'][:]
        else:
            datatmp = datanc.variables['latent_heat'][:]
    elif varname.lower() == 'hf':
        if choice == 0:
            datatmp = datanc.variables['H_F_MDS'][:]
        else:
            datatmp = datanc.variables['sensible_heat'][:]
    elif varname.lower() == 'gf':
        if choice == 0:
            datatmp = datanc.variables['G_F_MDS'][:]
        else:
            datatmp = datanc.variables['ground_heat'][:]
    elif varname.lower() == 'swnet':
        if choice == 0:
            datatmp = datanc.variables['SW_IN_F_MDS'][:]
            tmp = datanc.variables['SW_OUT'][:]
        else:
            datatmp = datanc.variables['sw_down'][:]
            tmp = datanc.variables['sw_up'][:]
        if len(tmp.shape) == 1:
            for lid in range(len(tmp)):
                datatmp[lid] -= tmp[lid]
        elif len(tmp.shape) == 2:
            for tid in range(len(tmp)):
                for lid in range(len(tmp[0])):
                    datatmp[tid][lid] -= tmp[tid,lid]
        elif len(tmp.shape) == 3:
            for tid in range(len(tmp)):
                for yid in range(len(tmp[0])):
                    for xid in range(len(tmp[0][0])):
                        datatmp[tid][yid][xid] -= tmp[tid,yid,xid]
        else:
            sys.exit("cannot handle variables with dimensions more than 3 (lon, lat, time or vector, time)")
    elif varname.lower() == 'lwnet':
        if choice == 0:
            datatmp = datanc.variables['LW_IN_F_MDS'][:]
            tmp = datanc.variables['LW_OUT'][:]
        else:
            datatmp = datanc.variables['lw_down'][:]
            tmp = datanc.variables['lw_up'][:]
        if len(tmp.shape) == 1:
            for lid in range(len(tmp)):
                datatmp[lid] -= tmp[lid]
        elif len(tmp.shape) == 2:
            for tid in range(len(tmp)):
                for lid in range(len(tmp[0])):
                    datatmp[tid][lid] -= tmp[tid,lid]
        elif len(tmp.shape) == 3:
            for tid in range(len(tmp)):
                for yid in range(len(tmp[0])):
                    for xid in range(len(tmp[0][0])):
                        datatmp[tid][yid][xid] -= tmp[tid,yid,xid]
        else:
            sys.exit("cannot handle variables with dimensions more than 3 (lon, lat, time or vector, time)")
    elif varname.lower() == 'lw_in':
        if choice == 0:
            datatmp = datanc.variables['LW_IN_F_MDS'][:]
        else:
            datatmp = datanc.variables['lw_down'][:]
    elif varname.lower() == 'lw_out':
        if choice == 0:
            datatmp = datanc.variables['LW_OUT'][:]
        else:
            datatmp = datanc.variables['lw_up'][:]
    elif varname.lower() == 'sw_in':
        if choice == 0:
            datatmp = datanc.variables['SW_IN_F_MDS'][:]
        else:
            datatmp = datanc.variables['sw_down'][:]
    elif varname.lower() == 'sw_out':
        if choice == 0:
            datatmp = datanc.variables['SW_OUT'][:]
        else:
            datatmp = datanc.variables['sw_up'][:]
    elif '-' == varname[0] and not '-' in varname[1:]:
        datatmp = datanc.variables[varname[1:]][:]
        if len(datatmp.shape) == 1:
            for lid in range(len(datatmp)):
                datatmp[lid] = -datatmp[lid]
        elif len(datatmp.shape) == 2:
            for tid in range(len(datatmp)):
                for lid in range(len(datatmp[0])):
                    datatmp[tid][lid] = -datatmp[tid][lid]
        elif len(datatmp.shape) == 3:
            for tid in range(len(datatmp)):
                for yid in range(len(datatmp[0])):
                    for xid in range(len(datatmp[0][0])):
                        datatmp[tid][yid][xid] = -datatmp[tid][yid][xid]
        else:
            sys.exit("cannot handle variables with dimensions more than 3 (lon, lat, time or vector, time)")
    elif '-' in varname:
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
            sys.exit("cannot handle variables with dimensions more than 3 (lon, lat, time or vector, time)")
    else:
        datatmp = datanc.variables[varname][:]
    timeFlag = 1
    try:
        timetmp = datanc.variables['time'][:]
    except:
        timeFlag = 0
    dims = len(datatmp.shape)
    tdim = "FIRST"
    if dims == 1:
        dataY = np.concatenate((dataY,datatmp))
        dataX = np.arange(len(dataY))
        timeX = np.arange(len(dataY))
        comment = ''
    elif dims == 2:
        if lpt != '' or ',' not in lonr+latr:
            if tdim.upper() == "LAST":
                tds = len(datatmp) - 1
                lds = len(datatmp[0])
            else:
                tds = len(datatmp[0]) - 1
                lds = len(datatmp)
            if lpt != '':
                lpts = lpt.split(',')
                lpt_index = int(lpts[0])
                lpt_index = max(min(lpt_index,lds-1),0)
            else:
                lpt_index = 0
                if len(lons) > 0 and len(lats) > 0:
                    lpt_index = find_nearest_point(lonr, latr, lons, lats)
            tstepBeg, tstepEnd = tstepS2I(tstep, tds)
            timesteps = tstepEnd - tstepBeg + 1
            if timesteps > 0:
                dataX = np.zeros(timesteps)
                dataY = np.zeros(timesteps)
                timeX = np.array([])
                for tid in range(timesteps):
                    i = tstepBeg + tid
                    cnt = 0.0
                    dataX[tid] = i
                    if timeFlag > 0:
                        timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i]))
                    else:
                        timeX = np.append(timeX,i)
                    if tdim.upper() == "LAST":
                        dataY[tid] = datatmp[i][lpt_index]
                    else:
                        dataY[tid] = datatmp[lpt_index][i]
                if len(lons) > 0 and len(lats) > 0:
                    comment = 'Location: %s, %s; Index: %s' % (lonF2S(lons[lpt_index]), latF2S(lats[lpt_index]), lpt_index)
                else:
                    comment = 'Index: %s' % (lpt_index)
            else:
                sys.exit("no data to plot")
        else:
            listWE = lonr.split(',')
            listSN = latr.split(',')
            lonBeg = lonS2F(listWE[0])
            lonEnd = lonS2F(listWE[1])
            latBeg = latS2F(listSN[0])
            latEnd = latS2F(listSN[1])
            tds = len(datatmp) - 1
            lds = len(datatmp[0])
            tstepBeg, tstepEnd = tstepS2I(tstep, tds)
            timesteps = tstepEnd - tstepBeg + 1
            if timesteps > 0:
                dataX = np.zeros(timesteps)
                dataY = np.zeros(timesteps)
                timeX = np.array([])
                for tid in range(timesteps):
                    i = tstepBeg + tid
                    cnt = 0.0
                    dataX[tid] = i
                    if timeFlag > 0:
                        timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i]))
                    else:
                        timeX = np.append(timeX,tid)
                    for lid in range(lds):
                        if lons[lid] >= lonBeg and lons[lid] <= lonEnd and lats[lid] >= latBeg and lats[lid] <= latEnd:
                            dataY[tid] += datatmp[i][lid]
                            cnt += 1.0
                    if cnt > 0.5:
                        dataY[tid] /= cnt
                    else:
                        dataY[tid] = -999999.99
                comment = 'Domain average: %s - %s, %s - %s' % (lonF2S(lonBeg), lonF2S(lonEnd), latF2S(latBeg), latF2S(latEnd))
            else:
                sys.exit("no data to plot")
    elif dims == 3:
        if lpt != '' or ',' not in lonr+latr:
            tds = len(datatmp) - 1
            yds = len(datatmp[0])
            xds = len(datatmp[0][0])
            if lpt != '':
                lpts = lpt.split(',')
                i_index = int(lpts[0])
                if ',' not in lpt:
                    print("warning: j index is set to 0 in case no j index")
                    j_index = 0
                else:
                    j_index = int(lpts[1])
                i_index = max(min(i_index,xds-1),0)
                j_index = max(min(j_index,yds-1),0)
            else:
                i_index,j_index = find_nearest_point2D(lonr, latr, lons, lats)
            tstepBeg, tstepEnd = tstepS2I(tstep, tds)
            timesteps = tstepEnd - tstepBeg + 1
            if timesteps > 0:
                dataX = np.zeros(timesteps)
                dataY = np.zeros(timesteps)
                timeX = np.array([])
                for tid in range(timesteps):
                    t_index = tstepBeg + tid
                    cnt = 0.0
                    dataX[tid] = t_index
                    if timeFlag > 0:
                        timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[t_index]))
                    else:
                        timeX = np.append(timeX,tid)
                    dataY[tid] = datatmp[t_index][j_index][i_index]
                if len(lons) > 0 and len(lats) > 0:
                    comment = 'Location: %s, %s; Index: %s, %s' % (lonF2S(lons[i_index]), latF2S(lats[j_index]), str(i_index), str(j_index))
                else:
                    comment = 'Index: %s, %s' % (str(i_index), str(j_index))
            else:
                sys.exit("no data to plot")
        else:
            listWE = lonr.split(',')
            listSN = latr.split(',')
            lonBeg = lonS2F(listWE[0])
            lonEnd = lonS2F(listWE[1])
            latBeg = latS2F(listSN[0])
            latEnd = latS2F(listSN[1])
            tds = len(datatmp) - 1
            yds = len(datatmp[0])
            xds = len(datatmp[0][0])
            tstepBeg, tstepEnd = tstepS2I(tstep, tds)
            timesteps = tstepEnd - tstepBeg + 1
            if timesteps > 0:
                dataX = np.zeros(timesteps)
                dataY = np.zeros(timesteps)
                timeX = np.array([])
                for tid in range(timesteps):
                    t_index = tstepBeg + tid
                    cnt = 0.0
                    dataX[tid] = t_index
                    if timeFlag > 0:
                        timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[t_index]))
                    else:
                        timeX = np.append(timeX,tid)
                    for yid in range(yds):
                        if lats[yid] >= latBeg and lats[yid] <= latEnd:
                            for xid in range(xds):
                                if lons[xid] >= lonBeg and lons[xid] <= lonEnd:
                                    dataY[tid] += datatmp[t_index][yid][xid]
                                    cnt += 1.0
                    if cnt > 0.5:
                        dataY[tid] /= cnt
                    else:
                        dataY[tid] = -999999.99
                comment = 'Domain average: %s - %s, %s - %s' % (lonF2S(lonBed), lonF2S(lonEnd), latF2S(latBeg), latF2S(latEnd))
            else:
                sys.exit("no data to plot")
    else:
        sys.exit("cannot handle variables with dimensions more than 2 (vector, time)")
    return dataX, timeX, dataY, comment

def read_var(datapath, geopath, varname, tstep, lonr, latr, lpt):
    obsfiles = glob.glob(datapath)
    geofile  = glob.glob(geopath)
    opath, obsfname = ntpath.split(datapath)
    gpath, geofname = ntpath.split(geopath)
    geofile_new  = geofile
    obsfiles_new = obsfiles
    choice = 0
    if opath == '' and gpath == '':
        cpath = os.getcwd()
        obsfiles_new = glob.glob(os.path.join(cpath,datapath))
        geofile_new = glob.glob(os.path.join(cpath,geofname))
    elif opath != '' and gpath == '':
        geofile_new = np.append(geofile_new, [os.path.join(opath,geofname)])
    elif opath == '' and gpath != '':
        obsfiles_new = glob.glob(os.path.join(gpath,datapath))
    lats = np.array([])
    lons = np.array([])
    dataX = np.array([[]])
    timeX = np.array([])
    dataY = np.array([])
    if not geopath == "":
        for g in geofile_new:
            geonc = nc.Dataset(g)
            if choice == 0:
                lattmp = geonc.variables['latitude'][:]
                lontmp = geonc.variables['longitude'][:]
            else:
                lattmp = geonc.variables['fluxnet_lat'][:]
                lontmp = geonc.variables['fluxnet_lon'][:]
            lats = np.concatenate((lats,lattmp))
            lons = np.concatenate((lons,lontmp))
            geonc.close()
    for f in obsfiles_new:
        datanc = nc.Dataset(f)
        if geopath == "":
            if choice == 0:
                if 'longitude' in datanc.variables.keys() and 'latitude' in datanc.variables.keys():
                    lattmp = datanc.variables['latitude'][:]
                    lontmp = datanc.variables['longitude'][:]
                    lats = np.concatenate((lats,lattmp))
                    lons = np.concatenate((lons,lontmp))
            else:
                if 'fluxnet_lon' in datanc.variables.keys() and 'fluxnet_lat' in datanc.variables.keys():
                    lattmp = datanc.variables['fluxnet_lat'][:]
                    lontmp = datanc.variables['fluxnet_lon'][:]
                    lats = np.concatenate((lats,lattmp))
                    lons = np.concatenate((lons,lontmp))
        varnames = varname.split(',')
        for vid in range(len(varnames)):
            X, T, Y, comment = read_single_var(datanc, varnames[vid], tstep, lonr, latr, lpt, lons, lats)
            if vid == 0:
                dataX = np.array([X])
                timeX = np.array([T])
                dataY = np.array([Y])
            else:
                dataX = np.insert(dataX,len(dataX),X,0)
                timeX = np.insert(timeX,len(timeX),T,0)
                dataY = np.insert(dataY,len(dataY),Y,0)
        datanc.close()
    return dataX, timeX, dataY, varnames, comment

#def gen_figure(inpath, geopath, outpath, screen, tstep, lonr, latr, lpt, colors, styles):
def gen_figure(inpath, geopath, outpath, screen, tstep, lonr, latr, lpt, varname, colors, styles):
#   varname = "swnet,lwdown,lwup,hflx,evap,gflux,snohf"
#   varname = "Rnet,hflx,evap,gflux,snohf"
#   varname = "Rnet,-hflx,-evap,gflux,snohf"
    dataX, timeX, dataY, varnames, comment = read_var(inpath, geopath, varname, tstep, lonr, latr, lpt)
    plotpath = outpath+'/timeseries_balance.png'
    metadata = {
                'var': varname
                }
    plot_time_series(dataX, timeX, dataY, metadata, plotpath, screen, varnames, comment, colors, styles)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-g', '--geo', help="path to the geographic info file", default="")
    ap.add_argument('-s', '--screen', help="no if plot to file", default="yes")
    ap.add_argument('-t', '--tstep', help="time step for plotting", default="0,-1")
    ap.add_argument('-p', '--point', help="location index for plotting", default="")
#   ap.add_argument('-v', '--variable', help="variable name to plot", default="rnet,le,hf,gf,lw_in,lw_out,sw_in,sw_out")
    ap.add_argument('-v', '--variable', help="variable name to plot", default="rnet,le,hf,gf,lwnet,swnet")
    ap.add_argument('-c', '--color', help="color for lines", default="")
    ap.add_argument('-l', '--linestyle', help="line styles", default="")
    ap.add_argument('-lon', '--longitude', help="longitude range for plotting", default="-180,180")
    ap.add_argument('-lat', '--latitude', help="latitude range for plotting", default="-90,90")
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.screen, MyArgs.tstep, MyArgs.longitude, MyArgs.latitude, MyArgs.point, MyArgs.variable, MyArgs.color, MyArgs.linestyle)
#   gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.screen, MyArgs.tstep, MyArgs.longitude, MyArgs.latitude, MyArgs.point, MyArgs.color, MyArgs.linestyle)
