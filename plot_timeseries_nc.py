#!/usr/bin/env python3
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
    # plot generic world map
    defaultColors = np.array(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    if screen.upper() == "NO":
        matplotlib.use('agg')
    fig = plt.figure(figsize=(12,8))
    if comment == '':
        plttitle = 'Variable: %s' % (metadata['var'])
    else:
        plttitle = comment
    colors = color.split(',')
    styles = style.split(',')
    for i in range(len(dataY)):
        time = timeX[i]
        x = dataX[i]
        y = dataY[i]
        c = defaultColors[min(i,len(defaultColors)-1)]
        s = 'solid'
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
        lpt_index = index
    else:
        print(strLon+" "+strLat+" "+str(dist)+" "+str(len(lons)))
        sys.exit("the nearest point is not found")
    return lpt_index

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
    dataX = np.array([])
    timeX = np.array([])
    dataY = np.array([])
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
            sys.exit("cannot handle variables with dimensions more than 3 (lon, lat, time or vector, time)")
    else:
        datatmp = datanc.variables[varname][:]
    timetmp = datanc.variables['time'][:]
    dims = len(datatmp.shape)
    if dims == 1:
        dataY = np.concatenate((dataY,datatmp))
        dataX = np.arange(len(dataY))
        timeX = np.arange(len(dataY))
        comment = ''
    elif dims == 2:
        if lpt != '' or ',' not in lonr+latr:
            tds = len(datatmp) - 1
            lds = len(datatmp[0])
            if lpt != '':
                lpts = lpt.split(',')
                lpt_index = int(lpts[0])
                lpt_index = max(min(lpt_index,lds-1),0)
            else:
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
                    #timeX[tid] = datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i])
                    timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i]))
                    dataY[tid] = datatmp[i][lpt_index]
                comment = 'Location: %s, %s; Index: %s' % (lonF2S(lons[lpt_index]), latF2S(lats[lpt_index]), lpt_index)
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
                    #timeX[tid] = datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i])
                    timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i]))
                    for lid in range(lds):
                        if lons[lid] >= lonBeg and lons[lid] <= lonEnd and lats[lid] >= latBeg and lats[lid] <= latEnd:
                            dataY[tid] += datatmp[i][lid]
                            cnt += 1.0
                    if cnt > 0.5:
                        dataY[tid] /= cnt
                    else:
                        dataY[tid] = -999999.99
                comment = 'Domain average: %s - %s, %s - %s' % (lonF2S(lonBed), lonF2S(lonEnd), latF2S(latBeg), latF2S(latEnd))
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
                    #timeX[tid] = datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i])
                    timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[t_index]))
                    dataY[tid] = datatmp[t_index][j_index][i_index]
                comment = 'Location: %s, %s; Index: %s, %s' % (lonF2S(lons[i_index]), latF2S(lats[j_index]), str(i_index), str(j_index))
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
                    #timeX[tid] = datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[i])
                    timeX = np.append(timeX,datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0) + timedelta(seconds=timetmp[t_index]))
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

def read_var(datapath, geopath, varname, tstep, lonr, latr, lpt, fact):
    obsfiles = glob.glob(datapath)
    geofile  = glob.glob(geopath)
    opath, obsfname = ntpath.split(datapath)
    gpath, geofname = ntpath.split(geopath)
    geofile_new  = geofile
    obsfiles_new = obsfiles
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
            lattmp = geonc.variables['latitude'][:]
            lontmp = geonc.variables['longitude'][:]
            lats = np.concatenate((lats,lattmp))
            lons = np.concatenate((lons,lontmp))
            geonc.close()
    for f in obsfiles_new:
        datanc = nc.Dataset(f)
        if geopath == "":
            lattmp = datanc.variables['latitude'][:]
            lontmp = datanc.variables['longitude'][:]
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
        if not fact == '1':
            factor = float(fact)
            for vid in range(len(dataY)):
                for tid in range(len(dataY[0])):
                    dataY[vid][tid] *= factor
            comment += ', Factor: %s'%(fact)
    return dataX, timeX, dataY, varnames, comment

def gen_figure(inpath, geopath, outpath, varname, screen, tstep, lonr, latr, lpt, fact, colors, styles):
    # read the files to get the 2D array to plot
    dataX, timeX, dataY, varnames, comment = read_var(inpath, geopath, varname, tstep, lonr, latr, lpt, fact)
    plotpath = outpath+'/timeseries_%s.png' % (varname)
    metadata = {
                'var': varname
                }
    plot_time_series(dataX, timeX, dataY, metadata, plotpath, screen, varnames, comment, colors, styles)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-g', '--geo', help="path to the geographic info file", default="")
    ap.add_argument('-v', '--variable', help="variable name to plot", required=True)
    ap.add_argument('-s', '--screen', help="no if plot to file", default="yes")
    ap.add_argument('-t', '--tstep', help="time step for plotting", default=0)
    ap.add_argument('-p', '--point', help="location index for plotting", default="")
    ap.add_argument('-c', '--color', help="color for lines", default="")
    ap.add_argument('-l', '--linestyle', help="line styles", default="")
    ap.add_argument('-lon', '--longitude', help="longitude range for plotting", default="-180,180")
    ap.add_argument('-lat', '--latitude', help="latitude range for plotting", default="-90,90")
    ap.add_argument('-factor', '--fact', help="factor for units conversion", default="1")
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.variable, MyArgs.screen, MyArgs.tstep, MyArgs.longitude, MyArgs.latitude, MyArgs.point, MyArgs.fact,MyArgs.color,MyArgs.linestyle)
