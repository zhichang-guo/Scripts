import netCDF4 as nc
import argparse
import numpy as np
import glob
import os
import sys

def read_var(geopath, option, varname):
    geofile = glob.glob(geopath)
    lats = np.array([])
    lons = np.array([])
    mask = np.array([])
    varnames = varname.split(',')
    geonc = nc.Dataset(geopath)
    lontmp = geonc.variables[varnames[0]][:]
    lattmp = geonc.variables[varnames[1]][:]
    lats = np.concatenate((lats,lattmp))
    lons = np.concatenate((lons,lontmp))
    if 'LAND' in option.upper():
        masktmp =  geonc.variables[varnames[2]][:]
        mask = np.concatenate((mask,masktmp))
    return lons, lats, mask
    
def lonlatS2F(strLocation):
    if ',' not in strLocation:
        sys.exit("Error: location is not in the format of longitude,latitude")
    else:
        location = strLocation.upper()
        location.replace('E','')
        location.replace('N','')
        locations = location.split(',')
        if 'W' in locations[0]:
            flon = -1.0*float(locations[0].replace('W',''))
        else:
            flon = float(locations[0])
        if 'S' in locations[1]:
            flat = -1.0*float(locations[1].replace('S',''))
        else:
            flat = float(locations[1])
        return flon,flat

def f2s(value):
    return str(int(value*100)/100.)
    
def lonlatF2S(flon,flat):
    flon = int(flon*100)/100.
    if flon > 180:
        flon = 360. - flon
        flon = int(flon*100)/100.
        strLon = str(flon) + 'W'
    elif flon < 0.0:
        strLon = str(-flon) + 'W'
    else:
        strLon = str(flon) + 'E'
    flat = int(flat*100)/100.
    if flat < 0.0:
        strLat = str(-flat) + 'S'
    else:
        strLat = str(flat) + 'N'
    return strLon+","+strLat

def find_nearest_point(option, strLocation, lons, lats, mask):
    xlon, xlat = lonlatS2F(strLocation)
    if xlon < 0.0:
        xlon += 360.
    distanceMin = 9999.9
    x_index = -1
    y_index = -1
    found = -1
    dims = len(mask.shape)
    if dims == 1:
        if not len(lons) == len(lats):
            sys.exit("Error: length of longitude and latitude differs")
        for lid in range(len(lons)):
            if option.upper() == 'ALL' or ('LAND' in option.upper() and mask[lid] > 0):
                if lons[lid] < 0.0:
                    lons[lid] += 360.
                distance = (xlon-lons[lid])*(xlon-lons[lid]) + (xlat-lats[lid])*(xlat-lats[lid])
                if distance <= distanceMin:
                    distanceMin = distance
                    x_index = lid
                    found = 1
#               print(str(lid)+": "+lonlatF2S(xlon,xlat)+" "+lonlatF2S(lons[lid],lats[lid])+" "+f2s(distance)+" "+f2s(distanceMin))
    elif dims == 2:
        for xid in range(len(lons)):
            for yid in range(len(lats)):
                if option.upper() == 'ALL' or ('LAND' in option.upper() and mask[xid,yid] > 0):
                    if lons[xid] < 0.0:
                        lons[xid] += 360.
                    distance = (xlon-lons[xid])*(xlon-lons[xid]) + (xlat-lats[yid])*(xlat-lats[yid])
                    if distance <= distanceMin:
                        distanceMin = distance
                        x_index = xid
                        y_index = yid
                        found = 1
    else:
        sys.exit("something is wrong with the dimensions")
    if found > 0:
        return x_index, y_index
    else:
        sys.exit("the nearest point is not found")

def read_and_find(inpath, option, varname, location):
    a = 10.
    lons, lats, mask = read_var(inpath, option, varname)
    i_index, j_index = find_nearest_point(option, location, lons, lats, mask)
    if j_index < 0:
        print("|-------------------------------------------------------")
        print("| The nearest point index for location ("+location+") is:      "+str(i_index))
        print("| Its longitude/latitude in the file:  ("+lonlatF2S(lons[i_index],lats[i_index])+")")
        print("|-------------------------------------------------------")
    else:
        print("|----------------------------------------------------------")
        print("| The nearest point i/j index for location("+location+"): "+str(i_index)+"/"+str(j_index))
        print("| Its longitude/latitude in the file: ("+lonlatF2S(lons[i_index],lats[j_index])+")")
        print("|----------------------------------------------------------")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input',    help="path to the input file", required=True)
    ap.add_argument('-g', '--gol',      help="over global or land", default='all')
    ap.add_argument('-v', '--variable', help="variable names for longitude/latitude", default="longitude,latitude,lsmask")
    ap.add_argument('-l', '--location', help="longitude/latitude of the location", required=True)
    MyArgs = ap.parse_args()
    read_and_find(MyArgs.input, MyArgs.gol, MyArgs.variable, MyArgs.location)
