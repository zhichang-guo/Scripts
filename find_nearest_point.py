#######################################################################################
## Find the nearest grid point index, examples:
##   1. For vector data, longitude/latitude are assumed to be functions of grid point.
##      This code will find the nearest grid point index with given lon/lat. 
##      python find_nearest_point.py -i file_name.nc -l 60W,20S
##   2. For field data, lon/lat are assumed to be functions of i/j index respectively.
##      This code will find the nearest i/j index with given lon/lat.
##      python find_nearest_point.py -i file_name.nc -l 60W,20S -f field
##   3. Find the nearest grid point only over land/ocean, land-sea mask is required. 
##      python find_nearest_point.py -i file_name.nc -l 60W,20S -g land           
##   4. The default variable names are assumed to be longitude/latitude/lsmask. 
##      Otherwise, specify them with -v
##      python find_nearest_point.py -i file_name.nc -l 60W,20S -v lon/lat/mask
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
#######################################################################################
import netCDF4 as nc
import argparse
import numpy as np
import glob
import os
import sys

def read_var(geopath, gol, varname):
    default_varnames = ["longitude", "latitude", "lsmask"]
    geofile = glob.glob(geopath)
    lats = np.array([])
    lons = np.array([])
    mask = np.array([])
    vnames = default_varnames
    varnames = varname.split(',')
    for v in range(len(varnames)):
        vnames[v] = varnames[v]
    geonc = nc.Dataset(geopath)
    if vnames[0] in geonc.variables.keys():
        lontmp = geonc.variables[vnames[0]][:]
        lons = np.concatenate((lons,lontmp))
    else:
        sys.exit("Error: variable "+vnames[0]+" cannot be found")
    if vnames[1] in geonc.variables.keys():
        lattmp = geonc.variables[vnames[1]][:]
        lats = np.concatenate((lats,lattmp))
    else:
        sys.exit("Error: variable "+vnames[1]+" cannot be found")
    if vnames[2] in geonc.variables.keys():
        masktmp = geonc.variables[vnames[2]][:]
        mask = np.concatenate((mask,masktmp))
    elif 'LAND' in gol.upper() or 'OCEAN' in gol.upper():
        print("Warning: land-sea mask cannot be found, all grid points are candidates")
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

def find_nearest_point(gol, strLocation, lons, lats, mask):
    xlon, xlat = lonlatS2F(strLocation)
    if xlon < 0.0:
        xlon += 360.
    xlona = xlon - 360.
    xlonb = xlon + 360.
    distanceMin = 9999.9
    lpt_index = -1
    found = -1
    if (not len(lons) == len(lats)) or ('ALL' not in gol.upper() and not len(lons) == len(mask)):
        sys.exit("Error: length of longitude/latitude/mask differs")
    for lid in range(len(lons)):
        if gol.upper() == 'ALL' or ('LAND' in gol.upper() and mask[lid] > 0) or ('OCEAN' in gol.upper() and mask[lid] < 1):
            if lons[lid] < 0.0:
                lons[lid] += 360.
            distLat = (xlat-lats[lid])*(xlat-lats[lid])
            dist  = (xlon-lons[lid])*(xlon-lons[lid]) + distLat
            distA = (xlona-lons[lid])*(xlona-lons[lid]) + distLat
            distB = (xlonb-lons[lid])*(xlonb-lons[lid]) + distLat
            if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
                distanceMin = min(dist,min(distA,distB))
                lpt_index = lid
                found = 1
    if found > 0:
        return lpt_index
    else:
        sys.exit("the nearest point is not found")

def find_nearest_point2D(gol, strLocation, lons, lats, mask):
    xlon, xlat = lonlatS2F(strLocation)
    if xlon < 0.0:
        xlon += 360.
    xlona = xlon - 360.
    xlonb = xlon + 360.
    distanceMin = 9999.9
    x_index = -1
    y_index = -1
    found = -1
    for xid in range(len(lons)):
        for yid in range(len(lats)):
            if gol.upper() == 'ALL' or ('LAND' in gol.upper() and mask[xid,yid] > 0) or ('OCEAN' in gol.upper() and mask[xid,yid] < 1):
                if lons[xid] < 0.0:
                    lons[xid] += 360.
                distLat = (xlat-lats[yid])*(xlat-lats[yid])
                dist  = (xlon-lons[xid])*(xlon-lons[xid]) + distLat
                distA = (xlona-lons[xid])*(xlona-lons[xid]) + distLat
                distB = (xlonb-lons[xid])*(xlonb-lons[xid]) + distLat
                if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
                    distanceMin = min(dist,min(distA,distB))
                    x_index = xid
                    y_index = yid
                    found = 1
    if found > 0:
        return x_index, y_index
    else:
        sys.exit("the nearest point is not found")

def read_and_find(inpath, gol, fov, varname, location):
    lons, lats, mask = read_var(inpath, gol, varname)
    if 'VECTOR' in fov.upper():
        lpt_index = find_nearest_point(gol, location, lons, lats, mask)
        print("|-------------------------------------------------------")
        print("| The nearest point index for location ("+location+") is:      "+str(lpt_index))
        print("| Its longitude/latitude in the file:  ("+lonlatF2S(lons[lpt_index],lats[lpt_index])+")")
        print("|-------------------------------------------------------")
    else:
        i_index, j_index = find_nearest_point2D(gol, location, lons, lats, mask)
        print("|----------------------------------------------------------")
        print("| The nearest point i/j index for location("+location+"): "+str(i_index)+"/"+str(j_index))
        print("| Its longitude/latitude in the file: ("+lonlatF2S(lons[i_index],lats[j_index])+")")
        print("|----------------------------------------------------------")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input',    help="path to the input file", required=True)
    ap.add_argument('-g', '--gol',      help="over all grids or only land", default='all')
    ap.add_argument('-f', '--vector',   help="field or vector", default='vector')
    ap.add_argument('-v', '--variable', help="variable names for longitude/latitude", default="longitude,latitude,lsmask")
    ap.add_argument('-l', '--location', help="longitude/latitude of the location", required=True)
    MyArgs = ap.parse_args()
    read_and_find(MyArgs.input, MyArgs.gol, MyArgs.vector, MyArgs.variable, MyArgs.location)
