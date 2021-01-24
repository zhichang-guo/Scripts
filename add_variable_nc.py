#!/usr/bin/env python3
######################################################################################################
## Add new variables to an exising netcdf file with constant values, examples:
##   1. python add_variable_nc.py -i file.nc -v temperatureAir@ObsError,temperatureAir@PreQC -d nlocs
##      Two vaiables temperatureAir@ObsError and temperatureAir@PreQC are added to the NetCDF file 
##      file.nc with dimension nlocs and default value zero
##   2. python add_variable_nc.py -i file.nc -v a,b -d nlocs -u m,K -c 0.1,300
##      Two vaiables a and b are added to the file with units m and K, values 0.1 and 300 respectively
##   3. python add_variable_nc.py -i in.nc -o out.nc -v a,b -d nlocs -u m,K -c 0.1,300
##      The input file (in.nc) is duplicated to a new file (out.nc) and two vaiables a,b are added
##      to the new file with units m and K, values 0.1 and 300 respectively
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
######################################################################################################
from netCDF4 import Dataset
import argparse
import numpy as np
import sys
import os
import shutil

def add_variable(inpath, outpath, varname, dname, group, unit, constant):
    if not os.path.isfile(inpath):
        sys.exit("The input file \""+inpath+"\" is not found")
    if not outpath == '' and not outpath == inpath:
        if os.path.isfile(outpath):
            print("\""+ outpath + "\" already exists. Input 'y' if you want to replace it, otherwise input 'n'")
            answer = input("Input your choice: \n")
            if not 'Y' in answer.upper():
                sys.exit("Action Canceled")
        filepath = outpath
        shutil.copy(inpath, filepath)
        print("The file "+inpath+" is copied to "+outpath)
    else:
        filepath = inpath
    rootgrp = Dataset(filepath, "a")
    if not dname in rootgrp.dimensions:
        sys.exit("Error: "+dname+" is not in dimension list")
    dobject = rootgrp.dimensions[dname]
    dsize = len(dobject)
    newgrp = rootgrp
    if not group == '':
        if group in rootgrp.groups:
            newgrp = rootgrp.groups[group]
        else:
            newgrp = rootgrp.createGroup(group)
    varnames = varname.split(',')
    units = unit.split(',')
    constants = constant.split(',')
    print("The following variable will be added to the file "+filepath)
    for vid in range(len(varnames)):
        if varnames[vid] in newgrp.variables.keys():
            print("***************************************************************************************")
            print("Error: Failed to insert the variable "+varnames[vid]+" since it is already in the file")
            print("***************************************************************************************")
            sys.exit()
        var = newgrp.createVariable(varnames[vid], "f4", (dname))
        comment = "    name: "+varnames[vid]+", dimension: "+dname
        unit = ''
        if vid < len(units) and not units[vid] == '':
            comment += ", units: "+units[vid]
            var.units = units[vid]
        const = ''
        if vid < len(constants) and not constants[vid] == '':
            comment += ", value: "+constants[vid]
            const = constants[vid]
        else:
            comment += ", value: 0.0"
        data = np.zeros(dsize)
        if not const == '':
            data[:] = float(const)
        var[:] = data[:]    
        print(comment)
    rootgrp.close()
    print("The script ended normally!")
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-o', '--output', help="path to the output file", default='')
    ap.add_argument('-v', '--variable', help="variable name to add", required=True)
    ap.add_argument('-d', '--dimension', help="dimension name for the variable", required=True)
    ap.add_argument('-g', '--group',  help="group name", default='')
    ap.add_argument('-u', '--unit', help="units", default='')
    ap.add_argument('-c', '--constant', help="initial constant value for the variable", default='')
    MyArgs = ap.parse_args()
    add_variable(MyArgs.input, MyArgs.output, MyArgs.variable, MyArgs.dimension, MyArgs.group, MyArgs.unit, MyArgs.constant)
