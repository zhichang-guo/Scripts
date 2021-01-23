#!/usr/bin/env python3
#####################################################################################################
## Add new variables into an exising netcdf file with constant values, examples:
##   1. python add_variable_nc.py -i file.nc -v temperatureAir@ObsError,temperatureAir@PreQC -d nlocs
##      two vaiables are added into the file with dimension nlocs and default value zero
##   2. python add_variable_nc.py -i file.nc -v a,b -d nlocs -u m,K -c 0.1,300
##      two vaiables are added into the file with units m and K, values 0.1 and 300 respectively
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
###################################################################################################
from netCDF4 import Dataset
import argparse
import numpy as np
import sys

def add_variable(inpath, varname, dname, group, unit, constant):
    rootgrp = Dataset(inpath, "a")
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
    for vid in range(len(varnames)):
        var = newgrp.createVariable(varnames[vid], "f4", (dname))
        unit = ''
        if vid < len(units) and not units[vid] == '':
            var.units = units[vid]
        const = ''
        if vid < len(constants) and not constants[vid] == '':
            const = constants[vid]
        data = np.zeros(dsize)
        if not const == '':
            data[:] = float(const)
        var[:] = data[:]    
    rootgrp.close()
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-v', '--variable', help="variable name to add", required=True)
    ap.add_argument('-d', '--dimension', help="dimension name for the variable", required=True)
    ap.add_argument('-g', '--group',  help="group name", default="")
    ap.add_argument('-u', '--unit', help="units", default='')
    ap.add_argument('-c', '--constant', help="initial constant value for the variable", default='')
    MyArgs = ap.parse_args()
    add_variable(MyArgs.input, MyArgs.variable, MyArgs.dimension, MyArgs.group, MyArgs.unit, MyArgs.constant)
