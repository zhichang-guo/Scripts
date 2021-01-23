#!/usr/bin/env python3
######################################################################################################
## Add new variables into an exising netcdf file with constant values, examples:
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

def rename_variable(inpath, outpath, varname, newname):
    if not os.path.isfile(inpath):
        sys.exit("The input file \""+inpath+"\" is not found")
    if not outpath == '' and not outpath == inpath:
        if os.path.isfile(outpath):
            print("\""+ outpath + "\" already exists. Input 'Y' if you want to replace it, otherwise input 'N'")
            answer = input("Input your choice: \n")
            if not 'Y' in answer.upper():
                sys.exit("Action Canceled")
        filepath = outpath
        shutil.copy(inpath, filepath)
        print("The file "+inpath+" is copied to "+outpath)
    else:
        filepath = inpath
    rootgrp = Dataset(filepath, "a")
    varnames = varname.split(',')
    newnames = newname.split(',')
    if not len(varnames) == len(newnames):
        sys.exit("The numbers of old names and new names do not match!")
    print("The following variable in the file \""+ filepath + "\" will be renamed")
    for vid in range(len(varnames)):
        if not varnames[vid] in rootgrp.variables.keys():
            sys.exit("The variable \"" + varnames[vid] + "\" is not found in the file \""+filepath+"\"")
        if newnames[vid] in rootgrp.variables.keys():
            sys.exit("The target variable name \"" + newnames[vid] + "\" already exists in the file \""+filepath+"\"")
        rootgrp.renameVariable(varnames[vid],newnames[vid])
        print("    " + varnames[vid] + " ---> " + newnames[vid])
    rootgrp.close()
    print("The script ended normally!")
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-o', '--output', help="path to the output file", default='')
    ap.add_argument('-v', '--variable', help="old variable name", required=True)
    ap.add_argument('-n', '--new', help="new variable name", required=True)
    MyArgs = ap.parse_args()
    rename_variable(MyArgs.input, MyArgs.output, MyArgs.variable, MyArgs.new)
