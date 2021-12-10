#!/usr/bin/env python3
######################################################################################################
## Rename variables in an exising netcdf file, examples:
##   1. python change_varatt_value_nc.py -i file.nc -v elevation -a checksum -n 999
##      The attribute checksum for the variable elevation is reset to 999 in the file "file.nc"
##   2. python change_varatt_value_nc.py -i in.nc -o out.nc -v elevation -a checksum -n 999
##      The input file "in.nc" is duplicated to a new file "out.nc" and the attribute checksum for
##      the variable elevation is reset to 999 in that file.
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
######################################################################################################
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-o', '--output', help="path to the output file", default='')
    ap.add_argument('-v', '--variable', help="variable name", required=True)
    ap.add_argument('-a', '--attname', help="variable attribute name", required=True)
    ap.add_argument('-n', '--attvalue', help="new variable attribute value", required=True)
from netCDF4 import Dataset
import argparse
import numpy as np
import sys
import os
import shutil

def change_varatt_value(inpath, outpath, varname, attname, attvalue):
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
    print("The attribute \""+attname+"\" for the variable \""+varname+"\" in the file \""+ filepath + "\" will be reset")
    if not varname in rootgrp.variables.keys():
        sys.exit("The variable \"" + varname + "\" is not found in the file \""+filepath+"\"")
    else:
        att_src = rootgrp.variables[varname]
        old_value = att_src.getncattr(attname)
        att_src.setncattr(attname,attvalue)
        print("The attribute \""+attname+": "+old_value+"\" for the variable \""+varname+"\" was reset to "+attvalue)
    rootgrp.close()
    print("The script ended normally!")
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-o', '--output', help="path to the output file", default='')
    ap.add_argument('-v', '--variable', help="variable name", required=True)
    ap.add_argument('-a', '--attname', help="variable attribute name", required=True)
    ap.add_argument('-n', '--attvalue', help="new variable attribute value", required=True)
    MyArgs = ap.parse_args()
    change_varatt_value(MyArgs.input, MyArgs.output, MyArgs.variable, MyArgs.attname, MyArgs.attvalue)
