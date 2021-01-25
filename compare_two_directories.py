#!/usr/bin/python
######################################################################################
## Compare two directories, examples:
##   1. python compare_two_directories.py dirA DirB false
##      compare two directories and find those different files
##   2. python compare_two_directories.py dirA DirB true
##      compare two directories, find those different files, and compare each file
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
######################################################################################
import os
import sys
import getopt
import glob
from filecmp import dircmp

def print_diff_files(dcmp, flag):
    count = 0
    if len(dcmp.diff_files) > 0:
        print("*********************************************************************")
        print("Following files are found in both directories and they are different:")
        for name in dcmp.diff_files:
            count += 1
            print("-------------------")
            print("    "+str(count) + ": "+name)
            if flag.upper() == 'TRUE':
                fileA = glob.glob(os.path.join(dcmp.left,name))
                fileB = glob.glob(os.path.join(dcmp.right,name))
                cmd = "diff "+str(fileA[0])+" "+str(fileB[0])
                print("Command: "+cmd)
                result = os.popen(cmd).read()
                print(result)
    if len(dcmp.common_files) > 0:
        print("*********************************************************************")
        print("Following files are found in both directories and they are identical:")
        count = 0
        for name in dcmp.common_files:
            count += 1
            print("    "+str(count) + ": "+name)
    if len(dcmp.left_only) > 0:
        print("*********************************************************************")
        print("Following files are only found in the left directory "+dcmp.left+" :")
        count = 0
        for name in dcmp.left_only:
            count += 1
            print("    "+str(count) + ": "+name)
    if len(dcmp.right_only) > 0:
        print("*********************************************************************")
        print("Following files are only found in the right directory "+dcmp.right+" :")
        count = 0
        for name in dcmp.right_only:
            count += 1
            print("    "+str(count) + ": "+name)
def compare_dirs(argv):
    # total arguments
    num = len(sys.argv)
    if not num == 4:
        print("*********************************************************************")
        print("*                                                                   *")
        print("*    Usage: python "+sys.argv[0] + " Dir_One Dir_Two Flag  *")
        print("*                                                                   *")
        print("*********************************************************************")
        sys.exit()
    dcmp = dircmp(sys.argv[1], sys.argv[2]) 
    print_diff_files(dcmp, sys.argv[3]) 
if __name__ == "__main__":
    compare_dirs(sys.argv[1:])

