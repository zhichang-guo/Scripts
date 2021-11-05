#!/usr/bin/python
######################################################################################
## Compare two directories, examples:
##   1. python compare_two_directories.py dirA DirB false
##      compare recursively two directories and find those different files
##   2. python compare_two_directories.py dirA DirB true
##      compare cursively two directories, find those different files, and compare 
##      each file
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
######################################################################################
import os
import sys
import getopt
import glob
from filecmp import dircmp

def print_diff_files(dcmp, flag, dirA, dirB, ldir):
    count = 0
    if len(dcmp.diff_files) > 0:
        print("*********************************************************************")
        print("Right directory: "+dirA)
        print("Left  directory: "+dirB)
        print("Local directory: "+ldir)
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
        if flag.upper() == 'VV':
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
    dirA = sys.argv[1]
    dirB = sys.argv[2]
    option = sys.argv[3]
#
    strlenA = len(dirA)
    strlenB = len(dirB)
    for subdir, dirs, files in os.walk(dirA):
        for dir in dirs:
            realpathA = subdir + os.sep + dir
            localpathA = realpathA[strlenA+1:]
            realpathB = dirB + os.sep + localpathA
            flag = 1
            if len(localpathA) > 1 and localpathA.startswith('.') and not localpathA.startswith('.'+os.sep):
                flag = 0
            if os.sep+'.' in localpathA:
                flag = 0
            if len(dir) > 1 and dir.startswith('.') and not dir.startswith('.'+os.sep):
                flag = 0
            if flag == 1:
                if os.path.isdir(realpathB):
                    dcmp = dircmp(realpathA, realpathB)
                    print_diff_files(dcmp, option, dirA, dirB, localpathA)
                else:
                    print("The directory "+localpathA+" exists in "+dirA+" but not in "+dirB)
            else:
                break

    for subdir, dirs, files in os.walk(dirB):
        for dir in dirs:
            realpathB = subdir + os.sep + dir
            localpathB = realpathA[strlenB+1:]
            realpathA = dirA + os.sep + localpathB
            flag = 1 
            if len(localpathB) > 1 and localpathB.startswith('.') and not localpathB.startswith('.'+os.sep):
                flag = 0 
            if os.sep+'.' in localpathB:
                flag = 0
            if len(dir) > 1 and dir.startswith('.') and not dir.startswith('.'+os.sep):
                flag = 0
            if flag == 1:
                if not os.path.isdir(realpathA):
                    print("The directory "+localpathB+" exists in "+dirB+" but not in "+dirA)
            else:
                break
if __name__ == "__main__":
    compare_dirs(sys.argv[1:])
