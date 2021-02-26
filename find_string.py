#!/usr/bin/python
######################################################################################
## Search for a string in all files under a directory, examples:
##   1. python find_string.py string rootDir false
##      Search for "string" in all files under the directory "rootDir"
##   2. python find_string.py string rootDir true
##      Search for "string" in all files under the directory "rootDir" regardless of 
##      the case of letters
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov
######################################################################################
import os
import sys

def find_string(argv):
    # total arguments
    num = len(sys.argv)
    if not num == 4:
        print("*********************************************************************")
        print("*                                                                   *")
        print("*      Usage: python "+sys.argv[0] + " String RootDirectory Flag       *")
        print("*                                                                   *")
        print("*********************************************************************")
        sys.exit()
    string = sys.argv[1]
    rootdir = sys.argv[2]
    flag = sys.argv[3]
    for folder, dirs, files in os.walk(rootdir):
        for file in files:
            if not file.startswith('.'):
                fullpath = os.path.join(folder, file)
                if flag.upper() == 'TRUE':
                    cmd = "grep -i \'"+str(string)+"\' \'"+str(fullpath) + "\'"
                else:
                    cmd = "grep \'"+str(string)+"\' "+str(fullpath) + "\'"
                result = os.popen(cmd).read()
                if result:
                    print("------------------------------------------------")
                    print("Found the string \""+str(string)+"\" in the file \""+str(fullpath)+"\"")
                    print(result)
                #try:
                #    with open(fullpath, "r") as f:
                #        for line in f:
                #            #if string in line:
                #            if string in line:
                #                print("file: ", ' ', line, ' ', fullpath)
                #except UnicodeDecodeError:
                #    pass # Fond non-text data
if __name__ == "__main__":
    find_string(sys.argv[1:])
