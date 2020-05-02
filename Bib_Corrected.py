#!/usr/bin/python
import re
class BIB_CORRECTING_TOOL:
    DIR_IN   = "./"
    FILE_IN  = "input.bib"
    FILE_OUT = "output.bib"

    def __init__ (self, DIR_IN, FILE_IN, FILE_OUT):
        FILE_BIB_IN  = open("{0:s}/{1:s}".format(DIR_IN, FILE_IN ), "r")
        FILE_BIB_OUT = open("{0:s}/{1:s}".format(DIR_IN, FILE_OUT), "w")
        ARR_ALL_IN_LINE = FILE_BIB_IN.readlines()
        
    def REMOVING_ROUND_BRACKET(self):
        for line in ARR_ALL_IN_LINE:    
            a = re.match("@",line)
            if a == None:
                if line[0] != "}":
                    line = line.replace(" = ",' = "')
                    line = line.replace("{","")
                    line = line.replace("}","")
                    if len(line) > 3:
                        if line[-2] == ",":
                            line = line.replace(',\n', '",\n')
                        else:
                            line = line.replace("\n",'"\n')
                    line = line.replace("\\textbackslash", "\\{textbackslash}")

            FILE_BIB_OUT.write("{0:s}".format( line ))
