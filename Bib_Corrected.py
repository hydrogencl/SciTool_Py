#!/usr/bin/python
import re
class BIB_CORRECTING_TOOL:
    DIR_IN   = "."
    FILE_IN  = "input.bib"
    FILE_OUT = "output.bib"
    ARR_OUT = []

    def __init__ (self, FILE_IN, FILE_OUT="out.bib", DIR_IN="./", DIR_OUT="./"):
        self.FILE_BIB_IN  = open("{0:s}/{1:s}".format(DIR_IN, FILE_IN ), "r")
        self.FILE_BIB_OUT = open("{0:s}/{1:s}".format(DIR_OUT, FILE_OUT), "w")
        
    def REMOVING_ROUND_BRACKET(self):
        ARR_ALL_IN_LINE = self.FILE_BIB_IN.readlines()
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
            self.ARR_OUT.append(line)

    def WRITEOUT_BIB(self):
        for line in self.ARR_OUT:
            self.FILE_BIB_OUT.write("{0:s}".format( line ))
if __name__ == "__main__":
    print("PLEASE state at least the input .bib file name:")
    try:
        STR_FILE_IN = raw_input("filename:")
    except:
        STR_FILE_IN = input("filename:")
    BCT = BIB_CORRECTING_TOOL(STR_FILE_IN)
    BCT.REMOVING_ROUND_BRACKET()
    BCT.WRITEOUT_BIB()
    print("Output as: {0:s}".format(BCT.FILE_OUT)) 
