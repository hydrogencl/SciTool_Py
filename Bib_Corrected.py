import re

DIR_IN   = "./"
FILE_IN  = "MYREF2.bib"
FILE_OUT = "MYREF2_Corrected.bib"

FILE_BIB_IN  = open("{0:s}/{1:s}".format(DIR_IN, FILE_IN ), "r")
FILE_BIB_OUT = open("{0:s}/{1:s}".format(DIR_IN, FILE_OUT), "w")
ARR_ALL_IN_LINE = FILE_BIB_IN.readlines()

for line in ARR_ALL_IN_LINE:
    #print(len(re.split("{",line)))
    a = re.match("@",line)
    b = re.search("}",line)
    c = re.search(",",line)

    if a == None:
        if line[0] != "}":
            line = line.replace(" = ",' = "')
            line = line.replace("{","")
            line = line.replace("}","")
            if len(line) > 3:
                if line[-2] == ",":
                   # print("FIND , in the end")
                    line = line.replace(',\n', '",\n')
                else:
                    line = line.replace("\n",'"\n')
            line = line.replace("\\textbackslash", "\\{textbackslash}")

    else:
        print(a.group(0))

    FILE_BIB_OUT.write("{0:s}".format( line ))
