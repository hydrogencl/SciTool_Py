#!/usr/bin/python3

import re,os,math
import argparse

# ------------------------------------------------------ #
#    Lib of icon profiling reader/analyser               #
#        Created on 26.09                                #
#                                                        #
#    Purpose:                                            #
#        To read the profiling of ICON smoothly,         #
#        regardless the different output from different  #
#        ICON compoments                                 #
#                                                        #
#    Issue:                                              #
#        The code is not optimized yet, scratch from     #
#        the jupyter notebook                            #
#                                                        #
#                                                        #
# ------------------------------------------------------ #



class icon_profiler:
    #strFileIn
    def __init__(self, fileIn):
        self.strFileIn = fileIn
        self.numStart = None
        self.numEnd   = None

        self.dicAllInOne = {}

        chkNum = 0
        try:
            with open(self.strFileIn) as fileIn:
                arrLines = fileIn.readlines()
                for i,m in enumerate(arrLines):
                    if len(re.findall("Timer report", m.strip())) == 1:
                        self.numStart = i
                    if len(re.findall("-{100}", m.strip())) == 1:
                        self.numEnd = i

        except:
            print("Error: Can't find the file")
        if self.numStart == None: chkNum += 2 ** 1  
        if self.numEnd   == None: chkNum += 2 ** 2 
        if chkNum == 2 ** 1 or chkNum == 2 ** 1 + 2 ** 2:
            print("Can't find the end line (----) of icon profile")
        elif chkNum == 2 ** 2 or chkNum == 2 ** 1 + 2 ** 2:
            print("Can't find the header of icon profile")
        else:
            print("find the time report correctly")
            self.arrLineIn = arrLines[self.numStart : self.numEnd + 1]

    def check_profiling_level(self, strLineIn, arrLevelHeader=["\ {1}(.*)\ {10}", "\ {2,}(L)\ ", "\ {5,}(L)\ ", "\ {8,}(L)\ ", "\ {11,}(L)\ ", "\ {14,}(L)\ " ]):
        numLevel = None
        for i,strRegex in enumerate(arrLevelHeader):
            if len(re.findall(strRegex, strLineIn)) == 1: numLevel = i
        return numLevel
    
        def get_time(strIn):
            numHour   = re.findall('(?<=)([0-9]+)(?=h|H)', strIn)
            numMinute = re.findall('(?<=)([0-9]+)(?=m|M)', strIn) 
            numSecond = re.findall('(?<=)([0-9]+)(?=s|S)', strIn)
            numHour   = int(numHour[0]) if len(numHour)   == 1 else  0
            numMinute = int(numMinute[0]) if len(numMinute) == 1 else 0
            numSecond = int(numSecond[0]) if len(numSecond) == 1 else 0
            return {"TotalSecond": numHour*3600 + numMinute*60 + numSecond, "Hour": numHour, "Minute": numMinute, "Second": numSecond }
    
    def get_icon_profiling(self, strIn, ifDebug=False, ifShow=False):
        #checking the string
        if len(re.findall("^[0-9]+$",strIn)) == 1:
            if ifShow: print("pure number: {}".format(re.findall("^[0-9]+$",strIn)))
            return {"Number": int(re.findall("^[0-9]*$",strIn)[0]), "Type": "Int"}
    
        elif len(re.findall("^[0-9]+\.[0-9]+$",strIn)) == 1:
            if ifShow: print("floating number: {}".format(re.findall("^[0-9]+$",strIn)))
            return {"Number": float(re.findall("^[0-9]+\.[0-9]+$",strIn)[0]), "Type": "Float"}
    
        elif len(re.findall("\[([0-9]+)\]",strIn)) == 1:
            if ifShow: print("ranks equals: {}".format(re.findall("\[([0-9]+)\]",strIn)[0]))
            return {"Rank": int(re.findall("\[([0-9]+)\]",strIn)[0]), "Type": "Rank"}
    
        elif len(re.findall("([0-9]*?\.?[0-9]+)s$",strIn)) >=1 or len(re.findall("([0-9]+)m",strIn)) >=1 or len(re.findall("^([0-9]+)h",strIn)) >=1:
            #print("find time equals: {}".format(re.findall("([0-9]*)s", strIn)[0]))
            #print("strIn: {}".format(strIn))
            hour=re.findall("^([0-9]+)h",strIn)  #if len(re.findall("^([0-9]*)h",strIn) == 0 else  0
            minute=re.findall("([0-9]+)m",strIn) #if len(re.findall("([0-9]*)m",strIn)  == 0 else  0
            second=re.findall("([0-9]*?\.?[0-9]+)s$",strIn) #if len(re.findall("([0-9]*)s",strIn)  == 0 else  0
    
            hour = int(hour[0]) if len(hour) ==1 else 0
            minute = int(minute[0]) if len(minute) ==1 else 0
            second = float(second[0]) if len(second) ==1 else 0
    
            #ifDebug: print("strIn: {}".format())
            if ifShow: print("in: {}, times: {} {} {}".format(strIn, hour, minute, second))
            return {"Hour": hour,"Minute": minute, "Second": second, "TotalSecond": hour*3600 + minute*60 + second, "Type": "Time"}
        else:
            if ifShow: print("Something else")
            return {"Text": strIn, "Type": "Text"}
    
    def get_dict_profiling(arrLines):
        dicOut    = {}
        chkLevel  = 0
        arrLevelStr = ["","","","","","",""]
    
        for i,m in enumerate(arrLines):
            chkLevel              = check_profiling_level(m, arrLevelHeader)
            dataOut               = get_profiling_info(m)
            arrLevelStr[chkLevel] = dataOut["strVar"]
        return

    def get_profiling_info_AiO(self, strLineIn, arrLevelHeader=["\ {1}(.*)\ {10}", "\ {2,}(L)\ ", "\ {5,}(L)\ ", "\ {8,}(L)\ ", "\ {11,}(L)\ ", "\ {14,}(L)\ " ]):
        ifNum = False
        numLevel = None
        #print("str as input:",strLineIn)
        for i,strRegex in enumerate(arrLevelHeader):
            if len(re.findall(strRegex, strLineIn)) == 1: numLevel = i
            
        #print("level: {}".format(numLevel))    
        arrLineIn = re.split("\s{2,}", strLineIn.strip() )    
        if numLevel == None:
            return {"Output": arrLineIn, "Type": "None", "Variable": None, "Level": None}
        elif len(arrLineIn) > 1 and arrLineIn[0]  == 'name':
            ifNum = False
            #print("item", len(arrLineIn))
            return {"Output": arrLineIn, "Type": "Header" }
        
        elif re.findall("\-{10,}", arrLineIn[0]):
            return {"Output": "-"*20, "Type": "Line" }        
        else:
            ifNum = True
            if len(re.findall("(?<=L\ )(.*)", arrLineIn[0])) >0:
                #print("other level")
                strVar = re.findall("(?<=L\ )(.*)", arrLineIn[0])
            elif len(re.findall("^(.*)$", arrLineIn[0])) >0:
                print("base level: {}".format(arrLineIn[0]))
                strVar = re.findall("^(.*)$", arrLineIn[0])
            else:
                print("can't find for {}".format(arrLineIn[0]))
                strVar = None
            return {"Output": arrLineIn, "Type": "Data", "Variable": strVar[0], "Level": numLevel}

    def get_correct_level(self, arrIn, null=""):
        arrOut = []
        for k in arrIn:
            if not k == null:
                arrOut.append(k)
        return arrOut
    
    def get_profiling_information_to_dict(self):
    
        arrLevelStr = ["","","","","",""]
        arrLevelChk = 0
        for i,m in enumerate(self.arrLineIn):
            dicInfo = self.get_profiling_info_AiO(m)
    
            if dicInfo["Type"] == "Header":
                arrVars = dicInfo["Output"]
            elif dicInfo["Type"] == "Data" and dicInfo["Level"] != None:
                #print(dicInfo)
                # Get the level right
                if dicInfo["Level"] < arrLevelChk:
                    for i in range(dicInfo["Level"], len(arrLevelStr)):
                        arrLevelStr[i] = ""
    
                self.dicAllInOne[dicInfo["Variable"]] = { }
                arrLevelStr[dicInfo["Level"]] = dicInfo["Variable"] 
                self.dicAllInOne[dicInfo["Variable"]]["Level"] = dicInfo["Level"]
                #print(dicInfo["Output"])
                for i,v in enumerate(arrVars):
                    dataOut = self.get_icon_profiling(dicInfo["Output"][i])
                    self.dicAllInOne[dicInfo["Variable"]][v] = dataOut
                
                self.dicAllInOne[dicInfo["Variable"]]["LevelStructure"] = self.get_correct_level(arrLevelStr)
                arrLevelChk = dicInfo["Level"]
            else:
                pass
        self.dicAllInOne

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reading input date')
    parser.add_argument("-f","--input-file",  
                         default=False, help='the input for icon.log')
    args = parser.parse_args()

    icon_log_in          = icon_profiler(args.input_file)
    dict_of_icon_profile = icon_log_in.get_profiling_information_to_dict()

    print(icon_log_in.dicAllInOne["total"]["total avg (s)"]["Number"])

