import sys, subprocess,  math, os, re, shutil, time
from subprocess import Popen
from subprocess import run

class SlurmController:
    def __init__(self, strProject="", strPartition="", numCoresPerNode=0, strMember=".mem", strRunFolder="run", strJobname="", strRootdir=".", ifVerbose=False, strMachine=None, ifServerLog=False, strServerLog="EFL_ServerLog"):
        self.str_project   = strProject
        self.str_partition = strPartition
        self.corespernode  = numCoresPerNode
        self.str_pre       = strMember
        self.folder_run    = strRunFolder
        self.if_verbose    = ifVerbose
        if strJobname != "":
            self.str_jobname  = strJobname 
        else:
            self.str_jobname  = "EnsemPy"
        self.str_outname   = "{0:s}-out".format(self.str_jobname)
        self.str_errname   = "{0:s}-err".format(self.str_jobname)
        self.str_rootdir   = strRootdir
        if strMachine == "JUWELS":
            self.corespernode  = 48
        if strMachine == "JURECA":
            self.corespernode  = 64
        self.if_serverlog      = ifServerLog
        if ifServerLog:
            self.strServerLog = "{0:s}/{1:s}".format(strRootdir, strServerLog)
            fileLog = open(self.strServerLog, "w")
            fileLog.close()
        self.num_waiting   = 5
        self.num_ping_file = 5
        self.StartTime     = time.time_ns() * 10 ** -9

    def ChangeServerLog(self, fileOutput):
        self.fileServerLog = open(fileOutput, "w")

    def InitEnsemble(self, members=1, UsingNodes=1, UsingCores=None ):
        if UsingCores == None:
            self.corespermember = int(UsingNodes*self.corespernode)
        else:
            self.corespermember = UsingCores
        self.nodespermember = UsingNodes
        self.nummembers     = members

        self.arr_hostpermember = [ "" for m in range(members)]
        try:
            len(self.arr_hostnames)
        except:
            print("No Nodelist, obtaining automatically")
            self.CheckNodelist()

        for m in range(members):
            for n in range(UsingNodes):
                if n == 0:
                    self.arr_hostpermember[m] = "{}".format(self.arr_hostnames[ m*UsingNodes + n ])
                elif n == UsingNodes-1:           
                    str_pre = self.arr_hostpermember[m]
                    self.arr_hostpermember[m] = "{},{}".format(str_pre, self.arr_hostnames[ m*UsingNodes + n ])
                else:
                    str_pre = self.arr_hostpermember[m]
                    self.arr_hostpermember[m] = "{},{}".format(str_pre, self.arr_hostnames[ m*UsingNodes + n ])

    def CheckNodelist(self):
        str_hostnames     = subprocess.run(["scontrol", "show", "hostnames"], capture_output=True).stdout
        arr_hostnames_tmp = str_hostnames.split(b"\n")
        self.arr_hostnames       = []
        for ind, item in enumerate(arr_hostnames_tmp):
            if item != b'':
                self.arr_hostnames.append(item.decode('UTF-8'))

    def CreateAMember(self, numMember, arrException=[] ):

        self.arr_runfolder_files = os.listdir("{0:s}/{1:s}".format(self.str_rootdir, self.folder_run))
        for ind, item in enumerate(self.arr_runfolder_files):
            for item2 in arrException:
                if item2 == item:
                    self.arr_runfolder_files.pop(ind)

        str_runfolder_out = "{0:s}/{1:s}{2:s}{3:04d}".format(self.str_rootdir, self.folder_run, self.str_pre, numMember)
        if not os.path.exists(str_runfolder_out): 
            os.mkdir("{0:s}".format(str_runfolder_out))
            os.chdir("{0:s}".format(str_runfolder_out))
            for item in self.arr_runfolder_files:
                subprocess.run(["ln","-s", "{0:s}/{1:s}/{2:s}".format(self.str_rootdir, self.folder_run, item), "."])

    def CreateMembers(self, arrException=[], if_force=True, if_hardlink=False):

        self.arr_runfolder_files = os.listdir("{0:s}/{1:s}".format(self.str_rootdir, self.folder_run))
        for ind, item in enumerate(self.arr_runfolder_files):
            for item2 in arrException:
                if item2 == item:
                    self.arr_runfolder_files.pop(ind)

        for ind in range(self.nummembers):
            str_runfolder_out = "{0:s}/{1:s}{2:s}{3:04d}".format(self.str_rootdir, self.folder_run, self.str_pre, ind)
            if if_force:
                if os.path.exists(str_runfolder_out):
                    shutil.rmtree(str_runfolder_out)

            if not os.path.exists(str_runfolder_out): 
                os.mkdir("{0:s}".format(str_runfolder_out))
                os.chdir("{0:s}".format(str_runfolder_out))
                for item in self.arr_runfolder_files:
                    strFileIn  = "{0:s}/{1:s}/{2:s}".format(self.str_rootdir, self.folder_run, item)
                    strFileOut = "{0:s}/{1:s}".format(str_runfolder_out, item)
                    if if_hardlink:
                        os.link(strFileIn, strFileOut)
                    else:
                        os.symlink(strFileIn, strFileOut)
                    #subprocess.run(["ln","-s", "{0:s}/{1:s}/{2:s}".format(self.str_rootdir, self.folder_run, item), "."])
                os.chdir(self.str_rootdir)

    def FileControl(self, strTargetFile, strAction, strSourceFile=""):
        for ind in range(self.nummembers):
            str_runfolder_out = "{3:s}/{0:s}{1:s}{2:04d}".format(self.folder_run, self.str_pre, ind, self.str_rootdir)
            strTargetPath = "{0:s}/{1:s}".format(str_runfolder_out, strTargetFile)
            if strAction == "remove":
                try:
                    os.remove(strTargetPath)
                except:
                    self.ServerLog("{0:s} is not existed. Skip".format(strTargetPath))
         
            elif strAction == "copy":
                strSourcePath = "{0:s}/{1:s}".format(str_runfolder_out, strSourceFile) 
                try:
                    os.copyfile(strTargetPath, strSourcePath)
                except:
                    self.ServerLog("{0:s} is not existed. Skip".format(strSourcePath))
            self.ServerLog("mem {0:04d} is done for the file control".format(ind))


    def RunMembers(self, strExecutor=""):
        if strExecutor == "":
            print("FATAL ERROR: you did not specific the executor's name") 
        for ind in range(self.nummembers):
            str_runfolder_out = "{3:s}/{0:s}{1:s}{2:04d}".format(self.folder_run, self.str_pre, ind, self.str_rootdir)
            os.chdir("{0:s}".format(str_runfolder_out))
            self.ServerLog("Mem: {3:04d}, Nodes: {0:4d}, Cores: {1:4d}, Nodelist: {2:s}"\
                           .format(self.nodespermember, self.corespermember, self.arr_hostpermember[ind], ind))
            Popen(["srun", \
                            "-N", "{0:d}".format(self.nodespermember),\
                            "-n", "{0:d}".format(self.corespermember),\
                            "--verbose",\
                            "--nodelist={0:s}".format(self.arr_hostpermember[ind]),\
                            "--output={0:s}".format(self.str_outname),\
                            "--error={0:s}".format(self.str_outname),\
                            "--job-name={0:s}".format(self.str_jobname),\
                            "{0:s}".format(strExecutor),\
                            "&"], 
                            stdin=None, stdout=None, stderr=None, close_fds=True)
            os.chdir("{0:s}".format(self.str_rootdir))

    def CheckMembersWRF(self, numWaitingTime=5 , ifExit=True, numWRFinitTime=10):
        time.sleep(numWRFinitTime)
        num_finished   = 0
        arrCheckEnding = [ 0 for n in range(self.nummembers)]
        arrCheckError  = [ 0 for n in range(self.nummembers)]
        ifRunning      = True
        while ifRunning:
            time.sleep(numWaitingTime)
            for ind in range(self.nummembers):
                str_runfolder_out = "{3:s}/{0:s}{1:s}{2:04d}".format(self.folder_run, self.str_pre, ind, self.str_rootdir)
                try:
                    with open("{0:s}/rsl.out.0000".format(str_runfolder_out), "r") as fileOut:
                        text_out = fileOut.readlines()[-1]
                    with open("{0:s}/rsl.error.0000".format(str_runfolder_out), "r") as fileOut:
                        text_err = fileOut.readlines()[-1]
                except:
                    self.ServerLog("Can not find the log files. Will retry after {0:d} seconds".format(numWRFinitTime))
                    self.ServerLog("Folder to log: {0:s}".format(str_runfolder_out))
                    time.sleep(numWRFinitTime)
                    ifContinueSearch = True
                    num_Retry        = 0
                    if ifContinueSearch and num_Retry < self.num_ping_file:
                        print("check if exist, retry {0:d}/{1:d}".format(num_Retry,self.num_ping_file))
                        chkOut = os.path.exists("{0:s}/rsl.out.0000".format(str_runfolder_out))
                        chkErr = os.path.exists("{0:s}/rsl.error.0000".format(str_runfolder_out))
                        if chkOut and chkErr:
                            ifContinueSearch = False
                        time.sleep(numWaitingTime)
                        num_Retry += 1
                    with open("{0:s}/rsl.out.0000".format(str_runfolder_out), "r") as fileOut:
                        text_out = fileOut.readlines()[-1]
                    with open("{0:s}/rsl.error.0000".format(str_runfolder_out), "r") as fileOut:
                        text_err = fileOut.readlines()[-1]
                    self.ServerLog("Find the log files of WRF.".format(numWRFinitTime))
    
                chk_out_tmp = self.CheckWRFexit(text_out)["OUT"]
                chk_err_tmp = self.CheckWRFexit(text_err)["ERROR"]
                arrCheckEnding[ind] = chk_out_tmp
                arrCheckError [ind] = chk_err_tmp
            if sum(arrCheckEnding) == 0:
                self.ServerLog("Still running")
            else:
                if num_finished == sum(arrCheckEnding):
                    self.ServerLog("Already Finished ({0:d}/{1:d})".format(num_finished, self.nummembers))
                else:
                    self.ServerLog("Adding New Finished: ")
                    for mem in range(self.nummembers):
                        if arrCheckEnding[mem] == 1:
                            self.ServerLog("Member: {0:04d}, Finished".format(mem))
                        else:
                            self.ServerLog("Member: {0:04d}, Still Running".format(mem))
                    num_finished = sum(arrCheckEnding)
            if num_finished == self.nummembers: 
                self.ServerLog("All members finished, stop listening")
                self.ServerLog("El Psy Congroo                      ")
                ifRunning = False
                if ifExit : sys.exit()

    def ServerLog(self, strDescription="", strRuntimeFormat="passtime"):
        """ Two time format: passtime or exacttime
        """

        time_now = time.gmtime()
        self.CheckTime = time.time_ns() * 10 ** -9
        self.PassTime  = self.CheckTime - self.StartTime
        str_time_now = "{0:04d}-{1:02d}-{2:02d}_{3:02d}:{4:02d}:{5:02d}".format(    \
                        time_now. tm_year  , time_now. tm_mon  , time_now. tm_mday ,\
                        time_now. tm_hour  , time_now. tm_min  , time_now. tm_sec )
        fileServerLog = open(self.strServerLog, "a")
        if strRuntimeFormat == "exacttime":
            fileServerLog.write("{0:10s} : {1:s}\n".format(str_time_now, strDescription))
        elif strRuntimeFormat == "passtime":
            fileServerLog.write("{0:10.1f} : {1:s}\n".format(self.PassTime, strDescription))
        else:
            fileServerLog.write("{0:10.1f} : {1:s}\n".format(self.PassTime, strDescription))
        fileServerLog.close()


    def CheckWRFexit(self, txt_in, strSuccess="SUCCESS COMPLETE WRF", strTerminate="terminated abnormally"):
        """Get the signal from rsl.out, check only if WRF is sucessfully done or fail
           Success log: wrf: SUCCESS COMPLETE WRF
        """
        find_out = re.findall(strSuccess,   txt_in)
        find_err = re.findall(strTerminate, txt_in)
        return {"OUT":  len(find_out) , "ERROR":  len(find_err) }




