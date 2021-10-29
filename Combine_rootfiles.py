import time
bigBang=time.time()
import glob,os
import pandas as pn
import uproot as up
import uproot3 as up3
import numpy as np

def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s)

print("I got libraries".format(wall_time(time.time()-bigBang)))
DataCombined = pn.DataFrame()
print("I got empty dataset".format(wall_time(time.time()-bigBang)))
n=0
num=0

for file in glob.glob("/afs/cern.ch/user/n/nstrautn/Results/2017_JPsi_MuonEG_RESULTS_*.root"):
    n+=1
    print("{0}/{1}\t{2}".format(n,len(glob.glob("/afs/cern.ch/user/n/nstrautn/Results/2017_JPsi_MuonEG_RESULTS_*.root")),
                                wall_time(time.time()-bigBang)))
    
    rootfile=file
    up3_tree = up3.open(rootfile)["Electrons"]
    
    hlt = []
    HLT = []
    hlt_list = []
    Electron_list = []
    
    DataTable = up3_tree.pandas.df(up3_tree.keys(), entrystop=-1) 
    
#    for x in up3_tree.keys():
#        if ('HLT_') in x :
#            hlt_list.append(x)
#        if ('Electron') in x :
#            Electron_list.append(x)
    
    
    
#    for z in hlt_list:
#        Check = []
#        for y in range(len(DataTable["Electron_pt_0"])):
#            if DataTable[z][y] == True:
#                Check.append(z)
#        if len(Check) > 1:
#            hlt.append(z)
#            if len(Check) > 0.1*len(DataTable["Electron_pt_0"]):
#                HLT.append(z)
                
#    print("Lenght of hlt list is -> "+str(len(hlt_list)))
#    print("We have used HLT_1 paths (>1) -> "+str(len(hlt)))
#    print("We have used 10% HLT_1 paths (>1) -> "+str(len(HLT)))
    
    
    
    DataCombined = pn.concat([DataCombined,DataTable],axis=0,ignore_index=True)
    
#    DataCombined = DataCombined.append(DataTable)
    
    print(" has {0}".format(len(DataTable["DiElectron_InvMass"])))
    num+=len(DataTable["DiElectron_InvMass"])

print("Started to rewrite NaN to False")
#for i in DataCombined.keys():
#    for y in range(len(DataCombined["Electron_pt_0"])):
#        if np.isnan(DataCombined[i][y]) == True:
#            DataCombined[i][y] = False
            
#DataCombined[np.isnan(DataCombined) == True] = False
nan_value=float("NaN")
DataCombined.replace(nan_value,False,inplace=True)

print("Started writing root file".format(wall_time(time.time()-bigBang)))
ROOT_File_RESULTS = up.recreate("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi/2017_JPsi_MuonEG_RESULTS_ALL.root")
ROOT_File_RESULTS["Electrons"] = DataCombined
ROOT_File_RESULTS.close 
    
print("{0} All done".format(wall_time(time.time()-bigBang)))
print("Had {0} or {1:.2f} M or {2:.2f} G".format(num,num/10**6,num/10**9))


























