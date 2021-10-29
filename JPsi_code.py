# Code for J/Psi with uproot 4 and 3
import time
bigBang=time.time()
import uproot as up
import uproot3 as up3
#import uproot_methods as up_me
from ROOT import TLorentzVector
import numpy as np
import pandas as pn
import sys 
import awkward as ak
from uproot3_methods import PtEtaPhiMassLorentzVector

#up.recreate("JPsi/2017_JPsi_MuonEG_RESULTS.root") # !!!!!! Always check this line so, I would not rewrite all the files!!!!!!!!!!

def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s)
def l_wall_time(seconds):
    m, s = divmod(seconds, 60)
    return "{0:.0f}min:{1:.0f}s".format(m,s)

File_Names = pn.read_csv("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/FileNames_2017_v2_JPsi_MuonEG.txt",sep='\t',names=["FileNames"])

print("Got all the libraries and filenames")

coeficient = 1
From = 0+coeficient*int(sys.argv[1])
To = 0+coeficient*int(sys.argv[1])+1

if To>len(File_Names["FileNames"]):
    To=len(File_Names["FileNames"])


#DataFF = pn.DataFrame()



for i in range(From,To):
    
    start=time.time()
    
    print("Working on file number {0}/{1}, {2}".format(i+1,To,wall_time(time.time()-bigBang)))
    
    rootfile="root://cmsxrootd.fnal.gov///"+File_Names["FileNames"][i]
    #tree = up.open(rootfile)["Events"]
    up3_tree = up3.open(rootfile)["Events"]
    
    #tree.arrays(tree.keys(), library="pd")
    
    print("Tree has been read in! {0}".format(l_wall_time(time.time()-start)))
    
    hlt_list = []
    Electron_list = []  
    
    for x in up3_tree.keys():
        if ('HLT_') in x :
            hlt_list.append(x)
        if ('Electron') in x :
            Electron_list.append(x)  
    
    Electron_list.remove('Jet_nElectrons')
    #branches = ["Electron_pt","nElectron","Electron_eta","Electron_charge","Electron_phi","Electron_cutBased","Electron_r9"]
    subsection = Electron_list + hlt_list
    
    
    
    DataTable = up3_tree.pandas.df(subsection, entrystop=-1) # Works for uproot 3
    #tree.arrays(subsection, library="pd")
  
    print("I have the branches! {0}".format(l_wall_time(time.time()-start)))
    
    two_electrons_mask = DataTable["nElectron"] == 2
    DF = DataTable[two_electrons_mask]
    #DF = DF.reset_index() - > Partaisa uz veco uproot 3 formatu, kur subentry un entry ir viena limeni, nevar izveleties tikai 1. vai 2. elektronu (query, protams, strada)
    
    gut_eta=DF["Electron_eta"][:,0].abs() < 1.4
    gut_pt=DF["Electron_pt"][:,0] > 7.0
    gut_ID=DF["Electron_cutBased"][:,0] >= 2
    gut_charge=DF["Electron_charge"][:,0] + DF["Electron_charge"][:,1] == 0
   
    fav=gut_charge&gut_pt&gut_eta&gut_ID
    df=DF[fav[DF.index.get_level_values('entry')].values] # Partaisa no 1 uz 2 per entry true/false stuff, tad var uzlikt masku uz pilnda DF
    
    #df = DataTable.query("nElectron==2")
    #df = df.reset_index()
    #df0 = df.query("subentry==0")
    #df1 = df.query("subentry==1")
    
    print("I have {0} diE events! {1}".format(len(df["nElectron"])/2,l_wall_time(time.time()-start)))
    
    #Electron_list.remove('nElectron') # Jautajums par so
    
    el_InvMass = []
    el_event_i = []
    el_data_filter = []    
    
    entry_N = df.index.get_level_values("entry")
    entry_df = entry_N[::2]
    #subentry_M = df.index.get_level_values("subentry")
    
    
    for y in entry_df:
            
        el1 = PtEtaPhiMassLorentzVector(pt=df["Electron_pt"][y][0], eta=df["Electron_eta"][y][0], phi=df["Electron_phi"][y][0], mass=0.000511) 
        el2 = PtEtaPhiMassLorentzVector(pt=df["Electron_pt"][y][1], eta=df["Electron_eta"][y][1], phi=df["Electron_phi"][y][1], mass=0.000511)

        if ((el1+el2).mass < 20):
            el_InvMass.append((el1+el2).mass)
            el_event_i.append(y)
            el_data_filter.append(True)
            el_data_filter.append(True)
        
        else:
            el_data_filter.append(False)
            el_data_filter.append(False)    
            
    dff = df[el_data_filter]
    df0=dff.xs(0,level=1)
    df1=dff.xs(1,level=1)
    DataF = df0.join(df1,lsuffix="_0",rsuffix="_1")
    DataF["DiElectron_InvMass"] = el_InvMass
    
    #DataFF = DataFF.append(DataF)

    
    print("Got invariant masses, hehe")
    print("I have {0} raw data events! {1}".format(len(el_InvMass),l_wall_time(time.time()-start)))
    print(" Each event took {0:.2f} ms".format(1000*(time.time()-start)/len(el_InvMass)))

print("ROOT file writing has started! {0}".format(wall_time(time.time()-bigBang)))     



# ROOT file for RESULT data
ROOT_File_RESULTS = up.recreate("/afs/cern.ch/user/n/nstrautn/Results/2017_JPsi_MuonEG_RESULTS_{0}_{1}.root".format(From,To))
#ROOT_File_RESULTS = up.update("JPsi/2017_JPsi_MuonEG_RESULTS.root")

ROOT_File_RESULTS["Electrons"] = DataF  # DataF
#ROOT_File_RESULTS["Electrons"] = up.newtree({DataF.columns:dataF.dtypes})       
#ROOT_File_RESULTS["Electrons"].extend({DataF.columns: DataF})


ROOT_File_RESULTS.close 

print("All done! {0}".format(wall_time(time.time()-bigBang)))  


