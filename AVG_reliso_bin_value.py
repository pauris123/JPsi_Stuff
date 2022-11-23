import numpy as np
import pandas as pn
import matplotlib.pyplot as plt
import uproot as up
import scipy.optimize as opt
from matplotlib.ticker import NullFormatter
import awkward as ak
import time
bigBang=time.time()

# CHECK IMINUIT LIBRARY AND MINIMIZER FOR FITTING THE PLOTS
def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s) #lxplus nepatik d
def l_wall_time(seconds):
    m, s = divmod(seconds, 60)
    return "{0:.0f}min:{1:.0f}s".format(m,s)



pt_bin  = []
eta_bin  = []

rel_iso_bin_MC = []
avg_reliso_MC  = []
std_avg_relsio_MC  = []

rel_iso_bin_DATA = [] = []
avg_reliso_DATA = []
std_avg_relsio_DATA = []



print("Libraries has been read in!")



path = "/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/"
file_list_MC = "Z_MC_TnP_ISO_file_UL_HZZ.root"
file_list_DATA = "Z_DATA_TnP_ISO_file_UL_HZZ.root"

#path = '/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_10/src/TnP_fits/TnP_LowPtElectrons/'
#file_list_MC = "TnPpairs_MC_HZZ_UL.root" # J/Psi
#file_list_DATA = "TnPpairs_DATA_HZZ_UL.root" # J/Psi


    
rootfile1=path+file_list_MC
rootfile2=path+file_list_DATA
tree1 = up.open(rootfile1)["Electrons;1"] # MC
tree2 = up.open(rootfile2)["Electrons;1"] # DATA

#tree1 = up.open(rootfile1)["Events"] # MC
#tree2 = up.open(rootfile2)["Events"] # DATA

print("Tree has been read in! {0}".format(l_wall_time(time.time()-bigBang)))

keys = ["el_pt","el_eta","el_Rel_ISO","pair_mass"] # for Z
#keys = ["ProbePt","ProbeEta","ProbeRelISO","Diele_mass"] # for JPsi

    
DataTable1 = tree1.arrays(keys, library="ak") # Z MC
DataTable2 = tree2.arrays(keys, library="ak") # Z DATA



print("Dataset has been read in! {0}".format(l_wall_time(time.time()-bigBang)))


    #pt_low = [5,7,10,15,20]
    #pt_high = [7,10,15,20,30]
pt_low = [5,7,10]
pt_high = [7,10,15]

eta_low = [0.0,0.8,1.56]
eta_high = [0.8,1.44,2.5]
section = ['Inner','Outer','Endcap']
iso_low = [-0.01,0.1,0.2,0.4]
iso_high = [0.1,0.2,0.4,1.0]





for l in range(len(eta_low)): # loop over barrel and encap ## does not plot electron ID's
    for k in range(len(pt_low)):
        for i in range(len(iso_low)):
        
        # For Z 
            cut_lowpt_0 = DataTable1["el_pt"] > pt_low[k]
            cut_highpt_0 = DataTable1["el_pt"] < pt_high[k]
            cut_loweta_0 = abs(DataTable1["el_eta"]) > eta_low[l]
            cut_higheta_0 = abs(DataTable1["el_eta"]) < eta_high[l]
            cut_lowiso_0 = DataTable1["el_Rel_ISO"] > iso_low[i]
            cut_highiso_0 = DataTable1["el_Rel_ISO"] < iso_high[i]

            cut_lowpt_1 = DataTable2["el_pt"] > pt_low[k]
            cut_highpt_1 = DataTable2["el_pt"] < pt_high[k]
            cut_loweta_1 = abs(DataTable2["el_eta"]) > eta_low[l]
            cut_higheta_1 = abs(DataTable2["el_eta"]) < eta_high[l]
            cut_lowiso_1 = DataTable2["el_Rel_ISO"] > iso_low[i]
            cut_highiso_1 = DataTable2["el_Rel_ISO"] < iso_high[i]
        
            data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0&cut_lowiso_0&cut_highiso_0
            data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1&cut_lowiso_1&cut_highiso_1

            DATA_1 = DataTable1[data_cut_0]
            DATA_2 = DataTable2[data_cut_1]

            print("Iso_bin {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(i,pt_low[k],pt_high[k],section[l],np.mean(DATA_1["el_Rel_ISO"]),np.std(DATA_1["el_Rel_ISO"])))
            print("Iso_bin {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(i,pt_low[k],pt_high[k],section[l],np.mean(DATA_2["el_Rel_ISO"]),np.std(DATA_2["el_Rel_ISO"])))
            print("\n")
            rel_iso_bin_MC.append(i)
            avg_reliso_MC.append(np.mean(DATA_1["el_Rel_ISO"]))
            std_avg_relsio_MC.append(np.std(DATA_1["el_Rel_ISO"]))

            rel_iso_bin_DATA.append(i)
            avg_reliso_DATA.append(np.mean(DATA_2["el_Rel_ISO"]))
            std_avg_relsio_DATA.append(np.std(DATA_2["el_Rel_ISO"]))

            pt_bin.append(pt_low[k])
            eta_bin.append(section[l])

        # J/Psi stuff
            #cut_lowpt_0 = DataTable1["ProbePt"] > pt_low[k]
            #cut_highpt_0 = DataTable1["ProbePt"] < pt_high[k]
            #cut_loweta_0 = abs(DataTable1["ProbeEta"]) > eta_low[l]
            #cut_higheta_0 = abs(DataTable1["ProbeEta"]) < eta_high[l]
            #cut_lowiso_0 = DataTable1["ProbeRelISO"] > iso_low[i]
            #cut_highiso_0 = DataTable1["ProbeRelISO"] < iso_high[i]

            #cut_lowpt_1 = DataTable2["ProbePt"] > pt_low[k]
            #cut_highpt_1 = DataTable2["ProbePt"] < pt_high[k]
            #cut_loweta_1 = abs(DataTable2["ProbeEta"]) > eta_low[l]
            #cut_higheta_1 = abs(DataTable2["ProbeEta"]) < eta_high[l]
            #cut_lowiso_1 = DataTable2["ProbeRelISO"] > iso_low[i]
            #cut_highiso_1 = DataTable2["ProbeRelISO"] < iso_high[i]
        
            #data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0&cut_lowiso_0&cut_highiso_0
            #data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1&cut_lowiso_1&cut_highiso_1

            #DATA_1 = DataTable1[data_cut_0]
            #DATA_2 = DataTable2[data_cut_1]

            #print("Iso_bin {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(i,pt_low[k],pt_high[k],section[l],np.mean(DATA_1["ProbeRelISO"]),np.std(DATA_1["ProbeRelISO"])))
            #print("Iso_bin {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(i,pt_low[k],pt_high[k],section[l],np.mean(DATA_2["ProbeRelISO"]),np.std(DATA_2["ProbeRelISO"])))
            #print("\n")
            #rel_iso_bin_MC.append(i)
            #avg_reliso_MC.append(np.mean(DATA_1["ProbeRelISO"]))
            #std_avg_relsio_MC.append(np.std(DATA_1["ProbeRelISO"]))

            #rel_iso_bin_DATA.append(i)
            #avg_reliso_DATA.append(np.mean(DATA_2["ProbeRelISO"]))
            #std_avg_relsio_DATA.append(np.std(DATA_2["ProbeRelISO"]))

            #pt_bin.append(pt_low[k])
            #eta_bin.append(section[l])

            #plt.hist(DATA_1["ProbeRelISO"],bins = 20, label = 'MC')
            #plt.hist(DATA_2["ProbeRelISO"], bins = 20, label = 'DATA')
            #plt.xlim(0,1)
            #plt.show()


#fo=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/rel_iso_bin_values_JPsi_MC_UL.txt","a")
fo=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/rel_iso_bin_values_Z_MC_UL.txt","a")
for i in range(len(avg_reliso_MC)):
    fo.write(str(rel_iso_bin_MC[i])+"\t"+str(eta_bin[i])+"\t"+str(pt_bin[i])+"\t"+str(avg_reliso_MC[i])+"\t"+str(std_avg_relsio_MC[i])+"\n")
fo.close()

#fi=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/rel_iso_bin_values_JPsi_DATA_UL.txt","a")
fi=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/rel_iso_bin_values_Z_DATA_UL.txt","a")
for i in range(len(avg_reliso_DATA)):
    fi.write(str(rel_iso_bin_DATA[i])+"\t"+str(eta_bin[i])+"\t"+str(pt_bin[i])+"\t"+str(avg_reliso_DATA[i])+"\t"+str(std_avg_relsio_DATA[i])+"\n")
fi.close()
    

        



                
#############################################################################################



