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

file_name_MC = []
avg_reliso_MC  = []
std_avg_relsio_MC  = []
med_reliso_MC = []

file_name_DATA = []
avg_reliso_DATA = []
std_avg_relsio_DATA = []
med_reliso_DATA = []



print("Libraries has been read in!")



path = "/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/"
file_list_MC = ["Z_MC_avg_iso_UL.root"]
file_list_DATA = ["Z_DATA_avg_iso_UL.root"]

#path = '/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_10/src/TnP_fits/TnP_LowPtElectrons/'
#file_list_MC = ["TnPpairs_MC_HZZ_UL.root"] 
#file_list_DATA = ["TnPpairs_DATA_HZZ_UL.root"]

for i in range(len(file_list_MC)):
    
    rootfile1=path+file_list_MC[i]
    rootfile2=path+file_list_DATA[i]
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





    for l in range(len(eta_low)): # loop over barrel and encap ## does not plot electron ID's
        for k in range(len(pt_low)):
        
        # For Z 
            cut_lowpt_0 = DataTable1["el_pt"] > pt_low[k]
            cut_highpt_0 = DataTable1["el_pt"] < pt_high[k]
            cut_loweta_0 = abs(DataTable1["el_eta"]) > eta_low[l]
            cut_higheta_0 = abs(DataTable1["el_eta"]) < eta_high[l]
            #cut_relIso_0 = DataTable1["el_Rel_ISO"] < 1.0

            cut_lowpt_1 = DataTable2["el_pt"] > pt_low[k]
            cut_highpt_1 = DataTable2["el_pt"] < pt_high[k]
            cut_loweta_1 = abs(DataTable2["el_eta"]) > eta_low[l]
            cut_higheta_1 = abs(DataTable2["el_eta"]) < eta_high[l]
            #cut_relIso_1 = DataTable2["el_Rel_ISO"] < 1.0
        
            data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0#&cut_relIso_0
            data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1#&cut_relIso_1

            plt.hist(DataTable2["el_Rel_ISO"], bins = 80,range = [0,4], label = 'DATA')
            plt.xlim(0,4)
            plt.show()

            DATA_1 = DataTable1[data_cut_0]
            DATA_2 = DataTable2[data_cut_1]

            print("File {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(file_list_MC[i],pt_low[k],pt_high[k],section[l],np.mean(DATA_1["el_Rel_ISO"]),np.std(DATA_1["el_Rel_ISO"])))
            print("File {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(file_list_DATA[i],pt_low[k],pt_high[k],section[l],np.mean(DATA_2["el_Rel_ISO"]),np.std(DATA_2["el_Rel_ISO"])))
            print("\n")
            file_name_MC.append(file_list_MC[i])
            avg_reliso_MC.append(np.mean(DATA_1["el_Rel_ISO"]))
            std_avg_relsio_MC.append(np.std(DATA_1["el_Rel_ISO"]))
            med_reliso_MC.append(np.quantile(DATA_1["el_Rel_ISO"],0.5))

            file_name_DATA.append(file_list_DATA[i])
            avg_reliso_DATA.append(np.mean(DATA_2["el_Rel_ISO"]))
            std_avg_relsio_DATA.append(np.std(DATA_2["el_Rel_ISO"]))
            med_reliso_DATA.append(np.quantile(DATA_2["el_Rel_ISO"],0.5))

            pt_bin.append(pt_low[k])
            eta_bin.append(section[l])

            #plt.hist(DATA_2["el_Rel_ISO"], bins = 80,range = [0,4], label = 'DATA')
            #plt.xlim(0,4)
            #plt.show()


            # For J/Psi
            #cut_lowpt_0 = DataTable1["ProbePt"] > pt_low[k]
            #cut_highpt_0 = DataTable1["ProbePt"] < pt_high[k]
            #cut_loweta_0 = abs(DataTable1["ProbeEta"]) > eta_low[l]
            #cut_higheta_0 = abs(DataTable1["ProbeEta"]) < eta_high[l]
            #cut_relIso_0 = DataTable1["ProbeRelISO"] < 1.0

            #cut_lowpt_1 = DataTable2["ProbePt"] > pt_low[k]
            #cut_highpt_1 = DataTable2["ProbePt"] < pt_high[k]
            #cut_loweta_1 = abs(DataTable2["ProbeEta"]) > eta_low[l]
            #cut_higheta_1 = abs(DataTable2["ProbeEta"]) < eta_high[l]
            #cut_relIso_1 = DataTable2["ProbeRelISO"] < 1.0

        
            #data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0&cut_relIso_0
            #data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1&cut_relIso_1

            #plt.hist(DataTable2["ProbeRelISO"], bins = 80,range = [0,4], label = 'DATA')
            #plt.xlim(0,4)
            #plt.show()

            #DATA_1 = DataTable1[data_cut_0]
            #DATA_2 = DataTable2[data_cut_1]

            #print("File {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(file_list_MC[i],pt_low[k],pt_high[k],section[l],np.mean(DATA_1["ProbeRelISO"]),np.std(DATA_1["ProbeRelISO"])))
            #print("File {0}, bin {1}-{2} pT, {3} , Avg rel_Iso {4} , Avg rel_Iso std {5}".format(file_list_DATA[i],pt_low[k],pt_high[k],section[l],np.mean(DATA_2["ProbeRelISO"]),np.std(DATA_2["ProbeRelISO"])))
            #print("\n")
            #file_name_MC.append(file_list_MC[i])
            #avg_reliso_MC.append(np.mean(DATA_1["ProbeRelISO"]))
            #std_avg_relsio_MC.append(np.std(DATA_1["ProbeRelISO"]))
            #med_reliso_MC.append(np.quantile(DATA_1["ProbeRelISO"],0.5))

            #file_name_DATA.append(file_list_DATA[i])
            #avg_reliso_DATA.append(np.mean(DATA_2["ProbeRelISO"]))
            #std_avg_relsio_DATA.append(np.std(DATA_2["ProbeRelISO"]))
            #med_reliso_DATA.append(np.quantile(DATA_2["ProbeRelISO"],0.5))

            #pt_bin.append(pt_low[k])
            #eta_bin.append(section[l])

            #plt.hist(DATA_1["ProbeRelISO"],bins = 80, range = [0,4], label = 'MC')
            #plt.xlim(0,2)
            #plt.show()

            #plt.hist(DATA_2["ProbeRelISO"], bins = 80,range = [0,4], label = 'DATA')
            #plt.xlim(0,4)
            #plt.show()


#fo=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/AVG_med_dataset_UL_rel_iso_values_JPsi_MC_relisoCut.txt","a")
fo=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/AVG_med_dataset_UL_rel_iso_values_Z_MC.txt","a")
for i in range(len(file_name_MC)):
    fo.write(file_name_MC[i]+"\t"+str(eta_bin[i])+"\t"+str(pt_bin[i])+"\t"+str(avg_reliso_MC[i])+"\t"+str(std_avg_relsio_MC[i])+"\t"+str(med_reliso_MC[i])+"\n")
fo.close()

#fi=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/AVG_med_dataset_UL_rel_iso_values_JPsi_DATA_relisoCut.txt","a")
fi=open("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/rel_iso_fits/AVG_med_dataset_UL_rel_iso_values_Z_DATA.txt","a")
for i in range(len(file_name_DATA)):
    fi.write(file_name_DATA[i]+"\t"+str(eta_bin[i])+"\t"+str(pt_bin[i])+"\t"+str(avg_reliso_DATA[i])+"\t"+str(std_avg_relsio_DATA[i])+"\t"+str(med_reliso_DATA[i])+"\n")
fi.close()
    

        



                
#############################################################################################




