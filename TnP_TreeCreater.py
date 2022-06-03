import numpy as np
import pandas as pn
import matplotlib.pyplot as plt
import uproot3 as up3
import uproot as up
from matplotlib.ticker import NullFormatter
import awkward as ak
import time
bigBang=time.time()


def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s) #lxplus nepatik d
def l_wall_time(seconds):
    m, s = divmod(seconds, 60)
    return "{0:.0f}min:{1:.0f}s".format(m,s)

print("Libraries has been read in!")

#rootfile1="/eos/user/n/nstrautn/Bparking_DATA/Bpark_DATA_MassCut_2018.root" # Data Bpark PF (parameters) full
rootfile1="/eos/user/n/nstrautn/Bparking_DATA/Bpark_MC_MassCut_FULL.root" # MC Bpark
tree1 = up.open(rootfile1)["Electrons;1"]
#tree1 = up.open(rootfile1)["electrons/Events;1"]
#tree2 = up.open(rootfile2)["Electrons;1"]
#tree2 = up.open(rootfile2)["electrons/Events;1"]

print("Tree has been read in! {0}".format(l_wall_time(time.time()-bigBang)))

DataTable1 = tree1.arrays(tree1.keys(), library="ak")
#DataTable2 = tree2.arrays(tree2.keys(), library="ak")

#DataTable3 = tree3.arrays(tree3.keys(), library="pd")

print("Dataset has been read in! {0}".format(l_wall_time(time.time()-bigBang)))

pair_list = []

ID_list = ["cutBasedElectronID_Fall17_94X_V2_loose_0","cutBasedElectronID_Fall17_94X_V2_medium_0","cutBasedElectronID_Fall17_94X_V2_tight_0","cutBasedElectronID_Fall17_94X_V2_veto_0","mvaEleID_Fall17_iso_V2_wp80_0","mvaEleID_Fall17_iso_V2_wp90_0","mvaEleID_Fall17_iso_V2_wpHZZ_unsopported_0","mvaEleID_Fall17_noIso_V2_wp80_0","mvaEleID_Fall17_noIso_V2_wp90_0","mvaEleID_Fall17_noIso_V2_wpLoose_unsopported_0"]


#eta_cut_11 = DataTable1["ele_eta_0"] < 1.444 
#eta_cut_22 = DataTable1["ele_eta_0"] > -1.444 
#eta_cut_33 = DataTable1["ele_eta_0"] > 1.566 
#eta_cut_44 = DataTable1["ele_eta_0"] < -1.566 

#eta_1 = eta_cut_11 & eta_cut_22
#eta_2 = eta_cut_11 and eta_cut_22

# Take 1st as tag
pt_cut_0 = DataTable1["ele_pt_0"] > 10
ID_cut_0 = DataTable1["mvaEleID_Fall17_noIso_V2_wp80_0"] == 1
ID_cut_0_0 = DataTable1["mvaEleID_Fall17_noIso_V2_wp80_1"] == 0
eta_cut_0_0 = abs(DataTable1["ele_eta_0"]) > 1.566
eta_cut_1_0 = abs(DataTable1["ele_eta_0"]) < 1.444

eta_cut_0 = eta_cut_0_0 | eta_cut_1_0 #
#eta_cut_00 = eta_cut_0_0 or eta_cut_1_0 #

DATA_0 = DataTable1[pt_cut_0&ID_cut_0&ID_cut_0_0&eta_cut_0]
tag_list_0 = [0 for i in range(len(DATA_0))]
DATA_0["tag_list"] = tag_list_0

# Take 2nd as tag
pt_cut_1 = DataTable1["ele_pt_1"] > 10
ID_cut_1 = DataTable1["mvaEleID_Fall17_noIso_V2_wp80_1"] == 1
ID_cut_1_1 = DataTable1["mvaEleID_Fall17_noIso_V2_wp80_0"] == 0
eta_cut_0_1 = abs(DataTable1["ele_eta_1"]) > 1.566
eta_cut_1_1 = abs(DataTable1["ele_eta_1"]) < 1.444

eta_cut_1 = eta_cut_0_1 | eta_cut_1_1

DATA_1 = DataTable1[pt_cut_1&ID_cut_1&ID_cut_1_1&eta_cut_1]
tag_list_1 = [1 for i in range(len(DATA_1))]
DATA_1["tag_list"] = tag_list_1

# Take random as tag

DATA_rand = DataTable1[pt_cut_0&ID_cut_0&eta_cut_0&pt_cut_1&ID_cut_1&eta_cut_1]
tag_list_rand = [np.random.randint(2) for i in range(len(DATA_rand))]
DATA_rand["tag_list"] = tag_list_rand

DATA = ak.concatenate([DATA_0,DATA_1,DATA_rand])

print("1st electrons as tags - {0}".format(len(DATA_0)))
print("2nd electrons as tags - {0}".format(len(DATA_1)))
print("Both electrons as tags (randomized) - {0}".format(len(DATA_rand)))

#for i in range(len(DataTable1)):
#
#    if (DataTable1["ele_pt_0"][i] > 10 & DataTable1["mvaEleID_Fall17_noIso_V2_wp80_0"][i] == 1 & (abs(DataTable1["ele_eta_0"][i]) > 1.566 or abs(DataTable1["ele_eta_0"][i]) < 1.444) & DataTable1[##"mvaEleID_Fall17_noIso_V2_wp80_1"][i] == 0) :
#        # Take 1st as tag
#        tag_list.append(0)
#        pair_list.append(True)
        
#    elif (DataTable1["ele_pt_1"][i] > 10 & DataTable1["mvaEleID_Fall17_noIso_V2_wp80_1"][i] == 1 & (abs(DataTable1["ele_eta_1"][i]) > 1.566 or abs(DataTable1["ele_eta_1"][i]) < 1.444) & DataTable1#["mvaEleID_Fall17_noIso_V2_wp80_0"][i] == 0) :
#        # Take 2nd as tag
#        tag_list.append(1)
#        pair_list.append(True)
        
#    elif (DataTable1["ele_pt_0"][i] > 10 & DataTable1["mvaEleID_Fall17_noIso_V2_wp80_0"][i] == 1 & (abs(DataTable1["ele_eta_0"][i]) > 1.566 or abs(DataTable1["ele_eta_0"][i]) < 1.444)) & #(DataTable1["ele_pt_1"][i] > 10 & DataTable1["mvaEleID_Fall17_noIso_V2_wp80_1"][i] == 1 & (abs(DataTable1["ele_eta_1"][i]) > 1.566 or abs(DataTable1["ele_eta_1"][i]) < 1.444)) :
#        # Take random one as tag
#        tag_list.append(np.random.randint(2))
#        pair_list.append(True)
#        
#    else:
#        pair_list.append(False)
    

dataframe = ak.to_pandas(DATA)
dataframe = dataframe.droplevel('subentry')
ROOT_File_RESULTS = up.recreate("/eos/user/n/nstrautn/Bparking_DATA/Bpark_MC_TnP_tree_test.root") #len(names["nos"]) # creating root file

ROOT_File_RESULTS["Electrons"] = dataframe #cut_data # adding dataframe with data to root file
           
ROOT_File_RESULTS.close  

