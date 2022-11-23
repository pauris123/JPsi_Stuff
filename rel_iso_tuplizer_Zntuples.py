import numpy as np
import pandas as pn
import uproot as up
import awkward as ak
import time
import math

bigBang=time.time()


def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s) #lxplus nepatik d
def l_wall_time(seconds):
    m, s = divmod(seconds, 60)
    return "{0:.0f}min:{1:.0f}s".format(m,s)

def Ele_eff_area(sc_eta):
    if( abs(sc_eta) >= 0.0 and abs(sc_eta) < 1.0 ):
        return 0.1752
        
    if( abs(sc_eta) >= 1.0 and abs(sc_eta) < 1.479 ):
        return 0.1862
        
    if( abs(sc_eta) >= 1.479 and abs(sc_eta) < 2.0 ):
        return 0.1411
        
    if( abs(sc_eta) >= 2.0 and abs(sc_eta) < 2.2 ):
        return 0.1534
        
    if( abs(sc_eta) >= 2.2 and abs(sc_eta) < 2.3 ):
        return 0.1903
        
    if( abs(sc_eta) >= 2.3 and abs(sc_eta) < 2.4 ):
        return 0.2243
        
    if( abs(sc_eta) >= 2.4 ):
        return 0.2687




print("Libraries has been read in!")

# rereco sets
#rootfile1="/eos/cms/store/group/phys_egamma/swmukher/rereco2018/ECAL_NOISE/DYJetsToLLmadgraphMLM.root" # Data Bpark PF (parameters) full
#rootfile1="/eos/cms/store/group/phys_egamma/swmukher/rereco2018/ECAL_NOISE/RunD.root"
#rootfile1="/eos/cms/store/group/phys_egamma/swmukher/rereco2018/ECAL_NOISE/DYToEEpowheg.root" # alt MC

# UL sets
#rootfile1="/eos/cms/store/group/phys_egamma/asroy/Tag-and-Probe_Tree/UL2018_MINIAOD_Nm1/DYJetsToLL_madgraphMLM.root" # MC
rootfile1="/eos/cms/store/group/phys_egamma/asroy/Tag-and-Probe_Tree/UL2018_MINIAOD_Nm1/EGamma_RunD.root" # DATA, run A,B,C,D
#rootfile1="/eos/cms/store/group/phys_egamma/asroy/Tag-and-Probe_Tree/UL2018_MINIAOD_Nm1/DYJetsToLL_amcatnloFXFX.root" # alt MC

tree1 = up.open(rootfile1)["tnpEleIDs/fitter_tree;1"]


print("Tree has been read in! {0}".format(l_wall_time(time.time()-bigBang)))

#DATA_0 = tree1.arrays(["el_chIso","el_eta","el_abseta","el_sc_abseta","el_pt","el_neuIso","el_phoIso","el_sc_eta","event_rho","tag_Ele_pt","tag_Ele_phi","tag_Ele_trigMVA","tag_sc_eta","el_3charge","event_met_pfmet","event_met_pfphi","tag_Ele_q","el_q","pair_mass","passingMVA94XwpHZZisoV2", "mcTrue","passingMVA94XwpLnoisoV2"], library="ak")
DATA_0 = tree1.arrays(["el_chIso","el_eta","el_abseta","el_sc_abseta","el_pt","el_neuIso","el_phoIso","el_sc_eta","event_rho","tag_Ele_pt","tag_Ele_phi","tag_Ele_trigMVA","tag_sc_eta","el_3charge","event_met_pfmet","event_met_pfphi","tag_Ele_q","el_q","pair_mass","passingMVA94XwpHZZisoV2","passingMVA94XwpLnoisoV2"], library="ak")


print("Dataset has been read in! {0}".format(l_wall_time(time.time()-bigBang)))




cut_highpt_0 = DATA_0["el_pt"] < 15
cut_higheta_0 = abs(DATA_0["el_eta"]) < 2.5      
cut_tag_pt = DATA_0["tag_Ele_pt"] > 30
cut_charge = DATA_0["el_3charge"] == 1
cut_tag_MVA = DATA_0["tag_Ele_trigMVA"] > 0.92
cut_tag_sceta = DATA_0["tag_sc_eta"] < 2.17     
cut_charge_1 = DATA_0["tag_Ele_q"]*DATA_0["el_q"] < 0
#cut_mc = DATA_0["mcTrue"] == 1

        
data_cut_0 = cut_highpt_0&cut_higheta_0&cut_tag_pt&cut_charge&cut_tag_MVA&cut_tag_sceta&cut_charge_1#&cut_mc

DATA_00 = DATA_0[data_cut_0]

DATA0 = ak.to_pandas(DATA_00)
        
calc_iso = []
probe_pass = []
probe_pass2 = []

print("Number of envets -> {0}".format(len(DATA0)))

for i in range(len(DATA0)):

    calc_iso.append((DATA0["el_chIso"][i] + max(0.0, DATA0["el_neuIso"][i] + DATA0["el_phoIso"][i] - DATA0["event_rho"][i] * Ele_eff_area(DATA0["el_sc_eta"][i])) )/DATA0["el_pt"][i])
    
    if DATA0["passingMVA94XwpLnoisoV2"][i] == 1:
        probe_pass.append(1)
    else:
        probe_pass.append(0)
    
    if DATA0["passingMVA94XwpHZZisoV2"][i] == 1:
        probe_pass2.append(1)
    else:
        probe_pass2.append(0)
        
DATA1 = DATA0.copy()
    
DATA0["el_Rel_ISO"] = calc_iso
DATA0["ProbePass"] = probe_pass

DATA1["el_Rel_ISO"] = calc_iso
DATA1["ProbePass"] = probe_pass2



cut_iso_1 = DATA0["el_Rel_ISO"] > -0.1
cut_iso_2 = DATA0["el_Rel_ISO"] > -0.1

cut_iso_3 = DATA1["el_Rel_ISO"] > -0.1
cut_iso_4 = DATA1["el_Rel_ISO"] > -0.1


#DATA = DATA0[cut_iso_1&cut_iso_2]

DATA2 = DATA1[cut_iso_3&cut_iso_4]

print("Number of envets written -> {0}".format(len(DATA2)))




#ROOT_File_RESULTS = up.recreate("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/Z_MC_TnP_ISO_file_DATA_010_030_wpL_2.root") #len(names["nos"]) # creating root file

#ROOT_File_RESULTS["Electrons"] = DATA # dataframe #cut_data # adding dataframe with data to root file#
           
#ROOT_File_RESULTS.close 




#ROOT_File_RESULTS = up.recreate("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/Z_altMC_TnP_ISO_file_UL_040_100_HZZ.root") #len(names["nos"]) # creating root file
ROOT_File_RESULTS = up.recreate("/afs/cern.ch/user/n/nstrautn/CMSSW_10_2_22/src/EgammaAnalysis/TnPTreeProducer/python/JPsi_Mini/Z_DATA_avg_iso_UL_4.root") #len(names["nos"]) # creating root file

ROOT_File_RESULTS["Electrons"] = DATA2 # dataframe #cut_data # adding dataframe with data to root file
           
ROOT_File_RESULTS.close  
                  
        

        





