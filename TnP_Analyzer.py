import numpy as np
import pandas as pn
import matplotlib.pyplot as plt
import uproot3 as up3
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

#def cBall_exp(x,alpha_l,alpha_h,n_l,n_h,mean,sigma,N,a,b,c,d): # signal fit
#    y=np.piecewise(x, 
#                 [((x-mean)/sigma >= -alpha_l) & ((x-mean)/sigma <= alpha_h), 
#                  (x-mean)/sigma < -alpha_l, 
#                  (x-mean)/sigma > alpha_h], 
#                 [lambda x : N*np.exp(-0.5*((x-mean)/sigma)**2)+(a*x**3+b*x**2+c*x+d),#+c*np.exp(-d*x) , #a*np.exp(-b*x)
#                  lambda x : N*np.exp(-0.5*alpha_l**2)*(alpha_l/n_l*(n_l/alpha_l-alpha_l-(x-mean)/sigma))**(-n_l)+(a*x**3+b*x**2+c*x+d),#+c*np.exp(-d*x), 
#                  lambda x : N*np.exp(-0.5*alpha_h**2)*(alpha_h/n_h*(n_h/alpha_h-alpha_h+(x-mean)/sigma))**(-n_h)+(a*x**3+b*x**2+c*x+d)]#+c*np.exp(-d*x)]
#                )
#    return y

def cBall_exp(x,alpha_l,alpha_h,n_l,n_h,mean,sigma,N): # signal fit
    y=np.piecewise(x, 
                 [((x-mean)/sigma >= -alpha_l) & ((x-mean)/sigma <= alpha_h), 
                  (x-mean)/sigma < -alpha_l, 
                  (x-mean)/sigma > alpha_h], 
                 [lambda x : N*np.exp(-0.5*((x-mean)/sigma)**2),#+c*np.exp(-d*x) , #a*np.exp(-b*x)
                  lambda x : N*np.exp(-0.5*alpha_l**2)*(alpha_l/n_l*(n_l/alpha_l-alpha_l-(x-mean)/sigma))**(-n_l),#+c*np.exp(-d*x), 
                  lambda x : N*np.exp(-0.5*alpha_h**2)*(alpha_h/n_h*(n_h/alpha_h-alpha_h+(x-mean)/sigma))**(-n_h)]#+c*np.exp(-d*x)]
                )
    return y
    
#def linear(x,a,b,c,d):
#    return a*x**3+b*x**2+c*x+d



print("Libraries has been read in!")

rootfile1="/eos/user/n/nstrautn/Bparking_DATA/Bpark_MC_TnP_tree_test.root" # Data Bpark PF (parameters) full
#rootfile1="/eos/user/n/nstrautn/Bparking_DATA/Bpark_MC_TnP_tree_tag7pT.root" # Data Bpark PF (parameters) full
#rootfile2="/eos/user/n/nstrautn/Bparking_DATA/Bpark_MC_MassCut_FULL.root" # MC Bpark
tree1 = up.open(rootfile1)["Electrons;1"]
#tree1 = up.open(rootfile1)["electrons/Events;1"]
#tree2 = up.open(rootfile2)["Electrons;1"]
#tree2 = up.open(rootfile2)["electrons/Events;1"]

print("Tree has been read in! {0}".format(l_wall_time(time.time()-bigBang)))

DataTable1 = tree1.arrays(tree1.keys(), library="ak")
#DataTable2 = tree2.arrays(tree2.keys(), library="ak")


DATA_0 = DataTable1[ DataTable1["tag_list"] == 0 ]
DATA_1 = DataTable1[ DataTable1["tag_list"] == 1 ]

print("Dataset has been read in! {0}".format(l_wall_time(time.time()-bigBang)))


#pt_low = [5,7,10,15,20]
#pt_high = [7,10,15,20,30]
pt_low = [5,6,7,8,9,10,11,12,13,14,15,20]
pt_high = [6,7,8,9,10,11,12,13,14,15,20,30]

eta_low = [0,1.56]
eta_high = [1.44,2.5]
section = ['Barrel','Endcap']

#k = 4 # 0-4, choose value for the bin
#l = 1 # 0 for barrel ; 1 for encap (choose)

# Probe cuts in bins
#cut_lowpt_0 = DATA_0["ele_pt_1"] > pt_low[k]
#cut_highpt_0 = DATA_0["ele_pt_1"] < pt_high[k]
#cut_loweta_0 = abs(DATA_0["ele_eta_1"]) > eta_low[l]
#cut_higheta_0 = abs(DATA_0["ele_eta_1"]) < eta_high[l]

#cut_lowpt_1 = DATA_1["ele_pt_0"] > pt_low[k]
#cut_highpt_1 = DATA_1["ele_pt_0"] < pt_high[k]
#cut_loweta_1 = abs(DATA_1["ele_eta_0"]) > eta_low[l]
#cut_higheta_1 = abs(DATA_1["ele_eta_0"]) < eta_high[l]
        
#data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0
#data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1

#DATA_00 = DATA_0[data_cut_0]
#DATA_11 = DATA_1[data_cut_1]

plot_list = tree1.keys()
plot_list.sort()
plot_list.remove("Diele_mass_0")
plot_list.remove("Diele_mass_1")
plot_list.remove("numele_0")
plot_list.remove("numele_1")
plot_list.remove("tag_list")
plot_list.remove("index")

ID_list = ["cutBasedElectronID_Fall17_94X_V2_loose_0","cutBasedElectronID_Fall17_94X_V2_medium_0","cutBasedElectronID_Fall17_94X_V2_tight_0","cutBasedElectronID_Fall17_94X_V2_veto_0","mvaEleID_Fall17_iso_V2_wp80_0","mvaEleID_Fall17_iso_V2_wp90_0","mvaEleID_Fall17_iso_V2_wpHZZ_unsopported_0","mvaEleID_Fall17_noIso_V2_wp80_0","mvaEleID_Fall17_noIso_V2_wp90_0","mvaEleID_Fall17_noIso_V2_wpLoose_unsopported_0"]

#list1=np.append(DATA_00[plot_list[0]],DATA_11[plot_list[1]])
#list1=DATA_00[plot_list[72]]
#list2=np.append(DATA_00[plot_list[1]],DATA_11[plot_list[0]])
#list2=DATA_11[plot_list[73]]
#plt.hist(list1,bins=100,  alpha=0.9, histtype=u'step',label='tags', color='b',normed = True)
#plt.hist(list2,bins=100,  alpha=0.9, histtype=u'step',label='probes', color='r',normed = True)
#plt.hist(DataTable2[i],bins=2,  alpha=0.9, histtype=u'step',label='MC', color='r') #u'step' normed=True
#plt.legend(loc='best')
#plt.xlabel("{0}".format(i), size = 15)
#plt.grid()
#plt.ylabel("Event count", size = 10)
#plt.title("{0} comparison".format(i), size = 15)
#plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_bin_comparison_tagWp80_10pT/Barrel/5_7_pT/index.php/TnP_comp_{0}_pT.png".format(i))
#plt.show()



# stock plotter, vajag nomainit destination un izveidot atlasi pa biniem before plotting

#for i in range(0,len(plot_list),2):
#
#    
#    
#    if type(DataTable1[plot_list[i]][0]) == bool:
#        plt.hist(np.append(DATA_00[plot_list[i]],DATA_11[plot_list[i+1]]),bins=2,  alpha=0.9, histtype=u'step',label='tags', color='b')
#        plt.hist(np.append(DATA_00[plot_list[i+1]],DATA_11[plot_list[i]]),bins=2,  alpha=0.9, histtype=u'step',label='probes', color='r') #u'step' normed=True
#        plt.legend(loc='best')
#        plt.xlabel("{0}".format(plot_list[i]), size = 15)
#        plt.grid()
#        plt.ylabel("Event count", size = 10)
#        plt.title("{0} comparison".format(plot_list[i]), size = 15)
            
            #plt.xticks(np.arange(0,6,0.5), size = 15)
#        plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_bin_comparison_tagWp80_10pT/{0}/{1}_{2}_pT/TnP_comp_{3}_pT.png".format(section[l],pt_low[k],pt_high[k],plot_list[i]))
#        plt.show()
#    else:
#        plt.hist(np.append(DATA_00[plot_list[i]],DATA_11[plot_list[i+1]]),bins=100,  alpha=0.9, histtype=u'step',label='tags', color='b')
#        plt.hist(np.append(DATA_00[plot_list[i+1]],DATA_11[plot_list[i]]),bins=100,  alpha=0.9, histtype=u'step',label='probes', color='r') #u'step' normed=True
#        plt.legend(loc='best')
#        plt.xlabel("{0}".format(plot_list[i]), size = 15)
#        plt.grid()
#        plt.ylabel("Event count", size = 10)
#        plt.title("{0} comparison".format(plot_list[i]), size = 15)
#        plt.xlim([np.percentile(DATA_00[plot_list[i]],1),np.percentile(DATA_00[plot_list[i]],99)])
        #plt.xticks(np.arange(0,6,0.5), size = 15)
#        plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_bin_comparison_tagWp80_10pT/{0}/{1}_{2}_pT/TnP_comp_{3}_pT.png".format(section[l],pt_low[k],pt_high[k],plot_list[i]))
#        plt.show()
            
    #print("{0} did not want to plot, shady boolean operator in plot function?".format(plot_list[i]))
############################################################################################
### Tags vs Probes parametru plotosanas dala
        
##for l in range(len(eta_low)): # loop over barrel and encap ## does not plot electron ID's
##    for k in range(len(pt_low)):
##          
##        cut_lowpt_0 = DATA_0["ele_pt_1"] > pt_low[k]
##        cut_highpt_0 = DATA_0["ele_pt_1"] < pt_high[k]
##        cut_loweta_0 = abs(DATA_0["ele_eta_1"]) > eta_low[l]
##        cut_higheta_0 = abs(DATA_0["ele_eta_1"]) < eta_high[l]
##
##        cut_lowpt_1 = DATA_1["ele_pt_0"] > pt_low[k]
##        cut_highpt_1 = DATA_1["ele_pt_0"] < pt_high[k]
##        cut_loweta_1 = abs(DATA_1["ele_eta_0"]) > eta_low[l]
##        cut_higheta_1 = abs(DATA_1["ele_eta_0"]) < eta_high[l]
##        
##        data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0
##        data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1
##
##        DATA_00 = DATA_0[data_cut_0]
##        DATA_11 = DATA_1[data_cut_1]
##
##
##        for i in range(0,len(plot_list),2):
##
##            try:
##    
##    
##                if type(DataTable1[plot_list[i]][0]) == bool:
##                    plt.hist(ak.concatenate([DATA_00[plot_list[i]],DATA_11[plot_list[i+1]]],axis=0),bins=2,  alpha=0.9, histtype=u'step',label='tags', color='b')
##                    plt.hist(ak.concatenate([DATA_00[plot_list[i+1]],DATA_11[plot_list[i]]],axis=0),bins=2,  alpha=0.9, histtype=u'step',label='probes', color='r') #u'step' normed=True
##                    plt.legend(loc='best')
##                    plt.xlabel("{0}".format(plot_list[i]), size = 15)
##                    plt.grid()
##                    plt.ylabel("Event count", size = 10)
##                    plt.title("{0} comparison".format(plot_list[i]), size = 15)
##            
##                    #plt.xticks(np.arange(0,6,0.5), size = 15)
##                    plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_bin_comparison_tagWp80_10pT/{0}_MC/{1}_{2}_pT/TnP_comp_{3}_MC.png".format(section[l],pt_low[k],pt_high[k],plot_list[i]))
##                    #plt.show()
##                    plt.clf()
##                else:
##                    plt.hist(ak.concatenate([DATA_00[plot_list[i]],DATA_11[plot_list[i+1]]],axis=0),bins=100, range=[0.9*min(ak.concatenate([DATA_00[plot_list[i]],DATA_11[plot_list[i+1]]],axis=0)),1.1*max(ak.concatenate([DATA_00[plot_list[i]],DATA_11[plot_list[i+1]]],axis=0))], alpha=0.9, histtype=u'step',label='tags', color='b')
##                    plt.hist(ak.concatenate([DATA_00[plot_list[i+1]],DATA_11[plot_list[i]]],axis=0),bins=100, range=[0.9*min(ak.concatenate([DATA_00[plot_list[i+1]],DATA_11[plot_list[i]]],axis=0)),1.1*max(ak.concatenate([DATA_00[plot_list[i+1]],DATA_11[plot_list[i]]],axis=0))], alpha=0.9, histtype=u'step',label='probes', color='r') #u'step' normed=True
##                    plt.legend(loc='best')
##                    plt.xlabel("{0}".format(plot_list[i]), size = 15)
##                    plt.grid()
##                    plt.ylabel("Event count", size = 10)
##                    plt.title("{0} comparison".format(plot_list[i]), size = 15)
##                    #plt.xlim([np.percentile(DATA_00[plot_list[i]],1),np.percentile(DATA_00[plot_list[i]],99)])
##                    #plt.xticks(np.arange(0,6,0.5), size = 15)
##                    plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_bin_comparison_tagWp80_10pT/{0}_MC/{1}_{2}_pT/TnP_comp_{3}_MC.png".format(section[l],pt_low[k],pt_high[k],plot_list[i]))
##                    #plt.show()
##                    plt.clf()
##                               
##            except:
##                print("{0} did not want to plot, shady boolean operator in plot function?".format(plot_list[i]))
##                continue
##                
#############################################################################################

efficiency_barrel = []
efficiency_endcap = []


for l in range(len(eta_low)): # loop over barrel and encap ## does not plot electron ID's
    for k in range(len(pt_low)):
          
        cut_lowpt_0 = DATA_0["ele_pt_1"] > pt_low[k]
        cut_highpt_0 = DATA_0["ele_pt_1"] < pt_high[k]
        cut_loweta_0 = abs(DATA_0["ele_eta_1"]) > eta_low[l]
        cut_higheta_0 = abs(DATA_0["ele_eta_1"]) < eta_high[l]

        cut_lowpt_1 = DATA_1["ele_pt_0"] > pt_low[k]
        cut_highpt_1 = DATA_1["ele_pt_0"] < pt_high[k]
        cut_loweta_1 = abs(DATA_1["ele_eta_0"]) > eta_low[l]
        cut_higheta_1 = abs(DATA_1["ele_eta_0"]) < eta_high[l]
        
        data_cut_0 = cut_lowpt_0&cut_highpt_0&cut_loweta_0&cut_higheta_0
        data_cut_1 = cut_lowpt_1&cut_highpt_1&cut_loweta_1&cut_higheta_1

        DATA_00 = DATA_0[data_cut_0]
        DATA_11 = DATA_1[data_cut_1]
        
        #Name the probe ID, for which we want to calculate effieciency
        probe_ID_0 = DATA_00["mvaEleID_Fall17_iso_V2_wp80_1"] == 1 # for passing
        probe_ID_1 = DATA_11["mvaEleID_Fall17_iso_V2_wp80_0"] == 1
        
        probe_ID_00 = DATA_00["mvaEleID_Fall17_iso_V2_wp80_1"] == 0 # for failing
        probe_ID_11 = DATA_11["mvaEleID_Fall17_iso_V2_wp80_0"] == 0
        
        DATA_000 = DATA_00[probe_ID_0] # passing
        DATA_111 = DATA_11[probe_ID_1] # passing
        
        DATA_0000 = DATA_00[probe_ID_00] # failing
        DATA_1111 = DATA_11[probe_ID_11] # failing
        

       # try:
        
        #yoo=(0.5,0.6,0.5,0.5,3.1,0.1,300,5,-10,-200,200)
        yoo=(0.5,0.6,0.5,0.5,3.1,0.1,300)
        y_qmc, edges_mc = np.histogram(ak.concatenate([DATA_000["Diele_mass_0"],DATA_111["Diele_mass_0"]],axis=0), bins=50, range=(2, 4)) # weights=df_mc["weights"]
        x_mc = (edges_mc[1:] + edges_mc[:-1])/2 
        y_mc=y_qmc
        yerr_mc=np.sqrt(y_mc) 
        s_mc_y=y_mc
        s_mc_x=x_mc 
        #parameters_mc, err_mc = opt.curve_fit(cBall_exp, s_mc_x, s_mc_y, p0=yoo ) #method = 'dogbox'
        #x_fit=np.linspace(2,4,300)
            
        #plt.figure()
            
        #plt.bar(s_mc_x,s_mc_y,color="k",   alpha=0.1, label="signal in mc", width=(s_mc_x[1]-s_mc_x[0])) #align='center'  # signal
            
        #plt.plot(x_fit,cBall_exp(x_fit-(edges_mc[1]-edges_mc[0])/2,*parameters_mc),label="DSCB\n a_l={0:.2f} ({1:.2f}),\n a_h={2:.2f} ({3:.2f}),\n n_l={4:.2f} ({5:.2f}),\n n_h={6:.2f} ({7:.2f}),\n mean={8:.2f} ({9:.2f}),\n sigma={10:.2f} ({11:.2f}),\n N={12:.2f} ({13:.2f})".format(parameters_mc[0],err_mc[0][0]**0.5,parameters_mc[1],err_mc[1][1]**0.5,parameters_mc[2],err_mc[2][2]**0.5,parameters_mc[3],err_mc[3][3]**0.5,parameters_mc[4],err_mc[4][4]**0.5,parameters_mc[5],err_mc[5][5]**0.5,parameters_mc[6],err_mc[6][6]**0.5)) #label="c"  # signal fit

        #plt.plot(x_fit,linear(x_fit,parameters_mc[-4],parameters_mc[-3],parameters_mc[-2],parameters_mc[-1]),label = 'bck') # background fit
            
        #plt.show()
        #plt.clf()            
            
        plt.hist(ak.concatenate([DATA_000["Diele_mass_0"],DATA_111["Diele_mass_0"]],axis=0),bins=100, alpha=0.9, histtype=u'step',label='passing', color='b')
        plt.legend(loc='best')
        plt.xlabel("Inv.Mass", size = 15)
        plt.grid()
        plt.ylabel("Event count", size = 10)
        plt.title("Passing probes {0} {1}_{2}".format(section[l],pt_low[k],pt_high[k]), size = 15)
            #plt.xlim([np.percentile(DATA_00[plot_list[i]],1),np.percentile(DATA_00[plot_list[i]],99)])
            #plt.xticks(np.arange(0,6,0.5), size = 15)
        #plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_JPsi_Plots_tag_7pT_wp80_noiso/iso_wp80/MC/MC_{0}/TnP_Passing_{1}_{2}_MC.png".format(section[l],pt_low[k],pt_high[k]))
        #plt.show()
        plt.clf()
        
        yoo_1=(0.6,0.7,1.5,1.5,3.1,0.1,200)
        y_qmc_1, edges_mc_1 = np.histogram(ak.concatenate([DATA_0000["Diele_mass_0"],DATA_1111["Diele_mass_0"]],axis=0), bins=50, range=(2, 4)) # weights=df_mc["weights"]
        x_mc_1 = (edges_mc_1[1:] + edges_mc_1[:-1])/2 
        y_mc_1=y_qmc_1
        yerr_mc_1=np.sqrt(y_mc_1) 
        s_mc_y_1=y_mc_1
        s_mc_x_1=x_mc_1 
        #parameters_mc_1, err_mc_1 = opt.curve_fit(cBall_exp, s_mc_x_1, s_mc_y_1, p0=yoo_1, maxfev = 5000 ) #method = 'dogbox'
        #x_fit_1=np.linspace(2,4,300)
            
        #plt.figure()
            
        #plt.bar(s_mc_x_1,s_mc_y_1,color="k",   alpha=0.1, label="signal in mc", width=(s_mc_x_1[1]-s_mc_x_1[0])) #align='center'  # signal
            
        #plt.plot(x_fit_1,cBall_exp(x_fit_1-(edges_mc_1[1]-edges_mc_1[0])/2,*parameters_mc_1),label="DSCB\n a_l={0:.2f} ({1:.2f}),\n a_h={2:.2f} ({3:.2f}),\n n_l={4:.2f} ({5:.2f}),\n n_h={6:.2f} ({7:.2f}),\n mean={8:.2f} ({9:.2f}),\n sigma={10:.2f} ({11:.2f}),\n N={12:.2f} ({13:.2f})".format(parameters_mc_1[0],err_mc_1[0][0]**0.5,parameters_mc_1[1],err_mc_1[1][1]**0.5,parameters_mc_1[2],err_mc_1[2][2]**0.5,parameters_mc_1[3],err_mc_1[3][3]**0.5,parameters_mc_1[4],err_mc_1[4][4]**0.5,parameters_mc_1[5],err_mc_1[5][5]**0.5,parameters_mc_1[6],err_mc_1[6][6]**0.5)) #label="c"  # signal fit

        #plt.plot(x_fit,linear(x_fit,parameters_mc[-4],parameters_mc[-3],parameters_mc[-2],parameters_mc[-1]),label = 'bck') # background fit        
        
        
        
        plt.hist(ak.concatenate([DATA_0000["Diele_mass_0"],DATA_1111["Diele_mass_0"]],axis=0),bins=100, alpha=0.9, histtype=u'step',label='failing', color='r')
        plt.legend(loc='best')
        plt.xlabel("Inv.Mass", size = 15)
        plt.grid()
        plt.ylabel("Event count", size = 10)
        plt.title("Failing probes {0} {1}_{2}".format(section[l],pt_low[k],pt_high[k]), size = 15)
            #plt.xlim([np.percentile(DATA_00[plot_list[i]],1),np.percentile(DATA_00[plot_list[i]],99)])
            #plt.xticks(np.arange(0,6,0.5), size = 15)
        #plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_JPsi_Plots_tag_7pT_wp80_noiso/iso_wp80/MC/MC_{0}/TnP_Failing_{1}_{2}_MC.png".format(section[l],pt_low[k],pt_high[k]))
        #plt.show()
        plt.clf()
        
       # except:
       #     print("Could not plot {0} {1}_{2}".format(section[l],pt_low[k],pt_high[k]))
       #     continue
        if l == 0:
            efficiency_barrel.append(float(sum(y_qmc))/float((sum(y_qmc)+sum(y_qmc_1))))
        else:
            efficiency_endcap.append(float(sum(y_qmc))/float((sum(y_qmc)+sum(y_qmc_1))))
                
plt.errorbar([5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,17.5,25],efficiency_barrel,xerr=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,2.5,5],label='MC_barrel_JPsi', color='b', marker = 'o', linestyle = '') #fmt = None
plt.errorbar([5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,17.5,25],efficiency_endcap,xerr=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,2.5,5],label='MC_endcap_JPsi', color='r', marker = 'o', linestyle = '')

plt.errorbar([6,8.5,12.5,17.5,25,35],[0.385,0.486,0.585,0.620,0.685,0.766],xerr=[1,1.5,2.5,2.5,5,5],yerr=[0.015,0.007,0.003,0.002,0.001,0.001],label='MC_barrel_Z', color='g', marker = 'o', linestyle = '') # iso wp80
plt.errorbar([6,8.5,12.5,17.5,25,35],[0.455,0.589,0.639,0.682,0.742,0.811],xerr=[1,1.5,2.5,2.5,5,5],yerr = [0.011,0.006,0.003,0.002,0.001,0.001],label='MC_endcap_Z', color='k', marker = 'o', linestyle = '') # iso wp80

#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.733,0.8145,0.833,0.838,0.866,0.909],xerr=[1,1.5,2.5,2.5,5,5],yerr=[0.009,0.005,0.003,0.002,0.001,0.001],label='MC_endcap_Z', color='k', marker = 'o', linestyle = '') # iso wp90
#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.667,0.738,0.781,0.788,0.826,0.884],xerr=[1,1.5,2.5,2.5,5,5],yerr = [0.014,0.006,0.002,0.002,0.001,0.001],label='MC_barrel_Z', color='g', marker = 'o', linestyle = '') # iso wp90

#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.676,0.846,0.814,0.903,0.959,0.984],xerr=[1,1.5,2.5,2.5,5,5],yerr=[0.01,0.004,0.003,0.001,0.001,0.001],label='MC_endcap_Z', color='k', marker = 'o', linestyle = '') # iso wp98 HZZ
#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.702,0.850,0.811,0.894,0.952,0.983],xerr=[1,1.5,2.5,2.5,5,5],yerr = [0.014,0.005,0.002,0.001,0.001,0.001],label='MC_barrel_Z', color='g', marker = 'o', linestyle = '') # iso wp98 HZZ

#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.612,0.815,0.886,0.951,0.981,0.993],xerr=[1,1.5,2.5,2.5,5,5],yerr=[0.010,0.005,0.002,0.001,0.001,0.001],label='MC_endcap_Z', color='k', marker = 'o', linestyle = '') # noiso wp98 wpL
#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.744,0.881,0.948,0.977,0.990,0.997],xerr=[1,1.5,2.5,2.5,5,5],yerr = [0.013,0.004,0.001,0.001,0.001,0.001],label='MC_barrel_Z', color='g', marker = 'o', linestyle = '') # noiso wp98 wpL

#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.727,0.816,0.806,0.824,0.862,0.906],xerr=[1,1.5,2.5,2.5,5,5],yerr=[0.009,0.005,0.003,0.002,0.001,0.001],label='MC_endcap_Z', color='k', marker = 'o', linestyle = '') # noiso wp90 HZZ
#plt.errorbar([6,8.5,12.5,17.5,25,35],[0.673,0.745,0.782,0.802,0.838,0.885],xerr=[1,1.5,2.5,2.5,5,5],yerr = [0.014,0.006,0.002,0.001,0.001,0.001],label='MC_barrel_Z', color='g', marker = 'o', linestyle = '') # noiso wp90 HZZ

plt.legend(loc='best')
plt.xlabel("pT, GeV", size = 15)
plt.grid()
plt.ylabel("Efficiency", size = 15)
plt.title("MC eff, MVA_Iso_wp80 ", size = 20)
plt.xlim(4,31)
plt.savefig("/eos/user/n/nstrautn/www/TnP_plots/TnP_JPsi_Plots_tag_10pT_wp80_noiso/iso_wp80/MC/JPsi_MC_eff_comparison_Iso_wp80.png")
plt.show()
