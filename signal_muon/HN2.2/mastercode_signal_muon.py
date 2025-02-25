import ROOT
from ROOT import TFile, gDirectory, TH1D, TH2D, TCanvas, TLorentzVector
import random
import math
from math import sqrt, pi, log10
from tqdm import tqdm
import os
import sys
import numpy as np
from copy import deepcopy
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1,0,1)

tau_decay_BR = 0.1739 #tau decay to muon and two neutrinos
tau_energy_thr = 0 #30 
muon_pt_thr = 1.3 #GeV
FD_res = 0.35 #mm
FD_thr = 0.7 #mm
z_thr = 15000 #mm

nBc = 5e10 # number of Bc per year

cutflow_labels = ["generated","LHCb geo","#mu_{1} trig","#mu_{2} trig", "#mu_{3} trig", "J/#psi mass", "B_{c} mass"]


real_ctau = {2.2:4628.421, 3.0:879.474, 4.0:161.646, 4.5:24.8, 5.0:41.739}

#Mass of HN [GeV]
M = float(sys.argv[1])

massPoint = int(sys.argv[2])

nInclusive = int(sys.argv[3])


ctau = float(sys.argv[4])

  

infile = TFile( '/eos/user/m/mmourad/HN_simulation/trees/signal_muon_mass'+sys.argv[1]+'_tree.root' )
tree = gDirectory.Get('DecayTree')
entries = tree.GetEntries()

# define Lorentz vectors and variables

mu1 = TLorentzVector()
mu2 = TLorentzVector()
mu3 = TLorentzVector()
tau = TLorentzVector()
neutrino = TLorentzVector()
B = TLorentzVector()
HN = TLorentzVector()

#define histograms

h_cutflow = TH1D("h_cutflow","#mu#mu#mu",len(cutflow_labels),0,len(cutflow_labels))
for i in range(len(cutflow_labels)):
  h_cutflow.GetXaxis().SetBinLabel(i+1,cutflow_labels[i])
h_cutflow.SetBinContent(1,nInclusive)

h_trigger_efficiency = TH1D("h_trigger_efficiency","trigger efficiency;muon p_{T} threshold [GeV];trigger efficiency",20,0,2)

h_FDeff_ctau = TH1D("h_FDeff_ctau","flight distance efficiency;c#tau [mm];flight distance efficiency",100,0,1)
h_FDeff_coupling2 = TH1D("h_FDeff_coupling2","FD efficiency vs coupling;-log(#alpha^{2});FD efficiency",80,0,8)

h_FDeff_times_coupling4_ctau = TH1D("h_FDeff_times_coupling4_ctau","FD efficiency #times #alpha^{4};c#tau [mm];FD efficiency #times #alpha^{4}",100,0,1)
h_FDeff_times_coupling4_coupling2 = TH1D("h_FDeff_times_coupling4_coupling2","FD efficiency #times #alpha^{4};-log(#alpha^{2});FD efficiency #times #alpha^{4}",80,0,8)


h_BR_Maj_coupling2 = TH1D("h_BR_Maj_coupling2","BR(Majorana);-log(#alpha^{2});BR(Majorana)",80,0,8)
h_BR_Dir_coupling2 = TH1D("h_BR_Dir_coupling2","BR(Dirac);-log(#alpha^{2});BR(Dirac)",80,0,8)

h_BR_Maj_couplingMuTau = TH2D("h_BR_Maj_couplingMuTau","BR(Majorana);-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)
h_BR_Dir_couplingMuTau = TH2D("h_BR_Dir_couplingMuTau","BR(Dirac);-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)

#h_GN_Maj_couplingMuTau = TH2D("h_GN_Maj_couplingMuTau","GN(Majorana);-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)
#h_GN_Dir_couplingMuTau = TH2D("h_GN_Dir_couplingMuTau","GN(Dirac);-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)

h_ctau_Maj_couplingMuTau = TH2D("h_ctau_Maj_couplingMuTau","c#tau(Majorana);-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)
h_ctau_Dir_couplingMuTau = TH2D("h_ctau_Dir_couplingMuTau","c#tau(Dirac);-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)


h_BR_ctau = TH1D("h_BR_ctau","BR(Majorana) for #alpha^{2}=;c#tau [mm];BR(Majorana)",100,0,20000)
h_BR_Maj_eff_coupling2 = TH1D("h_BR_Maj_eff_coupling2","Effective BR(Majorana);-log(#alpha^{2});Effective BR(Majorana)",80,0,8)
h_BR_Dir_eff_coupling2 = TH1D("h_BR_Dir_eff_coupling2","Effective BR(Dirac);-log(#alpha^{2});Effective BR(Dirac)",80,0,8)

h_nEvents_Maj_coupling2 = TH1D("h_nEvents_Maj_coupling2","N_{exp}^{Maj.} per year in HL-LHC;-log(#alpha^{2});N_{exp}^{Maj.} per year in HL-LHC",80,0,8)
h_nEvents_Dir_coupling2 = TH1D("h_nEvents_Dir_coupling2","N_{exp}^{Dir.} per year in HL-LHC;-log(#alpha^{2});N_{exp}^{Dir.} per year in HL-LHC",80,0,8)

h_nEvents_Maj_couplingMuTau = TH2D("h_nEvents_Maj_couplingMuTau","N_{exp}^{Maj.} per year in HL-LHC;-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)
h_nEvents_Dir_couplingMuTau = TH2D("h_nEvents_Dir_couplingMuTau","N_{exp}^{Dir.} per year in HL-LHC;-log(#alpha^{2}_{#mu});-log(#alpha^{2}_{#tau})",80,0,8,80,0,8)


h_FDx_HN = TH1D("h_FDx_HN"," Flight distance of HN in x; Flight distance of HN in x [mm]",50,0,4000)
h_FDy_HN = TH1D("h_FDy_HN"," Flight distance of HN in y; Flight distance of HN in y [mm]",50,0,4000)
h_FDt_HN = TH1D("h_FDt_HN","Transverse  Flight distance of HN;Transverse  Flight distance of HN[mm]",50,0,4000)
h_FDz_HN = TH1D("h_FDz_HN"," Flight distance of HN in z; Flight distance of HN in z [m]",20,0,10)
h_FD_HN = TH1D("h_FD_HN"," Flight distance of HN; Flight distance of HN [m]",20,0,10)
h_FD_HN_norm = TH1D("h_FD_HN_norm","Flight distance of HN in its RF / c#tau; Flight distance of HN in its RF / c#tau",25,0,1.5)

h_FD_mu1_mu3 = TH1D("h_FD_mu1_mu3","distance of #mu_{1} and #mu_{3};distance of #mu_{1} and #mu_{3} [mm]",100,0,100)


h_pt_mu1 = TH1D("h_pt_mu1","p_{T} of #mu_{1};p_{T} [GeV]",20,1.3,21.3)
h_pt_mu2 = TH1D("h_pt_mu2","p_{T} of #mu_{2};p_{T} [GeV]",20,1.3,21.3)
h_pt_mu3 = TH1D("h_pt_mu3","p_{T} of #mu_{3};p_{T} [GeV]",20,1.3,21.3)
h_pt_mu1mu2 = TH1D("h_pt_mu1mu2","p_{T} of #mu_{1}+#mu_{2};p_{T} [GeV]",30,0,30)
h_pt_mu1mu3 = TH1D("h_pt_mu1mu3","p_{T} of #mu_{1}+#mu_{3};p_{T} [GeV]",30,0,30)
h_pt_mu2mu3 = TH1D("h_pt_mu2mu3","p_{T} of #mu_{2}+#mu_{3};p_{T} [GeV]",30,0,30)
h_pt_HN = TH1D("h_pt_HN","p_{T} of HN;p_{T} [GeV]",40,0,40)
h_pt_B = TH1D("h_pt_B","p_{T} of B_{c};p_{T} [GeV]",40,0,40)
h_pt_neutrino = TH1D("h_pt_neutrino","p_{T} of #nu_{#mu3};p_{T} [GeV]",30,0,30)

h_pt_tau = TH1D("h_pt_tau","p_{T} of #tau;p_{T} [GeV]",40,0,40)

h_phi_mu1 = TH1D("h_phi_mu1","#phi of #mu_{1};#phi",50,-3.1415,3.1415)
h_phi_mu2 = TH1D("h_phi_mu2","#phi of #mu_{2};#phi",50,-3.1415,3.1415)
h_phi_mu3 = TH1D("h_phi_mu3","#phi of #mu_{3};#phi",50,-3.1415,3.1415)
h_phi_tau = TH1D("h_phi_tau","#phi of #tau;#phi",50,-3.1415,3.1415)
h_phi_neutrino = TH1D("h_phi_neutrino","#phi of #nu;#phi",50,-3.1415,3.1415)

h_deltaphi_mu1tau = TH1D("h_deltaphi_mu1tau","#Delta#phi #mu_{1} and #tau;#Delta#phi",30,0,1)
h_deltaphi_mu2tau = TH1D("h_deltaphi_mu2tau","#Delta#phi #mu_{2} and #tau;#Delta#phi",30,0,1)

h_deltaR_mu2tau = TH1D("h_deltaR_mu2tau","#Delta R between #mu_{2} and #tau;#Delta R", 45,0,1.5)

h_deltaphi_mu1mu2 = TH1D("h_deltaphi_mu1mu2","#Delta#phi #mu_{1} and #mu_{2};#Delta#phi",30,0,1)
h_deltaphi_mu1mu3 = TH1D("h_deltaphi_mu1mu3","#Delta#phi #mu_{1} and #mu_{3};#Delta#phi",30,0,1)
h_deltaphi_mu2mu3 = TH1D("h_deltaphi_mu2mu3","#Delta#phi #mu_{2} and #mu_{3};#Delta#phi",30,0,1)

h_deltaeta_mu1mu2 = TH1D("h_deltaeta_mu1mu2","#Delta#eta #mu_{1} and #mu_{2};#Delta#eta",30,0,1)
h_deltaeta_mu1mu3 = TH1D("h_deltaeta_mu1mu3","#Delta#eta #mu_{1} and #mu_{3};#Delta#eta",30,0,1)
h_deltaeta_mu2mu3 = TH1D("h_deltaeta_mu2mu3","#Delta#eta #mu_{2} and #mu_{3};#Delta#eta",30,0,1)


h_deltaR_mu1mu2 = TH1D("h_deltaR_mu1mu2","#Delta R between #mu_{1} and #mu_{2};#Delta R", 45,0,1.5)
h_deltaR_mu1mu3 = TH1D("h_deltaR_mu1mu3","#Delta R between #mu_{1} and #mu_{3};#Delta R", 45,0,1.5)
h_deltaR_mu2mu3 = TH1D("h_deltaR_mu2mu3","#Delta R between #mu_{2} and #mu_{3};#Delta R", 45,0,1.5)
#delta R is square root of delta phi squared + delta eta squared, 3d Lorentz invariant angle

h_theta_HN = TH1D("h_theta_HN","angle between visible HN and #gamma(visible B);#theta(visible HN, #gamma(visible B))", 40,0,pi)
h_true_theta_HN = TH1D("h_true_theta_HN","angle between true HN and #gamma(true B);#theta(true HN, #gamma(true B))", 40,0,pi)

h_energy_HN = TH1D("h_energy_HN","Energy of HN;E_{HN} [GeV]",50,0,1000)
h_energy_B = TH1D("h_energy_B","Energy of B_{c};E_{B} [GeV]",50,0,1000)
h_energy_tau = TH1D("h_energy_tau","Energy of #tau;E_{#tau} [GeV]",50, tau_energy_thr,1000)
h_energy_tau_HNRF = TH1D("h_energy_tau_HNRF","Energy of #tau in HN RF;E_{#tau} [GeV]",75, 1.5,3)
h_energy_mu3_HNRF = TH1D("h_energy_mu3_HNRF","Energy of #mu_{3} in visible HN RF;E_{#mu_{3}} [GeV]",30, 0,3)
h_energy_mu1 = TH1D("h_energy_mu1","Energy of #mu_{1};E_{#mu_{1}} [GeV]",25, 0,500)
h_energy_mu2 = TH1D("h_energy_mu2","Energy of #mu_{2};E_{#mu_{2}} [GeV]",25, 0,500)
h_energy_mu3 = TH1D("h_energy_mu3","Energy of #mu_{3};E_{#mu_{3}} [GeV]",25, 0,500)

#Dalitz plot for HN
h_dalitz_HN = TH2D("h_dalitz_HN ","Dalitzplot for B_{c}^{+} -> #mu+#mu+#mu-;m^{2}_{#mu_{1}#mu_{3}} [GeV];m^{2}_{#mu_{2}#mu_{3}} [GeV]" ,45 ,0 ,15, 45, 0 , 15)
#h_m_HNprime = TH1D("h_m_HNprime", "Generated Mass of HN; M [GeV]", 150, M - 0.1 , M + 0.1)
h_m_HN = TH1D("h_m_HN", "Visible Mass of HN(without neutrino) ; M [GeV]", 60, 0, 6)
h_m_HNprime = TH1D("h_m_HNprime", "Full Invariant Mass of HN; M [GeV]", 240, 0 , 8)

h_m_mu1mu2 = TH1D("h_m_mu1mu2", "invariant mass of #mu_{1} and #mu_{2}; invariant mass of #mu_{1} and #mu_{2} [GeV]", 40, 0, 6)
h_m_mu1mu3 = TH1D("h_m_mu1mu3", "invariant mass of #mu_{1} and #mu_{3}; invariant mass of #mu_{1} and #mu_{3} [GeV]", 40, 0, 6)
h_m_mu2mu3 = TH1D("h_m_mu2mu3", "invariant mass of #mu_{2} and #mu_{3}; invariant mass of #mu_{2} and #mu_{3} [GeV]", 60, 0, 6)

h_m_Jpsi = TH1D("h_m_Jpsi", "Invariant mass of #mu_{2} and #mu_{3}; invariant mass of #mu_{2} and #mu_{3} [GeV]", 100, 2.8 , 3.4)

h_m_B = TH1D("h_m_B", "Visible mass of B;visible mass of B [GeV]", 35, 0 , 7)

#2d plots of pt
h_pt_m1vsHN = TH2D("h_pt_m1vsHN","p_{T} of HN vs p_{T} of #mu_{1};p_{T} of HN [GeV]; p_{T} of #mu_{1} [GeV]", 40, 0, 100, 40, 0, 100)
h_pt_m1vssum = TH2D("h_pt_m1vssum","Sum of p_{T} of #mu_{2}, #tau, #bar#nu_{#tau} vs p_{T} of #mu_{1};p_{T} of sum [GeV]; p_{T} of #mu_{1} [GeV]", 40, 0, 100, 40, 0, 100)
h_pt_Bvsall = TH2D("h_pt_Bvsall","Sum of p_{T} of #mu_{1}, #mu_{2}, #tau, #bar#nu_{#tau} vs p_{T} of B_{c};p_{T} of sum [GeV]; p_{T} of B_{c} [GeV]", 40, 0, 100, 40, 0, 100)

#label axes of histograms
#either with semi-column method like above or:
#example: myHist.GetXaxis().SetTitle("theta")
#myHist.GetYaxis().SetTitle("K_{pT} [GeV]")



# get BR as a function of coupling squared
mass_index = int(10*(M-1.9))
#mass_index = massPoint
BR_Maj_mu = np.loadtxt('/afs/cern.ch/work/m/mmourad/master_arbeit/RapidSim/HN_simulation/signal_muon/BR_Majorana.txt', delimiter=',')
BR_Dir_mu = np.loadtxt('/afs/cern.ch/work/m/mmourad/master_arbeit/RapidSim/HN_simulation/signal_muon/BR_Dirac.txt', delimiter=',')
BR_Maj = np.loadtxt('/afs/cern.ch/work/m/mmourad/master_arbeit/RapidSim/HN_simulation/BR_calculation/BR_Majorana_mu_tau.txt', delimiter=',').reshape((81,81,42))
BR_Dir = np.loadtxt('/afs/cern.ch/work/m/mmourad/master_arbeit/RapidSim/HN_simulation/BR_calculation/BR_Dirac_mu_tau.txt', delimiter=',').reshape((81,81,42))
ctau_Maj = np.loadtxt('/afs/cern.ch/work/m/mmourad/master_arbeit/RapidSim/HN_simulation/BR_calculation/ctau_Majorana_mu_tau.txt', delimiter=',').reshape((81,81,42))
ctau_Dir = np.loadtxt('/afs/cern.ch/work/m/mmourad/master_arbeit/RapidSim/HN_simulation/BR_calculation/ctau_Dirac_mu_tau.txt', delimiter=',').reshape((81,81,42))


for bx in range(80):
  h_BR_Maj_coupling2.SetBinContent(bx+1,tau_decay_BR*BR_Maj_mu[bx][mass_index])
  h_BR_Dir_coupling2.SetBinContent(bx+1,tau_decay_BR*BR_Dir_mu[bx][mass_index])
  for by in range(80):
    h_BR_Maj_couplingMuTau.SetBinContent(bx+1,by+1,tau_decay_BR*BR_Maj[bx][by][mass_index])
    h_BR_Dir_couplingMuTau.SetBinContent(bx+1,by+1,tau_decay_BR*BR_Dir[bx][by][mass_index])
    #h_GN_Maj_couplingMuTau.SetBinContent(bx+1,by+1,GN_Maj[bx][by][mass_index])
    #h_GN_Dir_couplingMuTau.SetBinContent(bx+1,by+1,GN_Dir[bx][by][mass_index])
    h_ctau_Maj_couplingMuTau.SetBinContent(bx+1,by+1,ctau_Maj[bx][by][mass_index])
    h_ctau_Dir_couplingMuTau.SetBinContent(bx+1,by+1,ctau_Dir[bx][by][mass_index])

ctau_Maj /= ctau
ctau_Dir /= ctau

nEvents_Maj = tau_decay_BR*nBc*BR_Maj
nEvents_Dir = tau_decay_BR*nBc*BR_Dir



#loop through events
pbar = tqdm(total=entries, unit="")
accepted = 0
selected = 0

for event in range( entries ):
 pbar.update()
 tree.GetEntry(event)
 h_cutflow.Fill(1)
 
 for b in range(20):
   if (tree.mu1_PT > b/10. and tree.mu2_PT > b/10. and tree.mu3_PT > b/10.):
     h_trigger_efficiency.Fill(b/10.)
 
 
 
 
 #FDx_mu1_mu3 = tree.mu1_origX-tree.tau_origX
 #FDy_mu1_mu3 = tree.mu1_origY-tree.tau_origY
 #FDz_mu1_mu3 = tree.mu1_origZ-tree.tau_origZ
 
 FDx_mu1_mu3 = tree.mu1_origX-tree.mu3_origX
 FDy_mu1_mu3 = tree.mu1_origY-tree.mu3_origY
 FDz_mu1_mu3 = tree.mu1_origZ-tree.mu3_origZ
 FDt_mu1_mu3 = math.sqrt(FDx_mu1_mu3**2 + FDy_mu1_mu3**2)
 FD_mu1_mu3 = math.sqrt(FDx_mu1_mu3**2 + FDy_mu1_mu3**2 + FDz_mu1_mu3**2) # different from HN FD, because of the tau FD
 
 
 #print(tree.mu1_origX,tree.B_vtxX,tree.mu3_origX,tree.tau_origX,tree.HN_vtxX,tree.tau_vtxX)
 #print(tree.tau_vtxX-tree.tau_origX)
 
 
 
 
 # instead of regenerating a sample per each ctau, we scale the FD of our baseline ctau with a coefficient of new_ctau/baseline_ctau
 
 offset = random.gauss(0,FD_res) #mm
 for b in range(100):
   #FD = h_FDeff_ctau.GetBinCenter(b+1) # basically new_ctau
   FD = b/100. + 0.005 # basically new_ctau
   #reco_FD = abs((FD / ctau) * FD_mu1_mu3+offset)
   reco_FD = abs((FD / ctau) * FD_mu1_mu3 + offset)
   if (reco_FD > FD_thr):
     h_FDeff_ctau.Fill(FD)
     h_FDeff_times_coupling4_ctau.Fill(FD,(10E-5*ctau/FD)**2) #alpha^2 = (10E-5*ctau/FD)
 
 
 if (event%10 == 0):
   for b in range(80):
     #logAlpha2 = -1*h_FDeff_coupling2.GetBinCenter(b+1)
     logAlpha2 = -1*((b/10. +0.05))
     reco_FD = abs( pow(10,-5-logAlpha2) * FD_mu1_mu3 + offset) #mm     #alpha^2 = (10E-5*ctau/FD)
     if (reco_FD > FD_thr and reco_FD < z_thr):
       h_FDeff_coupling2.Fill(-1*logAlpha2)
       h_FDeff_times_coupling4_coupling2.Fill(-1*logAlpha2,pow(10,logAlpha2)**2)
       h_nEvents_Maj_coupling2.Fill(-1*logAlpha2,h_BR_Maj_coupling2.GetBinContent(b+1)*nBc)
       h_nEvents_Dir_coupling2.Fill(-1*logAlpha2,h_BR_Dir_coupling2.GetBinContent(b+1)*nBc)
       h_BR_Maj_eff_coupling2.Fill(-1*logAlpha2,h_BR_Maj_coupling2.GetBinContent(b+1))
       h_BR_Dir_eff_coupling2.Fill(-1*logAlpha2,h_BR_Dir_coupling2.GetBinContent(b+1))
     
     for by in range(80):
       reco_FD_Maj = abs( ctau_Maj[b][by][mass_index] * FD_mu1_mu3 + offset) #mm
       reco_FD_Dir = abs( ctau_Dir[b][by][mass_index] * FD_mu1_mu3 + offset) #mm
       if (reco_FD_Maj > FD_thr and reco_FD_Maj < z_thr):
         h_nEvents_Maj_couplingMuTau.Fill((b/10. +0.05),(by/10. +0.05),nEvents_Maj[b][by][mass_index])
       if (reco_FD_Dir > FD_thr and reco_FD_Dir < z_thr):
         h_nEvents_Dir_couplingMuTau.Fill((b/10. +0.05),(by/10. +0.05),nEvents_Dir[b][by][mass_index])
 
 
 if (M == 4.5):
   continue
 
 
 if (tree.mu1_PT < muon_pt_thr):
   continue
 h_cutflow.Fill(2)
 if (tree.mu2_PT < muon_pt_thr):
   continue
 h_cutflow.Fill(3)
 if (tree.mu3_PT < muon_pt_thr):
   continue
 h_cutflow.Fill(4)
 accepted += 1
 
 
 
 mu1.SetPtEtaPhiE(tree.mu1_PT,tree.mu1_eta, tree.mu1_phi, tree.mu1_E) 
 mu2.SetPtEtaPhiE(tree.mu2_PT,tree.mu2_eta, tree.mu2_phi, tree.mu2_E) 
 mu3.SetPtEtaPhiE(tree.mu3_PT,tree.mu3_eta, tree.mu3_phi, tree.mu3_E) 
 tau.SetPtEtaPhiE(tree.tau_PT,tree.tau_eta, tree.tau_phi, tree.tau_E)
 HN.SetPtEtaPhiE(tree.HN_PT,tree.HN_eta, tree.HN_phi, tree.HN_E)
 B.SetPtEtaPhiE(tree.B_PT,tree.B_eta, tree.B_phi, tree.B_E)
 neutrino.SetPtEtaPhiE(tree.neutrino_PT,tree.neutrino_eta, tree.neutrino_phi, tree.neutrino_E)
 #neutrino = B - mu1 - mu2 - tau
 
 
 
 
 #calculate flight distance of heavy neutrino from decay vertices of B and HN
 HN_FDx = tree.HN_vtxX-tree.B_vtxX
 HN_FDy = tree.HN_vtxY-tree.B_vtxY
 HN_FDz = tree.HN_vtxZ-tree.B_vtxZ
 HN_FDt = math.sqrt(HN_FDx**2 + HN_FDy**2)
 HN_FD = math.sqrt(HN_FDx**2 + HN_FDy**2 + HN_FDz**2) 
 
 
 
 visible_B_boost = (mu1+mu2+mu3).BoostVector()
 boosted_HN_BRF = deepcopy(mu2+mu3)
 boosted_HN_BRF.Boost(-visible_B_boost)
 
 true_B_boost = B.BoostVector()
 true_boosted_HN_BRF = deepcopy(HN)
 true_boosted_HN_BRF.Boost(-true_B_boost)
 
 true_HN_boost = HN.BoostVector()
 true_boosted_tau_HNRF = deepcopy(tau)
 true_boosted_tau_HNRF.Boost(-true_HN_boost)
 #print(true_boosted_tau_HNRF.E(),true_boosted_tau_HNRF.P())
 #print(B.Pz(),mu1.Pz(),HN.Pz(),tau.Pz(),mu2.Pz(),neutrino.Pz(),(tau+mu2+neutrino).M(),(mu1+HN).M())
 
 visible_HN_boost = (mu2+mu3).BoostVector()
 boosted_mu3_HNRF = deepcopy(mu3)
 boosted_mu3_HNRF.Boost(-visible_HN_boost)
 
 
 
 #fill histograms
 h_pt_B.Fill(B.Pt())
 h_energy_B.Fill(B.E())
 h_pt_HN.Fill(HN.Pt())
 h_energy_HN.Fill(HN.E())
 h_energy_tau.Fill(tau.E())
 h_energy_mu1.Fill(mu1.E())
 h_energy_mu2.Fill(mu2.E())
 h_energy_mu3.Fill(mu3.E())
 h_pt_mu1.Fill(mu1.Pt())
 h_pt_mu2.Fill(mu2.Pt())
 h_pt_mu1mu2.Fill((mu1+mu2).Pt())
 h_pt_mu1mu3.Fill((mu1+mu3).Pt())
 h_pt_mu2mu3.Fill((mu2+mu3).Pt())
 h_pt_mu3.Fill(mu3.Pt())
 h_pt_tau.Fill(tau.Pt())
 h_pt_neutrino.Fill(neutrino.Pt())
 h_phi_mu1.Fill(mu1.Phi())
 h_phi_mu2.Fill(mu2.Phi())
 h_phi_mu3.Fill(mu3.Phi())
 h_phi_tau.Fill(tau.Phi())
 h_phi_neutrino.Fill(neutrino.Phi())
 h_deltaphi_mu1tau.Fill(abs(mu1.DeltaPhi(tau)))
 h_deltaphi_mu2tau.Fill(abs(mu2.DeltaPhi(tau)))
 h_deltaR_mu2tau.Fill(abs(mu2.DeltaR(tau))) 
 
 h_deltaphi_mu1mu2.Fill(abs(mu1.DeltaPhi(mu2)))
 h_deltaphi_mu1mu3.Fill(abs(mu1.DeltaPhi(mu3)))
 h_deltaphi_mu2mu3.Fill(abs(mu2.DeltaPhi(mu3)))
 h_deltaeta_mu1mu2.Fill(abs(mu1.Eta()-mu2.Eta()))
 h_deltaeta_mu1mu3.Fill(abs(mu1.Eta()-mu3.Eta()))
 h_deltaeta_mu2mu3.Fill(abs(mu2.Eta()-mu3.Eta()))
 h_deltaR_mu1mu2.Fill(abs(mu1.DeltaR(mu2))) 
 h_deltaR_mu1mu3.Fill(abs(mu1.DeltaR(mu3))) 
 h_deltaR_mu2mu3.Fill(abs(mu2.DeltaR(mu3))) 
 
 h_theta_HN.Fill(abs(boosted_HN_BRF.Angle((mu1+mu2+mu3).Vect())))
 h_true_theta_HN.Fill(abs(true_boosted_HN_BRF.Angle(B.Vect())))
 h_energy_tau_HNRF.Fill(true_boosted_tau_HNRF.E())
 h_energy_mu3_HNRF.Fill(boosted_mu3_HNRF.E())
 
 h_FDx_HN.Fill(HN_FDx*real_ctau[M]/ctau)
 h_FDy_HN.Fill(HN_FDy*real_ctau[M]/ctau)
 h_FDt_HN.Fill(HN_FDt*real_ctau[M]/ctau)
 h_FDz_HN.Fill(0.001*HN_FDz*real_ctau[M]/ctau)
 h_FD_HN.Fill(0.001*HN_FD*real_ctau[M]/ctau)
 h_FD_mu1_mu3.Fill(FD_mu1_mu3*real_ctau[M]/ctau)
 
 
 #caculate boost factor of HN 
 gamma_HN = HN.E() / HN.M()
 #invariant mass
 msquared_mu1tau = (mu1 + tau).M() **2
 msquared_mu2tau = (mu2 + tau).M() **2
 msquared_mu1mu3 = (mu1 + mu3).M() **2
 msquared_mu2mu3 = (mu2 + mu3).M() **2
 #m_HN = (mu2 + tau).M()
 m_HN = (mu2 + mu3).M()
 m_HNprime = (mu2 + tau + neutrino).M()
 
 #fill histogram with FD of HN in its rest frame
 h_FD_HN_norm.Fill(HN_FD/(gamma_HN*ctau))
 
 h_dalitz_HN.Fill(msquared_mu1mu3, msquared_mu2mu3)
 
 h_m_HN.Fill(m_HN)
 h_m_HNprime.Fill(m_HNprime)
 h_m_Jpsi.Fill((mu1+mu2).M())
 h_m_mu1mu2.Fill((mu1+mu2).M())
 h_m_mu1mu3.Fill((mu1+mu3).M())
 h_m_mu2mu3.Fill((mu2+mu3).M())
 h_m_B.Fill((mu1+mu2+mu3).M())
 #2d plots
 h_pt_m1vsHN.Fill(HN.Pt(), mu1.Pt())
 h_pt_m1vssum.Fill((mu2 + tau + neutrino).Pt(), mu1.Pt())
 h_pt_Bvsall.Fill((mu2 + mu1 + tau + neutrino).Pt(), B.Pt())
 
 
 if ((mu1+mu2).M() > 3.0 and (mu1+mu2).M() < 3.2):
   continue
 h_cutflow.Fill(5)
 if ((mu1+mu2+mu3).M() > 5.8):
   continue
 h_cutflow.Fill(6)
 
 selected += 1
 
pbar.close()

print("geo_accepted: ", entries)
print("analysis_accepted: ", accepted)
print("analysis selection efficiency: ", accepted/entries)

h_cutflow.Scale(1./nInclusive)
h_trigger_efficiency.Scale(1./entries)
h_FDeff_ctau.Scale(1./entries)
h_FDeff_times_coupling4_ctau.Scale(1./entries)
h_FDeff_coupling2.Scale(1./entries)
h_FDeff_times_coupling4_coupling2.Scale(1./entries)
h_nEvents_Maj_coupling2.Scale(1./nInclusive)
h_nEvents_Dir_coupling2.Scale(1./nInclusive)
h_nEvents_Maj_couplingMuTau.Scale(1./nInclusive)
h_nEvents_Dir_couplingMuTau.Scale(1./nInclusive)
h_BR_Maj_eff_coupling2.Scale(1./nInclusive)
h_BR_Dir_eff_coupling2.Scale(1./nInclusive)

#fit histograms with normalized FD
exp_func = ROOT.TF1("exp_func", "[0]*TMath::Exp([1]*x)", 0.01, 1.5)
exp_func.SetParameters(1, -1) #set initial parameters
exp_func.SetLineColor(1+massPoint)
if (h_FD_HN_norm.Integral() != 0):
  h_FD_HN_norm.Scale(1./h_FD_HN_norm.Integral())
h_FD_HN_norm.Fit("exp_func", "R")
decay_constant = exp_func.GetParameter(1)
print("Normalized FD of HN =", 1/decay_constant)

''' 
#define canvas and draw and save histograms

c_cutflow = TCanvas("c_cutflow","canvas title", 800, 800) 
#h_cutflow.SetMinimum(0)
h_cutflow.Draw("h")
c_cutflow.SaveAs("plots/cutflow.png")

c_FDx_HN = TCanvas("c_FDx_HN","canvas title", 800, 800) 
#c_FDx_HN.SetLogy(1)
h_FDx_HN.Draw()
c_FDx_HN.SaveAs("plots/FDx_HN.png")

c_FDy_HN = TCanvas("c_FDy_HN","canvas title", 800, 800) 
#c_FDy_HN.SetLogy(1)
h_FDy_HN.Draw()
c_FDy_HN.SaveAs("plots/FDy_HN.png")

c_FDt_HN = TCanvas("c_FDt_HN","canvas title", 800, 800) 
c_FDt_HN.SetLogy(1)
h_FDt_HN.Draw()
c_FDt_HN.SaveAs("plots/FDt_HN.png")

c_FDz_HN = TCanvas("c_FDz_HN","canvas title", 800, 800) 
c_FDz_HN.SetLogy(1)
h_FDz_HN.Draw()
c_FDz_HN.SaveAs("plots/FDz_HN.png")

c_FD_HN = TCanvas("c_FD_HN","canvas title", 800, 800) 
c_FD_HN.SetLogy(1)
h_FD_HN.Draw()
c_FD_HN.SaveAs("plots/FD_HN.png")
c_FD_HN.SaveAs("plots/FD_HN.root")

c_FD_HN_norm = TCanvas("c_FD_HN_norm","canvas title", 800, 800) 
c_FD_HN_norm.SetLogy(1)
h_FD_HN_norm.Draw()
exp_func.Draw("same")
c_FD_HN_norm.SaveAs("plots/FD_HN_norm.png")
c_FD_HN_norm.SaveAs("plots/FD_HN_norm.root")

#energy canvas
c_energy_B = TCanvas("c_energy_B","canvas title", 800, 800) 
c_energy_B.SetLogy(1)
h_energy_B.Draw()
c_energy_B.SaveAs("plots/energy_B.png")

c_energy_HN = TCanvas("c_energy_HN","canvas title", 800, 800) 
c_energy_HN.SetLogy(1)
h_energy_HN.Draw()
c_energy_HN.SaveAs("plots/energy_HN.png")

c_energy_tau = TCanvas("c_energy_tau","canvas title", 800, 800) 
c_energy_tau.SetLogy(1)
h_energy_tau.Draw()
c_energy_tau.SaveAs("plots/energy_tau.png")

c_energy_mu1 = TCanvas("c_energy_mu1","canvas title", 800, 800) 
c_energy_mu1.SetLogy(1)
h_energy_mu1.Draw()
c_energy_mu1.SaveAs("plots/energy_mu1.png")

c_energy_mu2 = TCanvas("c_energy_mu2","canvas title", 800, 800) 
c_energy_mu2.SetLogy(1)
h_energy_mu2.Draw()
c_energy_mu2.SaveAs("plots/energy_mu2.png")

c_energy_mu3 = TCanvas("c_energy_mu3","canvas title", 800, 800) 
c_energy_mu3.SetLogy(1)
h_energy_mu3.Draw()
c_energy_mu3.SaveAs("plots/energy_mu3.png")

#pt canvas
c_pt_B = TCanvas("c_pt_B","canvas title", 800, 800) 
c_pt_B.SetLogy(1)
h_pt_B.Draw()
c_pt_B.SaveAs("plots/pt_B.png")

c_pt_HN = TCanvas("c_pt_HN","canvas title", 800, 800) 
c_pt_HN.SetLogy(1)
h_pt_HN.Draw()
c_pt_HN.SaveAs("plots/pt_HN.png")

c_pt_mu1 = TCanvas("c_pt_mu1","canvas title", 800, 800) 
c_pt_mu1.SetLogy(1)
h_pt_mu1.Draw()
c_pt_mu1.SaveAs("plots/pt_mu1.png")

c_pt_mu2 = TCanvas("c_pt_mu2","canvas title", 800, 800) 
c_pt_mu2.SetLogy(1)
h_pt_mu2.Draw()
c_pt_mu2.SaveAs("plots/pt_mu2.png")

c_pt_mu3 = TCanvas("c_pt_mu3","canvas title", 800, 800) 
c_pt_mu3.SetLogy(1)
h_pt_mu3.Draw()
c_pt_mu3.SaveAs("plots/pt_mu3.png")

c_pt_tau = TCanvas("c_pt_tau","canvas title", 800, 800) 
c_pt_tau.SetLogy(1)
h_pt_tau.Draw()
c_pt_tau.SaveAs("plots/pt_tau.png")


c_pt_neutrino = TCanvas("c_pt_neutrino","canvas title", 800, 800) 
c_pt_neutrino.SetLogy(1)
h_pt_neutrino.Draw()
c_pt_neutrino.SaveAs("plots/pt_neutrino.png")



c_phi_mu1 = TCanvas("c_phi_mu1","canvas title", 800, 800) 
h_phi_mu1.Draw()
c_phi_mu1.SaveAs("plots/phi_mu1.png")

c_phi_mu2 = TCanvas("c_phi_mu2","canvas title", 800, 800) 
h_phi_mu2.Draw()
c_phi_mu2.SaveAs("plots/phi_mu2.png")

c_phi_mu3 = TCanvas("c_phi_mu3","canvas title", 800, 800) 
h_phi_mu3.Draw()
c_phi_mu3.SaveAs("plots/phi_mu3.png")

c_phi_tau = TCanvas("c_phi_tau","canvas title", 800, 800) 
h_phi_tau.Draw()
c_phi_tau.SaveAs("plots/phi_tau.png")

c_phi_neutrino = TCanvas("c_phi_neutrino","canvas title", 800, 800) 
h_phi_neutrino.Draw()
c_phi_neutrino.SaveAs("plots/phi_neutrino.png")

c_deltaphi_mu1mu2 = TCanvas("c_deltaphi_mu1mu2","canvas title", 800, 800) 
h_deltaphi_mu1mu2.Draw()
c_deltaphi_mu1mu2.SaveAs("plots/deltaphi_mu1mu2.png")

c_deltaphi_mu1tau = TCanvas("c_deltaphi_mu1tau","canvas title", 800, 800) 
h_deltaphi_mu1tau.Draw()
c_deltaphi_mu1tau.SaveAs("plots/deltaphi_mu1tau.png")

c_deltaphi_mu2tau = TCanvas("c_deltaphi_mu2tau","canvas title", 800, 800) 
h_deltaphi_mu2tau.Draw()
c_deltaphi_mu2tau.SaveAs("plots/deltaphi_mu2tau.png")

c_deltaphi_mu1mu3 = TCanvas("c_deltaphi_mu1mu3","canvas title", 800, 800) 
h_deltaphi_mu1mu3.Draw()
c_deltaphi_mu1mu3.SaveAs("plots/deltaphi_mu1mu3.png")

c_deltaphi_mu2mu3 = TCanvas("c_deltaphi_mu2mu3","canvas title", 800, 800) 
h_deltaphi_mu2mu3.Draw()
c_deltaphi_mu2mu3.SaveAs("plots/deltaphi_mu2mu3.png")

c_deltaR_mu2tau = TCanvas("c_deltaR_mu2tau", "canvas title", 800, 800)
h_deltaR_mu2tau.Draw()
c_deltaR_mu2tau.SaveAs("plots/deltaR_mu2tau.png")

c_deltaR_mu2mu3 = TCanvas("c_deltaR_mu2mu3", "canvas title", 800, 800)
h_deltaR_mu2mu3.Draw()
c_deltaR_mu2mu3.SaveAs("plots/deltaR_mu2mu3.png")

c_dalitz_HN = TCanvas("c_dalitz_HN", "canvas title", 800, 800)
#c_dalitz_HN.SetLogz(1)
h_dalitz_HN.Draw("colz")
c_dalitz_HN.SaveAs("plots/dalitz_HN.png")


c_pt_m1vsHN = TCanvas("c_ptmu1vsHN", "canvastitle", 800, 800)
c_pt_m1vsHN.SetLogz(1)
h_pt_m1vsHN.Draw("colz")
c_pt_m1vsHN.SaveAs("plots/pt_m1vsHN.png")

c_pt_m1vssum = TCanvas("c_ptmu1vssum", "canvastitle", 800, 800)
c_pt_m1vssum.SetLogz(1)
h_pt_m1vssum.Draw("colz")
c_pt_m1vssum.SaveAs("plots/pt_m1vssum.png")


c_pt_Bvsall = TCanvas("c_ptBvsall", "canvastitle", 800, 800)
c_pt_Bvsall.SetLogz(1)
h_pt_Bvsall.Draw("colz")
c_pt_Bvsall.SaveAs("plots/pt_Bvsall.png")


c_m_HN = TCanvas("c_m_HN", "canvastitle", 800, 800)
c_m_HN.SetLogy(1)
h_m_HN.Draw()
c_m_HN.SaveAs("plots/m_HN.png")

c_m_HNprime = TCanvas("c_m_HNprime", "canvastitle", 800, 800)
c_m_HNprime.SetLogy(1)
h_m_HNprime.Draw()
c_m_HNprime.SaveAs("plots/m_HNprime.png")
'''


# Saving TH1 histograms
outFile=TFile("histograms_mass_"+str(M)+".root","RECREATE")
outFile.cd()
outFile.mkdir("HN_mass_"+str(M))
outFile.cd("HN_mass_"+str(M))

h_cutflow.Write()
h_trigger_efficiency.Write()
h_FDeff_ctau.Write()
h_FDeff_times_coupling4_ctau.Write()
h_FDeff_coupling2.Write()
h_FDeff_times_coupling4_coupling2.Write()
h_nEvents_Maj_coupling2.Write()
h_nEvents_Dir_coupling2.Write()
h_nEvents_Maj_couplingMuTau.Write()
h_nEvents_Dir_couplingMuTau.Write()
h_BR_Maj_coupling2.Write()
h_BR_Dir_coupling2.Write()
h_BR_Maj_couplingMuTau.Write()
h_BR_Dir_couplingMuTau.Write()
h_BR_Maj_eff_coupling2.Write()
h_BR_Dir_eff_coupling2.Write()
h_FDx_HN.Write()
h_FDy_HN.Write()
h_FDt_HN.Write()
h_FDz_HN.Write()
h_FD_HN.Write()
h_FD_HN_norm.Write()
h_FD_mu1_mu3.Write()
h_pt_mu1.Write()
h_pt_mu2.Write()
h_pt_mu1mu2.Write()
h_pt_mu1mu3.Write()
h_pt_mu2mu3.Write()
h_pt_mu3.Write()
h_pt_tau.Write()
h_pt_HN.Write()
h_pt_B.Write()
h_pt_neutrino.Write()
h_phi_mu1.Write()
h_phi_mu2.Write()
h_phi_mu3.Write()
h_phi_tau.Write()
h_phi_neutrino.Write()
h_deltaphi_mu1mu2.Write()
h_deltaphi_mu1mu3.Write()
h_deltaphi_mu2mu3.Write()
h_deltaphi_mu1tau.Write()
h_deltaphi_mu2tau.Write()
h_deltaR_mu2tau.Write()
h_deltaeta_mu1mu2.Write()
h_deltaeta_mu1mu3.Write()
h_deltaeta_mu2mu3.Write()
h_deltaR_mu1mu2.Write()
h_deltaR_mu1mu3.Write()
h_deltaR_mu2mu3.Write()
h_theta_HN.Write()
h_true_theta_HN.Write()
h_energy_HN.Write()
h_energy_B.Write()
h_energy_tau.Write()
h_energy_tau_HNRF.Write()
h_energy_mu1.Write()
h_energy_mu2.Write()
h_energy_mu3.Write()
h_energy_mu3_HNRF.Write()
h_m_HN.Write()
h_m_HNprime.Write()
h_m_Jpsi.Write()
h_m_mu1mu2.Write()
h_m_mu1mu3.Write()
h_m_mu2mu3.Write()
h_m_B.Write()
h_dalitz_HN.Write()

outFile.Close()