import ROOT
from ROOT import TFile, gDirectory
from ROOT import TH1D, TH2D, TCanvas, TLorentzVector
import math
from math import sqrt
from tqdm import tqdm
import os
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1,0,1)

tau_energy_thr = 0 #30
muon_pt_thr = 5 #5

cutflow_labels = ["generated","LHCb geo","mu1 trig","mu2 trig", "mu3 trig", "J/#psi mass", "B_{c} mass"]


backID = 2

nInclusive = int(sys.argv[1])


infile = TFile( '/eos/user/m/mmourad/HN_simulation/trees/gamma2_tree.root' )
tree = gDirectory.Get('DecayTree')
entries = tree.GetEntries()

# define Lorentz vectors and variables

mu1 = TLorentzVector()
mu2 = TLorentzVector()
mu3 = TLorentzVector()
neutrino = TLorentzVector()
B = TLorentzVector()
HN = TLorentzVector()

#define histograms

h_cutflow = TH1D("h_cutflow","#mu#mu#mu",len(cutflow_labels),0,len(cutflow_labels))
for i in range(len(cutflow_labels)):
  h_cutflow.GetXaxis().SetBinLabel(i+1,cutflow_labels[i])
h_cutflow.SetBinContent(1,nInclusive)

h_FDx_HN = TH1D("h_FDx_HN"," Flight distance of HN in x; Flight distance of HN in x [mm]",50,0,200)
h_FDy_HN = TH1D("h_FDy_HN"," Flight distance of HN in y; Flight distance of HN in y [mm]",50,0,200)
h_FDt_HN = TH1D("h_FDt_HN","Transverse  Flight distance of HN;Transverse  Flight distance of HN[mm]",50,0,200)
h_FDz_HN = TH1D("h_FDz_HN"," Flight distance of HN in z; Flight distance of HN in z [mm]",50,0,500)
h_FD_HN = TH1D("h_FD_HN"," Flight distance of HN; Flight distance of HN [mm]",100,0,500)
h_FD_HN_norm = TH1D("h_FD_HN_norm"," Flight distance of HN in its RF; Flight distance of HN in its RF [mm]",100,0,4)

h_mu12_vtxZ = TH1D("h_mu12_vtxZ","vertex distance of mu1 and mu2 in Z;vertex distance of mu1 and mu2 in Z [mm]",50,0,10)
h_mu12_vtxT = TH1D("h_mu12_vtxT","transverse vertex distance of mu1 and mu2;transverse vertex distance of mu1 and mu2 [mm]",50,0,10)
h_mu12_vtx = TH1D("h_mu12_vtx","total vertex distance of mu1 and mu2;total vertex distance of mu1 and mu2 [mm]",50,0,10)

h_mu1mu3_vtxZ = TH1D("h_mu1mu3_vtxZ","vertex distance of mu1 and mu3 in Z;vertex distance of mu1 and mu3 in Z [mm]",50,0,10)
h_mu1mu3_vtxT = TH1D("h_mu1mu3_vtxT","transverse vertex distance of mu1 and mu3;transverse vertex distance of mu1 and mu3 [mm]",50,0,10)
h_mu1mu3_vtx = TH1D("h_mu1mu3_vtx","total vertex distance of mu1 and mu3;total vertex distance of mu1 and mu3 [mm]",50,0,10)


h_pt_mu1 = TH1D("h_pt_mu1","p_{T} of #mu_{1};p_{T} [GeV]",50,muon_pt_thr,100)
h_pt_mu2 = TH1D("h_pt_mu2","p_{T} of #mu_{2};p_{T} [GeV]",50,muon_pt_thr,100)
h_pt_mu1mu2 = TH1D("h_pt_mu1mu2","p_{T} of #mu_{1}+#mu_{2};p_{T} [GeV]",50,0,100)
h_pt_mu3 = TH1D("h_pt_mu3","p_{T} of #mu_{3};p_{T} [GeV]",50,0,100)
h_pt_HN = TH1D("h_pt_HN","p_{T} of HN;p_{T} [GeV]",50,0,100)
h_pt_B = TH1D("h_pt_B","p_{T} of B_{c};p_{T} [GeV]",50,0,100)
h_pt_neutrino = TH1D("h_pt_neutrino","p_{T} of #nu_{#mu3};p_{T} [GeV]",50,0,30)

h_phi_mu1 = TH1D("h_phi_mu1","#phi of #mu_{1};#phi",50,-3.1415,3.1415)
h_phi_mu2 = TH1D("h_phi_mu2","#phi of #mu_{2};#phi",50,-3.1415,3.1415)
h_phi_mu3 = TH1D("h_phi_mu3","#phi of #mu_{3};#phi",50,-3.1415,3.1415)
h_phi_neutrino = TH1D("h_phi_neutrino","#phi of #nu;#phi",50,-3.1415,3.1415)

h_deltaphi_mu1mu2 = TH1D("h_deltaphi_mu1mu2","#Delta#phi #mu_{1} and #mu_{2};#Delta#phi",30,0,1)
h_deltaphi_mu1mu3 = TH1D("h_deltaphi_mu1mu3","#Delta#phi #mu_{1} and #mu_{3};#Delta#phi",30,0,1)
h_deltaphi_mu2mu3 = TH1D("h_deltaphi_mu2mu3","#Delta#phi #mu_{2} and #mu_{3};#Delta#phi",30,0,1)


h_deltaR_mu2mu3 = TH1D("h_deltaR_mu2mu3","#Delta R between #mu_{2} and #mu3;#Delta R", 45,0,1.5)
#delta R is square root of delta phi squared + delta eta squared, 3d Lorentz invariant angle

h_energy_HN = TH1D("h_energy_HN","Energy of HN;E_{HN} [GeV]",100,0,1000)
h_energy_B = TH1D("h_energy_B","Energy of B_{c};E_{B} [GeV]",100,0,1000)
h_energy_mu3 = TH1D("h_energy_mu3","Energy of #mu_{3};E_{#mu_{3}} [GeV]",100, 0,1000)
h_energy_mu1 = TH1D("h_energy_mu1","Energy of #mu_{1};E_{#mu_{1}} [GeV]",100, 0,1000)
h_energy_mu2 = TH1D("h_energy_mu2","Energy of #mu_{2};E_{#mu_{2}} [GeV]",100, 0,1000)

#Dalitz plot for HN
h_dalitz_HN = TH2D(" h_dalitz_HN ","Dalitzplot for B_{c}+ -> #mu+#mu+#mu3-;m^{2}_{#mu_{1}#mu3} [GeV];m^{2}_{#mu_{2}#mu_{3}} [GeV]" ,100 ,0 ,50, 100, 0 , 10)
#h_m_HN = TH1D("h_m_HN", "Reconstructed Mass of HN(without neutrino) ; M [GeV]", 150, 0, 6)
h_m_HN = TH1D("h_m_HN", "Reconstructed Mass of HN(without neutrino) ; M [GeV]", 150, 0, 5)
h_m_HNprime = TH1D("h_m_HNprime", "Generated Mass of HN; M [GeV]", 150, 0 , 6)

h_m_Jpsi = TH1D("h_m_Jpsi", "Invariant mass of mu1 and mu2; invariant mass of mu1 and mu2 [GeV]", 100, 2.8 , 3.4)

h_m_B = TH1D("h_m_B", "Visible mass of B_{c};visible mass of B [GeV]", 100, 0 , 7)

#2d plots of pt
h_pt_m1vsHN = TH2D("h_pt_m1vsHN","p_{T} of HN vs p_{T} of #mu_{1};p_{T} of HN [GeV]; p_{T} of #mu_{1} [GeV]", 40, 0, 100, 40, 0, 100)
h_pt_m1vssum = TH2D("h_pt_m1vssum","Sum of p_{T} of #mu_{2}, #mu_{3}, #bar#nu_{#mu} vs p_{T} of #mu_{1};p_{T} of sum [GeV]; p_{T} of #mu_{1} [GeV]", 40, 0, 100, 40, 0, 100)
h_pt_Bvsall = TH2D("h_pt_Bvsall","Sum of p_{T} of #mu_{1}, #mu_{2}, #mu_{3}, #bar#nu_{#mu} vs p_{T} of B_{c};p_{T} of sum [GeV]; p_{T} of B_{c} [GeV]", 40, 0, 100, 40, 0, 100)

#label axes of histograms
#either with semi-column method like above or:
#example: myHist.GetXaxis().SetTitle("theta")
#myHist.GetYaxis().SetTitle("K_{pT} [GeV]")

#loop through events
pbar = tqdm(total=entries, unit="")
accepted = 0

for event in range( entries ):
 pbar.update()
 tree.GetEntry(event)
 #cuts
 h_cutflow.Fill(1)
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
 #HN.SetPtEtaPhiE(tree.HN_PT,tree.HN_eta, tree.HN_phi, tree.HN_E)
 B.SetPtEtaPhiE(tree.B_PT,tree.B_eta, tree.B_phi, tree.B_E)
 neutrino.SetPtEtaPhiE(tree.neutrino_PT,tree.neutrino_eta, tree.neutrino_phi, tree.neutrino_E)
 #neutrino = B - mu1 - mu2 - mu3  
 #fill histograms
 h_pt_B.Fill(B.Pt())
 h_energy_B.Fill(B.E())
 h_pt_HN.Fill(-1)
 h_energy_HN.Fill(-1)
 h_energy_mu3.Fill(mu3.E())
 h_energy_mu1.Fill(mu1.E())
 h_energy_mu2.Fill(mu2.E())
 h_pt_mu1.Fill(mu1.Pt())
 h_pt_mu2.Fill(mu2.Pt())
 h_pt_mu1mu2.Fill((mu1+mu2).Pt())
 h_pt_mu3.Fill(mu3.Pt())
 h_pt_neutrino.Fill(neutrino.Pt())
 h_phi_mu1.Fill(mu1.Phi())
 h_phi_mu2.Fill(mu2.Phi())
 h_phi_mu3.Fill(mu3.Phi())
 h_phi_neutrino.Fill(neutrino.Phi())
 h_deltaphi_mu1mu2.Fill(abs(mu1.DeltaPhi(mu2)))
 h_deltaphi_mu1mu3.Fill(abs(mu1.DeltaPhi(mu3)))
 h_deltaphi_mu2mu3.Fill(abs(mu2.DeltaPhi(mu3)))
 h_deltaR_mu2mu3.Fill(abs(mu2.DeltaR(mu3))) 
 #calculate flight distance of heavy neutrino from decay vertices of B and HN
 #HN_FDx = tree.HN_vtxX-tree.B_vtxX
 #HN_FDy = tree.HN_vtxY-tree.B_vtxY
 #HN_FDz = tree.HN_vtxZ-tree.B_vtxZ
 HN_FDx = -1
 HN_FDy = -1
 HN_FDz = -1
 HN_FDt = -1
 HN_FD = -1
 h_FDx_HN.Fill(HN_FDx)
 h_FDy_HN.Fill(HN_FDy)
 h_FDt_HN.Fill(HN_FDt)
 h_FDz_HN.Fill(HN_FDz)
 h_FD_HN.Fill(HN_FD)
 
 mu12_vtxX = tree.mu1_origX - tree.mu2_origX
 mu12_vtxY = tree.mu1_origY - tree.mu2_origY
 mu12_vtxZ = tree.mu1_origZ - tree.mu2_origZ
 mu12_vtxT = sqrt(mu12_vtxX**2 + mu12_vtxY**2)
 mu12_vtx = sqrt(mu12_vtxX**2 + mu12_vtxY**2 + mu12_vtxZ**2)
 
 h_mu12_vtxZ.Fill(mu12_vtxZ)
 h_mu12_vtxT.Fill(mu12_vtxT)
 h_mu12_vtx.Fill(mu12_vtx)
 
 mu1mu3_vtxX = tree.mu1_origX - tree.mu2_origX
 mu1mu3_vtxY = tree.mu1_origY - tree.mu2_origY
 mu1mu3_vtxZ = tree.mu1_origZ - tree.mu2_origZ
 mu1mu3_vtxT = sqrt(mu1mu3_vtxX**2 + mu1mu3_vtxY**2)
 mu1mu3_vtx = sqrt(mu1mu3_vtxX**2 + mu1mu3_vtxY**2 + mu1mu3_vtxZ**2)
 
 h_mu1mu3_vtxZ.Fill(mu1mu3_vtxZ)
 h_mu1mu3_vtxT.Fill(mu1mu3_vtxT)
 h_mu1mu3_vtx.Fill(mu1mu3_vtx)
 
 #caculate boost factor of HN 
 #gamma_HN = HN.E() / HN.M()
 gamma_HN = 1
 #fill histogram with FD of HN in its rest frame
 h_FD_HN_norm.Fill(HN_FD/gamma_HN)
 #invariant mass
 m_Jpsi = (mu1+mu2).M()
 h_m_Jpsi.Fill(m_Jpsi)
 msquared_mu1mu3 = (mu1 + mu3).M() **2
 msquared_mu2mu3 = (mu2 + mu3).M() **2
 h_dalitz_HN.Fill(msquared_mu1mu3, msquared_mu2mu3)
 m_HN = (mu2 + mu3).M()
 h_m_HN.Fill(m_HN)
 m_HNprime = (mu2 + mu3 + neutrino).M()
 h_m_HNprime.Fill(m_HNprime)
 h_m_B.Fill((mu1+mu2+mu3).M())
 #2d plots
 #h_pt_m1vsHN.Fill(HN.Pt(), mu1.Pt())
 h_pt_m1vssum.Fill((mu2 + mu3 + neutrino).Pt(), mu1.Pt())
 h_pt_Bvsall.Fill((mu2 + mu1 + mu3 + neutrino).Pt(), B.Pt())
 
 
 if ((mu1+mu2).M() > 3.0 and (mu1+mu2).M() < 3.2):
   continue
 h_cutflow.Fill(5)
 if ((mu1+mu2+mu3).M() > 5.8):
   continue
 h_cutflow.Fill(6)
 
pbar.close()

print("geo_accepted: ", entries)
print("analysis_accepted: ", accepted)
print("analysis selection efficiency: ", accepted/entries)

h_cutflow.Scale(1./nInclusive)


#fit histograms with normalized FD
exp_func = ROOT.TF1("exp_func", "[0]*TMath::Exp([1]*x)", 0.01, 15)
exp_func.SetParameters(100000, -0.06) #set initial parameters
exp_func.SetLineColor(10+backID)
h_FD_HN_norm.Fit("exp_func", "R")
decay_constant = exp_func.GetParameter(1)
print("Normalized FD of HN =", 1/decay_constant)

 
#define canvas and draw and save histograms

c_cutflow = TCanvas("c_cutflow","canvas title", 800, 800) 
#h_cutflow.SetMinimum(0)
c_cutflow.SetLogy(1)
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



c_mu12_vtxZ = TCanvas("c_mu12_vtxZ","canvas title", 800, 800) 
c_mu12_vtxZ.SetLogy(1)
h_mu12_vtxZ.Draw()
c_mu12_vtxZ.SaveAs("plots/mu12_vtxZ.png")

c_mu12_vtxT = TCanvas("c_mu12_vtxT","canvas title", 800, 800) 
c_mu12_vtxT.SetLogy(1)
h_mu12_vtxT.Draw()
c_mu12_vtxT.SaveAs("plots/mu12_vtxT.png")

c_mu12_vtx = TCanvas("c_mu12_vtx","canvas title", 800, 800) 
c_mu12_vtx.SetLogy(1)
h_mu12_vtx.Draw()
c_mu12_vtx.SaveAs("plots/mu12_vtx.png")

c_mu1mu3_vtxZ = TCanvas("c_mu1mu3_vtxZ","canvas title", 800, 800) 
c_mu1mu3_vtxZ.SetLogy(1)
h_mu1mu3_vtxZ.Draw()
c_mu1mu3_vtxZ.SaveAs("plots/mu1mu3_vtxZ.png")

c_mu1mu3_vtxT = TCanvas("c_mu1mu3_vtxT","canvas title", 800, 800) 
c_mu1mu3_vtxT.SetLogy(1)
h_mu1mu3_vtxT.Draw()
c_mu1mu3_vtxT.SaveAs("plots/mu1mu3_vtxT.png")

c_mu1mu3_vtx = TCanvas("c_mu1mu3_vtx","canvas title", 800, 800) 
c_mu1mu3_vtx.SetLogy(1)
h_mu1mu3_vtx.Draw()
c_mu1mu3_vtx.SaveAs("plots/mu1mu3_vtx.png")

#energy canvas
c_energy_B = TCanvas("c_energy_B","canvas title", 800, 800) 
c_energy_B.SetLogy(1)
h_energy_B.Draw()
c_energy_B.SaveAs("plots/energy_B.png")

c_energy_HN = TCanvas("c_energy_HN","canvas title", 800, 800) 
c_energy_HN.SetLogy(1)
h_energy_HN.Draw()
c_energy_HN.SaveAs("plots/energy_HN.png")

c_energy_mu3 = TCanvas("c_energy_mu3","canvas title", 800, 800) 
c_energy_mu3.SetLogy(1)
h_energy_mu3.Draw()
c_energy_mu3.SaveAs("plots/energy_mu3.png")

c_energy_mu1 = TCanvas("c_energy_mu1","canvas title", 800, 800) 
c_energy_mu1.SetLogy(1)
h_energy_mu1.Draw()
c_energy_mu1.SaveAs("plots/energy_mu1.png")

c_energy_mu2 = TCanvas("c_energy_mu2","canvas title", 800, 800) 
c_energy_mu2.SetLogy(1)
h_energy_mu2.Draw()
c_energy_mu2.SaveAs("plots/energy_mu2.png")


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

c_phi_neutrino = TCanvas("c_phi_neutrino","canvas title", 800, 800) 
h_phi_neutrino.Draw()
c_phi_neutrino.SaveAs("plots/phi_neutrino.png")

c_deltaphi_mu1mu2 = TCanvas("c_deltaphi_mu1mu2","canvas title", 800, 800) 
h_deltaphi_mu1mu2.Draw()
c_deltaphi_mu1mu2.SaveAs("plots/deltaphi_mu1mu2.png")

c_deltaphi_mu1mu3 = TCanvas("c_deltaphi_mu1mu3","canvas title", 800, 800) 
h_deltaphi_mu1mu3.Draw()
c_deltaphi_mu1mu3.SaveAs("plots/deltaphi_mu1mu3.png")

c_deltaphi_mu2mu3 = TCanvas("c_deltaphi_mu2mu3","canvas title", 800, 800) 
h_deltaphi_mu2mu3.Draw()
c_deltaphi_mu2mu3.SaveAs("plots/deltaphi_mu2mu3.png")

c_deltaR_mu2mu3 = TCanvas("c_deltaR_mu2mu3", "canvas title", 800, 800)
h_deltaR_mu2mu3.Draw()
c_deltaR_mu2mu3.SaveAs("plots/deltaR_mu2mu3.png")

c_dalitz_HN = TCanvas("c_dalitz_HN", "canvas title", 800, 800)
c_dalitz_HN.SetLogz(1)
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

c_m_Jpsi = TCanvas("c_m_Jpsi", "canvastitle", 800, 800)
h_m_Jpsi.Draw()
c_m_Jpsi.SaveAs("plots/m_Jpsi.png")



# Saving TH1 histograms
outFile=TFile("histograms.root","RECREATE")
outFile.cd()
outFile.mkdir("gamma"+str(backID))
outFile.cd("gamma"+str(backID))

h_cutflow.Write()
h_FDx_HN.Write()
h_FDy_HN.Write()
h_FDt_HN.Write()
h_FDz_HN.Write()
h_FD_HN.Write()
h_FD_HN_norm.Write()
h_pt_mu1.Write()
h_pt_mu2.Write()
h_pt_mu1mu2.Write()
h_pt_mu3.Write()
h_pt_HN.Write()
h_pt_B.Write()
h_pt_neutrino.Write()
h_phi_mu1.Write()
h_phi_mu2.Write()
h_phi_mu3.Write()
h_phi_neutrino.Write()
h_deltaphi_mu1mu2.Write()
h_deltaphi_mu1mu3.Write()
h_deltaphi_mu2mu3.Write()
h_deltaR_mu2mu3.Write()
h_energy_HN.Write()
h_energy_B.Write()
h_energy_mu3.Write()
h_energy_mu1.Write()
h_energy_mu2.Write()
h_m_HN.Write()
h_m_HNprime.Write()
h_m_Jpsi.Write()
h_m_B.Write()

outFile.Close()