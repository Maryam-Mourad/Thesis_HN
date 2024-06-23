import ROOT
from ROOT import TFile, gDirectory, TH1D, TH2D, TCanvas, TLegend
import os
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1,0,1)


mass_list=sys.argv[1:]


merged_file = TFile("histograms.root","READ")


legend = TLegend(0.7, 0.75, 0.88, 0.88) # x0, y0, x1, y1
legend.SetFillStyle(0)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.SetShadowColor(0)
legend.SetTextSize(0.035)


c_cutflow = TCanvas("c_cutflow","canvas title", 800, 800) 
c_cutflow.SetLogy(1)
for massID in range(len(mass_list)):
    h_cutflow=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_cutflow")
    h_cutflow.SetLineColor(1+massID)
    h_cutflow.SetLineWidth(3)
    legend.AddEntry(h_cutflow," "+mass_list[massID]+" GeV","l")
    h_cutflow.Draw("hesame")
legend.Draw()
c_cutflow.SaveAs("plots/cutflow.png")

c_FD_HN = TCanvas("c_FD_HN","canvas title", 800, 800) 
c_FD_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_FD_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FD_HN")
    h_FD_HN.SetLineColor(1+massID)
    h_FD_HN.SetLineWidth(3)
    h_FD_HN.Draw("hesame")
legend.Draw()
c_FD_HN.SaveAs("plots/FD_HN.png")




c_FDx_HN = TCanvas("c_FDx_HN","canvas title", 800, 800) 
#c_FDx_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_FDx_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FDx_HN")
    h_FDx_HN.SetLineColor(1+massID)
    h_FDx_HN.SetLineWidth(3)
    h_FDx_HN.Draw("hesame")
legend.Draw()
c_FDx_HN.SaveAs("plots/FDx_HN.png")

c_FDy_HN = TCanvas("c_FDy_HN","canvas title", 800, 800) 
#c_FDy_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_FDy_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FDy_HN")
    h_FDy_HN.SetLineColor(1+massID)
    h_FDy_HN.SetLineWidth(3)
    h_FDy_HN.Draw("hesame")
legend.Draw()
c_FDy_HN.SaveAs("plots/FDy_HN.png")

c_FDt_HN = TCanvas("c_FDt_HN","canvas title", 800, 800) 
c_FDt_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_FDt_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FDt_HN")
    h_FDt_HN.SetLineColor(1+massID)
    h_FDt_HN.SetLineWidth(3)
    h_FDt_HN.Draw("hesame")
legend.Draw()
c_FDt_HN.SaveAs("plots/FDt_HN.png")

c_FDz_HN = TCanvas("c_FDz_HN","canvas title", 800, 800) 
c_FDz_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_FDz_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FDz_HN")
    h_FDz_HN.SetLineColor(1+massID)
    h_FDz_HN.SetLineWidth(3)
    h_FDz_HN.Draw("hesame")
legend.Draw()
c_FDz_HN.SaveAs("plots/FDz_HN.png")

c_FD_HN_norm = TCanvas("c_FD_HN_norm","canvas title", 800, 800) 
c_FD_HN_norm.SetLogy(1)
minY=-1
maxY=-1
for massID in range(len(mass_list)):
    h_FD_HN_norm=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FD_HN_norm")
    if(h_FD_HN_norm.GetMaximum()>maxY or maxY==-1):
      maxY = h_FD_HN_norm.GetMaximum()
    if(h_FD_HN_norm.GetMinimum()<minY or minY==-1):
      minY = h_FD_HN_norm.GetMinimum()
      
if(minY<=0):
  minY=0.1

for massID in range(len(mass_list)):
    h_FD_HN_norm=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_FD_HN_norm")
    h_FD_HN_norm.SetLineColor(1+massID)
    h_FD_HN_norm.SetLineWidth(3)
    if (massID == 0):
      h_FD_HN_norm.GetYaxis().SetRangeUser(minY,maxY)
    h_FD_HN_norm.Draw("hesame")
legend.Draw()
c_FD_HN_norm.SaveAs("plots/FD_HN_norm.png")

#energy canvas
c_energy_B = TCanvas("c_energy_B","canvas title", 800, 800) 
c_energy_B.SetLogy(1)
for massID in range(len(mass_list)):
    h_energy_B=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_energy_B")
    h_energy_B.SetLineColor(1+massID)
    h_energy_B.SetLineWidth(3)
    h_energy_B.Draw("hesame")
legend.Draw()
c_energy_B.SaveAs("plots/energy_B.png")

c_energy_HN = TCanvas("c_energy_HN","canvas title", 800, 800) 
c_energy_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_energy_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_energy_HN")
    h_energy_HN.SetLineColor(1+massID)
    h_energy_HN.SetLineWidth(3)
    h_energy_HN.Draw("hesame")
legend.Draw()
c_energy_HN.SaveAs("plots/energy_HN.png")

c_energy_pi = TCanvas("c_energy_pi","canvas title", 800, 800) 
c_energy_pi.SetLogy(1)
minY=-1
maxY=-1
for massID in range(len(mass_list)):
    h_energy_pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_energy_pi")
    if(h_energy_pi.GetMaximum()>maxY or maxY==-1):
      maxY = h_energy_pi.GetMaximum()
    if(h_energy_pi.GetMinimum()<minY or minY==-1):
      minY = h_energy_pi.GetMinimum()
      
if(minY<=0):
  minY=0.1

for massID in range(len(mass_list)):
    h_energy_pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_energy_pi")
    h_energy_pi.SetLineColor(1+massID)
    h_energy_pi.SetLineWidth(3)
    if (massID == 0):
      h_energy_pi.GetYaxis().SetRangeUser(minY,maxY)
    h_energy_pi.Draw("hesame")
legend.Draw()
c_energy_pi.SaveAs("plots/energy_pi.png")

c_energy_mu1 = TCanvas("c_energy_mu1","canvas title", 800, 800) 
c_energy_mu1.SetLogy(1)
for massID in range(len(mass_list)):
    h_energy_mu1=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_energy_mu1")
    h_energy_mu1.SetLineColor(1+massID)
    h_energy_mu1.SetLineWidth(3)
    h_energy_mu1.Draw("hesame")
legend.Draw()
c_energy_mu1.SaveAs("plots/energy_mu1.png")

c_energy_mu2 = TCanvas("c_energy_mu2","canvas title", 800, 800) 
c_energy_mu2.SetLogy(1)
for massID in range(len(mass_list)):
    h_energy_mu2=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_energy_mu2")
    h_energy_mu2.SetLineColor(1+massID)
    h_energy_mu2.SetLineWidth(3)
    h_energy_mu2.Draw("hesame")
legend.Draw()
c_energy_mu2.SaveAs("plots/energy_mu2.png")


#pt canvas
c_pt_B = TCanvas("c_pt_B","canvas title", 800, 800) 
c_pt_B.SetLogy(1)
for massID in range(len(mass_list)):
    h_pt_B=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_B")
    h_pt_B.SetLineColor(1+massID)
    h_pt_B.SetLineWidth(3)
    h_pt_B.Draw("hesame")
legend.Draw()
c_pt_B.SaveAs("plots/pt_B.png")

c_pt_HN = TCanvas("c_pt_HN","canvas title", 800, 800) 
c_pt_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_pt_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_HN")
    h_pt_HN.SetLineColor(1+massID)
    h_pt_HN.SetLineWidth(3)
    h_pt_HN.Draw("hesame")
legend.Draw()
c_pt_HN.SaveAs("plots/pt_HN.png")

c_pt_mu1 = TCanvas("c_pt_mu1","canvas title", 800, 800) 
c_pt_mu1.SetLogy(1)
minY=-1
maxY=-1
for massID in range(len(mass_list)):
    h_pt_mu1=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_mu1")
    if(h_pt_mu1.GetMaximum()>maxY or maxY==-1):
      maxY = h_pt_mu1.GetMaximum()
    if(h_pt_mu1.GetMinimum()<minY or minY==-1):
      minY = h_pt_mu1.GetMinimum()
      
if(minY<=0):
  minY=0.1

for massID in range(len(mass_list)):
    h_pt_mu1=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_mu1")
    h_pt_mu1.SetLineColor(1+massID)
    h_pt_mu1.SetLineWidth(3)
    if (massID == 0):
      h_pt_mu1.GetYaxis().SetRangeUser(minY,maxY)
    h_pt_mu1.Draw("hesame")
legend.Draw()
c_pt_mu1.SaveAs("plots/pt_mu1.png")

c_pt_mu2 = TCanvas("c_pt_mu2","canvas title", 800, 800) 
c_pt_mu2.SetLogy(1)
minY=-1
maxY=-1
for massID in range(len(mass_list)):
    h_pt_mu2=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_mu2")
    if(h_pt_mu2.GetMaximum()>maxY or maxY==-1):
      maxY = h_pt_mu2.GetMaximum()
    if(h_pt_mu2.GetMinimum()<minY or minY==-1):
      minY = h_pt_mu2.GetMinimum()
      
if(minY<=0):
  minY=0.1

for massID in range(len(mass_list)):
    h_pt_mu2=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_mu2")
    h_pt_mu2.SetLineColor(1+massID)
    h_pt_mu2.SetLineWidth(3)
    if (massID == 0):
      h_pt_mu2.GetYaxis().SetRangeUser(minY,maxY)
    h_pt_mu2.Draw("hesame")
legend.Draw()
c_pt_mu2.SaveAs("plots/pt_mu2.png")

c_pt_pi = TCanvas("c_pt_pi","canvas title", 800, 800) 
c_pt_pi.SetLogy(1)
for massID in range(len(mass_list)):
    h_pt_pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_pi")
    h_pt_pi.SetLineColor(1+massID)
    h_pt_pi.SetLineWidth(3)
    h_pt_pi.Draw("hesame")
legend.Draw()
c_pt_pi.SaveAs("plots/pt_pi.png")


c_pt_neutrino = TCanvas("c_pt_neutrino","canvas title", 800, 800) 
c_pt_neutrino.SetLogy(1)
for massID in range(len(mass_list)):
    h_pt_neutrino=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_pt_neutrino")
    h_pt_neutrino.SetLineColor(1+massID)
    h_pt_neutrino.SetLineWidth(3)
    h_pt_neutrino.Draw("hesame")
legend.Draw()
c_pt_neutrino.SaveAs("plots/pt_neutrino.png")



c_phi_mu1 = TCanvas("c_phi_mu1","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_mu1=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_phi_mu1")
    h_phi_mu1.SetLineColor(1+massID)
    h_phi_mu1.SetLineWidth(3)
    h_phi_mu1.Draw("hesame")
legend.Draw()
c_phi_mu1.SaveAs("plots/phi_mu1.png")

c_phi_mu2 = TCanvas("c_phi_mu2","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_mu2=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_phi_mu2")
    h_phi_mu2.SetLineColor(1+massID)
    h_phi_mu2.SetLineWidth(3)
    h_phi_mu2.Draw("hesame")
legend.Draw()
c_phi_mu2.SaveAs("plots/phi_mu2.png")

c_phi_pi = TCanvas("c_phi_pi","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_phi_pi")
    h_phi_pi.SetLineColor(1+massID)
    h_phi_pi.SetLineWidth(3)
    h_phi_pi.Draw("hesame")
legend.Draw()
c_phi_pi.SaveAs("plots/phi_pi.png")

c_phi_neutrino = TCanvas("c_phi_neutrino","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_neutrino=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_phi_neutrino")
    h_phi_neutrino.SetLineColor(1+massID)
    h_phi_neutrino.SetLineWidth(3)
    h_phi_neutrino.Draw("hesame")
legend.Draw()
c_phi_neutrino.SaveAs("plots/phi_neutrino.png")

c_deltaphi_mu1mu2 = TCanvas("c_deltaphi_mu1mu2","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_mu1mu2=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_deltaphi_mu1mu2")
    h_phi_mu1mu2.SetLineColor(1+massID)
    h_phi_mu1mu2.SetLineWidth(3)
    h_phi_mu1mu2.Draw("hesame")
legend.Draw()
c_deltaphi_mu1mu2.SaveAs("plots/deltaphi_mu1mu2.png")

c_deltaphi_mu1pi = TCanvas("c_deltaphi_mu1pi","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_mu1pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_deltaphi_mu1pi")
    h_phi_mu1pi.SetLineColor(1+massID)
    h_phi_mu1pi.SetLineWidth(3)
    h_phi_mu1pi.Draw("hesame")
legend.Draw()
c_deltaphi_mu1pi.SaveAs("plots/deltaphi_mu1pi.png")

c_deltaphi_mu2pi = TCanvas("c_deltaphi_mu2pi","canvas title", 800, 800) 
for massID in range(len(mass_list)):
    h_phi_mu2pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_deltaphi_mu2pi")
    h_phi_mu2pi.SetLineColor(1+massID)
    h_phi_mu2pi.SetLineWidth(3)
    h_phi_mu2pi.Draw("hesame")
legend.Draw()
c_deltaphi_mu2pi.SaveAs("plots/deltaphi_mu2pi.png")

c_deltaR_mu2pi = TCanvas("c_deltaR_mu2pi", "canvas title", 800, 800)
for massID in range(len(mass_list)):
    h_R_mu2pi=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_deltaR_mu2pi")
    h_R_mu2pi.SetLineColor(1+massID)
    h_R_mu2pi.SetLineWidth(3)
    h_R_mu2pi.Draw("hesame")
legend.Draw()
c_deltaR_mu2pi.SaveAs("plots/deltaR_mu2pi.png")


c_m_HN = TCanvas("c_m_HN", "canvastitle", 800, 800)
c_m_HN.SetLogy(1)
for massID in range(len(mass_list)):
    h_m_HN=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_m_HN")
    h_m_HN.SetLineColor(1+massID)
    h_m_HN.SetLineWidth(3)
    h_m_HN.Draw("hesame")
legend.Draw()
c_m_HN.SaveAs("plots/m_HN.png")

c_m_HNprime = TCanvas("c_m_HNprime", "canvastitle", 800, 800)
c_m_HNprime.SetLogy(1)
for massID in range(len(mass_list)):
    h_m_HNprime=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_m_HNprime")
    h_m_HNprime.SetLineColor(1+massID)
    h_m_HNprime.SetLineWidth(3)
    h_m_HNprime.Draw("hesame")
legend.Draw()
c_m_HNprime.SaveAs("plots/m_HNprime.png")