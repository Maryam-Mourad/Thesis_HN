import ROOT
from ROOT import TFile, gDirectory, TH1D, TH2D, TCanvas, TLegend
import os
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1,0,1)


mass_list=sys.argv[1:]
backID = [2]
background_processes = ["B_{c}#rightarrowJ/#psi+#mu+#nu_{#mu}"]

hists = {'cutflow':1,'FD_HN':1,'FDx_HN':0,'FDy_HN':0,'FDt_HN':1,'FDz_HN':1,'FD_HN_norm':1,
        'energy_B':1,'energy_HN':1,'energy_mu3':1,'energy_mu1':1,'energy_mu2':1,
        'pt_B':1,'pt_HN':1,'pt_mu1':1,'pt_mu2':1,'pt_mu3':1,'pt_neutrino':1,'pt_mu1mu2':1,
        'phi_mu1':0,'phi_mu2':0,'phi_mu3':0,'phi_neutrino':0,
        'deltaphi_mu1mu2':0,'deltaphi_mu1mu3':0,'deltaphi_mu2mu3':0,'deltaR_mu2mu3':0,
        'm_HN':1,'m_HNprime':1,'m_Jpsi':1,'m_B':1}


merged_file = TFile("signal_muon/histograms_muon.root","READ")


keys = list(hists.keys())
values = list(hists.values())

legend = TLegend(0.65, 0.7, 0.88, 0.88) # x0, y0, x1, y1
legend.SetFillStyle(0)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.SetShadowColor(0)
legend.SetTextSize(0.035)




for h in range(len(hists)):
    hName = keys[h]
    c_temp = TCanvas("name"+str(h),"canvas title", 800, 800) 
    c_temp.SetLogy(values[h])
    minY=-1
    maxY=-1
    for massID in range(len(mass_list)):
        temp=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_"+hName)
        if(temp.GetMaximum()>maxY or maxY==-1):
          maxY = temp.GetMaximum()
        if(temp.GetMinimum()<minY or minY==-1):
          minY = temp.GetMinimum()
    for b in backID:
        temp=merged_file.Get("gamma"+str(b)+"/h_"+hName)
        if(temp.GetMaximum()>maxY or maxY==-1):
          maxY = temp.GetMaximum()
        if(temp.GetMinimum()<minY or minY==-1):
          minY = temp.GetMinimum()
    
    if(minY<=0):
      minY=0.0001
    
    if (values[h]):
      maxY *= 10
    else:
      maxY *= 1.2
    
    for massID in range(len(mass_list)):
        temp=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_"+hName)
        temp.SetLineColor(1+massID)
        temp.SetLineWidth(3)
        if (h == 0):
            legend.AddEntry(temp,"HN "+mass_list[massID]+" GeV","l")
        if (massID == 0):
          temp.GetYaxis().SetRangeUser(minY,maxY)
        temp.Draw("hesame")
    for b in backID:
        temp=merged_file.Get("gamma"+str(b)+"/h_"+hName)
        temp.SetLineColor(10+b)
        temp.SetLineWidth(3)
        temp.SetLineStyle(2)
        if (h == 0):
            #legend.AddEntry(temp,"gamma"+str(b),"l")
            legend.AddEntry(temp,background_processes[0],"l")
        temp.Draw("hesame")
    legend.Draw()
    c_temp.SaveAs("plots_muon/"+hName+".png")

