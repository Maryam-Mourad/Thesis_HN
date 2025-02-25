import ROOT
from ROOT import TFile, gDirectory, TH1D, TH2D, TCanvas, TLegend, TPad
import os
import sys
import numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPalette(1,0,1)
ROOT.gStyle.SetPalette(57) #kBird
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)


mass_list=sys.argv[1:]
backID = [2,42,9]
background_processes = {2:"B_{c}#rightarrowJ/#psi+#mu+#nu_{#mu}",42:"B^{+}#rightarrow#mu#mu#mu+#nu_{#mu}",9:"B^{+}#rightarrowD*+#mu+#nu_{#mu}"}
background_colors = {2:2, 42:3, 9:1}

# meaning of digits: = exclude error bars - not normalize to integral of 1 - only for signal - logarithmic y axis
hists = {'cutflow':101,'trigger_efficiency':101,'FD_HN':11,'FDx_HN':11,'FDy_HN':11,'FDt_HN':11,'FDz_HN':11,'FD_HN_norm':11,'FD_mu1_mu3':1,
        'FDeff_ctau':1110,'FDeff_times_coupling4_ctau':1111,'FDeff_coupling2':1110,'FDeff_times_coupling4_coupling2':1111,
        'BR_Maj_coupling2':1111,'BR_Dir_coupling2':1111,'BR_Maj_eff_coupling2':1111,'BR_Dir_eff_coupling2':1111,
        'nEvents_Maj_coupling2':1111,'nEvents_Dir_coupling2':1111,
        'energy_B':10,'energy_HN':10,'energy_mu3':0,'energy_mu1':0,'energy_mu2':0, 'energy_mu3_HNRF':0, 'energy_tau_HNRF':10,
        'pt_B':10,'pt_HN':10,'pt_mu1':0,'pt_mu2':0,'pt_mu3':0,'pt_neutrino':10,'pt_mu1mu2':0,'pt_mu1mu3':0,'pt_mu2mu3':0,
        'phi_mu1':0,'phi_mu2':0,'phi_mu3':0,'phi_neutrino':0,
        'deltaphi_mu1mu2':0,'deltaphi_mu1mu3':0,'deltaphi_mu2mu3':0,
        'deltaeta_mu1mu2':0,'deltaeta_mu1mu3':0,'deltaeta_mu2mu3':0,
        'deltaR_mu1mu2':0,'deltaR_mu1mu3':0,'deltaR_mu2mu3':0, 'theta_HN':0, 'true_theta_HN':10,
        'm_HN':0,'m_HNprime':10,'m_B':0,'m_mu1mu2':0,'m_mu1mu3':0,'m_mu2mu3':0}


merged_file = TFile("signal_muon/histograms_muon.root","READ")


keys = list(hists.keys())
values = list(hists.values())

legend = TLegend(0.63, 0.5, 0.83, 0.87) # x0, y0, x1, y1
legend.SetFillStyle(0)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.SetShadowColor(0)
legend.SetTextSize(0.035)

noBackLegend = TLegend(0.63, 0.72, 0.83, 0.87) # x0, y0, x1, y1
noBackLegend.SetFillStyle(0)
noBackLegend.SetFillColor(0)
noBackLegend.SetLineColor(0)
noBackLegend.SetShadowColor(0)
noBackLegend.SetTextSize(0.035)




for h in range(len(hists)):
    hName = keys[h]
    log = values[h]%10
    onlySig = int(values[h]/10)%10
    normalize = int(values[h]/100)%10
    noError = int(values[h]/1000)%10
    c_temp = TCanvas("name"+str(h),"canvas title", 800, 800) 
    #pad_temp = TPad("pad_temp", "pad", 0.04, 0.04, 0.99, 0.99)
    #pad_temp.cd()
    c_temp.SetLogy(log)
    
    
    c_temp.SetTopMargin(0.09)
    c_temp.SetBottomMargin(0.1)
    c_temp.SetLeftMargin(0.13)
    c_temp.SetRightMargin(0.025)
    
    minY=-1
    maxY=-1
    for massID in range(len(mass_list)):
        if(float(mass_list[massID]) == 4.5):
          continue
        temp=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_"+hName)
        if not normalize and temp.Integral()!=0:
          temp.Scale(1./temp.Integral())
        if(temp.GetMaximum()>maxY or maxY==-1):
          maxY = temp.GetMaximum()
        if(temp.GetMinimum()<minY or minY==-1):
          minY = temp.GetMinimum()
    for b in backID:
        if (onlySig):
          continue
        temp=merged_file.Get("gamma"+str(b)+"/h_"+hName)
        if not normalize and temp.Integral()!=0:
          temp.Scale(1./temp.Integral())
        if(temp.GetMaximum()>maxY or maxY==-1):
          maxY = temp.GetMaximum()
        if(temp.GetMinimum()<minY or minY==-1):
          minY = temp.GetMinimum()
    
    if(minY<=0):
      minY=1e-18
    
    #if(hName == 'nEvents_Maj_coupling2' or hName == 'nEvents_Dir_coupling2'):
      #minY=1e-10
    
    if (log):
      maxY *= 10
    else:
      maxY *= 1.3
    
    for massID in range(len(mass_list)):
        if(float(mass_list[massID]) == 4.5):
          continue
        temp=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_"+hName)
        if (hName == 'm_HN'):
          temp.SetTitle('invariant mass of #mu_{2} and #mu_{3}; invariant mass of #mu_{2} and #mu_{3} [GeV]')
        temp.SetLineColor(1+massID)
        temp.SetLineWidth(3)
        if (h == 0):
            legend.AddEntry(temp,"HN "+mass_list[massID]+" GeV","l")
            noBackLegend.AddEntry(temp,"HN "+mass_list[massID]+" GeV","l")
        if (massID == 0):
          temp.GetYaxis().SetRangeUser(minY,maxY)
        if (noError):
          temp.Draw("histsame")
        else:
          temp.Draw("hesame")
    for b in backID:
        if (onlySig):
          continue
        temp=merged_file.Get("gamma"+str(b)+"/h_"+hName)
        temp.SetLineColor(background_colors[b])
        temp.SetLineWidth(3)
        temp.SetLineStyle(2)
        if (h == 0):
            #legend.AddEntry(temp,"gamma"+str(b),"l")
            legend.AddEntry(temp,background_processes[b],"l")
        if (noError):
          temp.Draw("histsame")
        else:
          temp.Draw("hesame")
    if (onlySig):
        noBackLegend.Draw()
    else:
        legend.Draw()
    c_temp.SaveAs("plots_muon/"+hName+".png")
    c_temp.SaveAs("plots_muon/kinematics/"+hName+".pdf")

'''
c_dalitz_HN = TCanvas("c_dalitz_HN","c_dalitz_HN", (len(mass_list)+len(backID))*800, 800)
c_dalitz_HN.Divide((len(mass_list)+len(backID)),1)
for bID in range(len(backID)):
  b = backID[bID]
  temp=merged_file.Get("gamma"+str(b)+"/h_dalitz_HN")
  temp.SetTitle(background_processes[b])
  c_dalitz_HN.cd(bID+1)
  temp.Draw("colz")
for massID in range(len(mass_list)):
  temp=merged_file.Get("HN_mass_"+mass_list[massID]+"/h_dalitz_HN")
  temp.SetTitle("HN "+mass_list[massID])
  c_dalitz_HN.cd(len(backID)+massID+1)
  temp.Draw("colz")
c_dalitz_HN.SaveAs("plots_muon/dalitz_HN.png")
'''







hists2D = ['nEvents_Maj_couplingMuTau','nEvents_Dir_couplingMuTau','BR_Maj_couplingMuTau','BR_Dir_couplingMuTau']

contours = np.array([1.])

for h in range(len(hists2D)):
  hName = hists2D[h]
  for massID in range(len(mass_list)):
    c_temp = TCanvas('c_'+hName+'_m'+mass_list[massID],"canvas title", 800, 800) 
    c_temp.SetTopMargin(0.09)
    c_temp.SetBottomMargin(0.1)
    c_temp.SetLeftMargin(0.1)
    c_temp.SetRightMargin(0.14)
    c_temp.SetLogz(1)
    temp=merged_file.Get("HN_mass_"+mass_list[massID]+'/h_'+hName)
    temp.SetTitle('HN '+mass_list[massID]+': '+temp.GetTitle())
    if (h<2):
      temp.DrawCopy("colz")
      temp.SetContour(1,contours)
      temp.Draw("cont3 same")
      temp.SetLineColor(1)
      temp.SetLineWidth(3)
    else:
      temp.Draw("colz")
    c_temp.SaveAs('plots_muon/2D/'+hName+'_m'+mass_list[massID]+'.png')
    c_temp.SaveAs('plots_muon/kinematics/'+hName+'_m'+mass_list[massID]+'.pdf')
