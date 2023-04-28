#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import numpy as np
import pandas as pd


# In[2]:


ROOT.ROOT.EnableImplicitMT(64)


# In[245]:

ROOT.gInterpreter.Declare('''

bool singleMuCheck(float muPt, float muEta, int muQual, float minPt, float maxEta, int minQual){
  if (muPt < minPt) return false;
  if (abs(muEta) > maxEta) return false;
  if (muQual < minQual) return false;
  return true;
}

float L1deltaR(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2){
    const Double_t muonMass = 0.105658;
    const ROOT::Math::PtEtaPhiMVector mu1(pt1, eta1, phi1, muonMass);
    const ROOT::Math::PtEtaPhiMVector mu2(pt2, eta2, phi2, muonMass);
    return ROOT::Math::VectorUtil::DeltaR(mu1, mu2);
}

bool flag(ROOT::VecOps::RVec<float> L1mu_pt, ROOT::VecOps::RVec<float> L1mu_eta, ROOT::VecOps::RVec<float> L1mu_phi, ROOT::VecOps::RVec<int> L1mu_charge, ROOT::VecOps::RVec<int> L1mu_quality,
                float minPt, float maxEta, int minQuality, int sign, float deltaRmax){
int sum = 0;
for (int i = 0; i < L1mu_pt.size(); i++) {
  if (!singleMuCheck(L1mu_pt[i], L1mu_eta[i], L1mu_quality[i], minPt, maxEta, minQuality)) continue;
  for (int j = i+1; j < L1mu_pt.size(); j++) {
    if (L1mu_charge[i]*L1mu_charge[j] != sign) continue;
    if (!singleMuCheck(L1mu_pt[j], L1mu_eta[j], L1mu_quality[j], minPt, maxEta, minQuality)) continue;
    if (L1deltaR(L1mu_pt[i], L1mu_eta[i], L1mu_phi[i], L1mu_pt[j], L1mu_eta[j], L1mu_phi[j]) > deltaRmax) continue;
    sum++;
  }
}
return (sum > 0);
}

float dPtL1offline(ROOT::VecOps::RVec<float> L1mu_pt, ROOT::VecOps::RVec<float> L1mu_eta, ROOT::VecOps::RVec<float> L1mu_phi, ROOT::VecOps::RVec<int> L1mu_charge,
                   float offMu_px, float offMu_py, float offMu_pz, float offMu_charge){
  const Double_t muonMass = 0.105658;
  const ROOT::Math::PxPyPzMVector mu1(offMu_px, offMu_py, offMu_pz, muonMass);
  float dRmin = 999;
  int index = -1;
  for (int i = 0; i < L1mu_pt.size(); i++) {
    if (L1mu_charge[i]*offMu_charge < 0) continue;
    const ROOT::Math::PtEtaPhiMVector mu2(L1mu_pt[i], L1mu_eta[i], L1mu_phi[i], muonMass);
    float dR = ROOT::Math::VectorUtil::DeltaR(mu1, mu2);
    if (dR < dRmin) {
      dRmin = dR;
      index = i;
    }
  }
  float val = 0.;
  if (index > -1) val = (mu1.Pt() - L1mu_pt[index])/mu1.Pt();
  else val = -999; 
  return val;
}

float dRminL1offline(ROOT::VecOps::RVec<float> L1mu_pt, ROOT::VecOps::RVec<float> L1mu_eta, ROOT::VecOps::RVec<float> L1mu_phi, ROOT::VecOps::RVec<int> L1mu_charge,
                   float offMu_px, float offMu_py, float offMu_pz, float offMu_charge){
  const Double_t muonMass = 0.105658;
  const ROOT::Math::PxPyPzMVector mu1(offMu_px, offMu_py, offMu_pz, muonMass);
  float dRmin = 999;
  for (int i = 0; i < L1mu_pt.size(); i++) {
    if (L1mu_charge[i]*offMu_charge < 0) continue;
    const ROOT::Math::PtEtaPhiMVector mu2(L1mu_pt[i], L1mu_eta[i], L1mu_phi[i], muonMass);
    float dR = ROOT::Math::VectorUtil::DeltaR(mu1, mu2);
    if (dR < dRmin) {
      dRmin = dR;
    }
  }
  return dRmin;
}

'''
)

#!ls /lustre/cms/store/user/vmastrap/ParkingDoubleMuonLowMass?/MuMuTrigger_?_12feb23/*/*/*.root > listMuMuTrigger_2022G_DMLM_12feb23.txt


# In[247]:

with open('test.txt') as f:
    allFilesG = f.readlines()
for i in range(len(allFilesG)):
    allFilesG[i] = allFilesG[i].replace("\n", "")

# In[254]:

data = ROOT.RDataFrame("rootuple/ntuple", allFilesG)

# In[5]:

data = data.Define("J_pt", "sqrt(J_px*J_px + J_py*J_py)")
data = data.Define("J_p", "sqrt(J_px*J_px + J_py*J_py + J_pz*J_pz)")
data = data.Define("J_p1", "sqrt(J_px1*J_px1 + J_py1*J_py1 + J_pz1*J_pz1)")
data = data.Define("J_p2", "sqrt(J_px2*J_px2 + J_py2*J_py2 + J_pz2*J_pz2)")
data = data.Define("J_eta", "atanh(J_pz/J_p)")

data = data.Define("ptdiff1", "dPtL1offline(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, J_px1, J_py1, J_pz1, J_charge1)")
data = data.Define("ptdiff2", "dPtL1offline(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, J_px2, J_py2, J_pz2, J_charge2)")
data = data.Define("drmin1", "dRminL1offline(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, J_px1, J_py1, J_pz1, J_charge1)")
data = data.Define("drmin2", "dRminL1offline(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, J_px2, J_py2, J_pz2, J_charge2)")
data = data.Define("l1pt1", "J_pt1 - ptdiff1")
data = data.Define("l1pt2", "J_pt2 - ptdiff2")

data = data.Define("J_deltaEta", "abs(J_eta1 - J_eta2)")
data = data.Define("J_deltaPhi", "Double_t angle = abs(J_phi1 - J_phi2); if (angle > ROOT::Math::Pi()) angle = 2*ROOT::Math::Pi() - angle; return angle;")
data = data.Define("J_deltaR", "sqrt(J_deltaEta*J_deltaEta + J_deltaPhi*J_deltaPhi)")

# In[22]:

data = data.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
data = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
data = data.Filter("HLTLowMassInclusive == 1")


hist1 = data.Histo1D(("", ";(off-L1)/off p_{T}(#mu);Candidates / 0.01", 200, -1, 1), "ptdiff1")
hist2 = data.Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff2")

hist10 = data.Filter("l1pt1 < 4.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff1")
hist20 = data.Filter("l1pt2 < 4.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff2")
hist11 = data.Filter("l1pt1 < 7.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff1")
hist21 = data.Filter("l1pt2 < 7.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff2")
hist12 = data.Filter("l1pt1 < 10.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff1")
hist22 = data.Filter("l1pt2 < 10.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff2")
#hist13 = data.Filter("l1pt1 < 3.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff1")
#hist23 = data.Filter("l1pt2 < 3.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff2")
#hist14 = data.Filter("l1pt1 < 4.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff1")
#hist24 = data.Filter("l1pt2 < 4.1").Histo1D(("", ";#Delta p_{T} [GeV];Candidates / 50 MeV", 200, -1, 1), "ptdiff2")

hist1.Add(hist2.GetPtr())
hist10.Add(hist20.GetPtr())
hist11.Add(hist21.GetPtr())
hist12.Add(hist22.GetPtr())
#hist13.Add(hist23.GetPtr())
#hist14.Add(hist24.GetPtr())

c0 = ROOT.TCanvas("c0", "c0", 1200, 1200)

legend = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
legend.AddEntry(hist1.GetPtr(), "Any p_{T,L1}(#mu)", "l")
legend.AddEntry(hist10.GetPtr(), "p_{T,L1}(#mu) #leq 4", "f")
legend.AddEntry(hist11.GetPtr(), "p_{T,L1}(#mu) #leq 7", "f")
legend.AddEntry(hist12.GetPtr(), "p_{T,L1}(#mu) #leq 10", "f")
#legend.AddEntry(hist13.GetPtr(), "p_{T}^{L1} #le 3", "f")
#legend.AddEntry(hist14.GetPtr(), "p_{T}^{L1} #le 4", "f")

hist1.SetStats(0)

hist10.SetFillColor(2)
hist11.SetFillColor(3)
hist12.SetFillColor(4)
#hist13.SetFillColor(5)
#hist14.SetFillColor(6)

hist1.Draw("")
#hist14.Draw("same")
#hist13.Draw("same")
hist12.Draw("same")
hist11.Draw("same")
hist10.Draw("same")

legend.Draw("same")

c0.SetLeftMargin(0.15)
c0.Print("ptResolutionTrigger.png")
c0.Clear()
c0.Close()


#

hist2d = data.Histo2D(("", ";p_{T}(reco) [GeV];p_{T}(L1) [GeV]", 100, 0, 20, 100, 0, 20), "J_pt1", "l1pt1")
hist2dB = data.Histo2D(("", ";p_{T}(reco) [GeV];p_{T}(L1) [GeV]", 100, 0, 20, 100, 0, 20), "J_pt2", "l1pt2")

hist2d.Add(hist2dB.GetPtr())

c0 = ROOT.TCanvas("", "", 1200, 1200)

hist2d.SetStats(0)
hist2d.Draw("colz")
c0.SetLeftMargin(0.15)
c0.SetRightMargin(0.15)
c0.Print("ptL1vsOffTrigger.png")
c0.Close()
c0.Clear()
exit()
#
hist1 = data.Histo1D(("", ";#DeltaR;Candidates / 0.01", 100, 0., 1.0), "drmin1")
hist2 = data.Histo1D(("", ";#DeltaR;Candidates / 0.01", 100, 0., 1.0), "drmin2")

hist1.Add(hist2.GetPtr())

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

hist1.SetStats(0)
hist1.Scale(1./hist1.Integral())
hist1.Draw("hist")

c1.SetLeftMargin(0.15)
c1.Print("drminBarrel.png")
c1.Clear()
c1.Close()

