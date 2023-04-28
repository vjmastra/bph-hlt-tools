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

float L1deltaEta(float eta1, float eta2){
    return abs(eta1 - eta2);
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

bool flagEta(ROOT::VecOps::RVec<float> L1mu_pt, ROOT::VecOps::RVec<float> L1mu_eta, ROOT::VecOps::RVec<float> L1mu_phi, ROOT::VecOps::RVec<int> L1mu_charge, ROOT::VecOps::RVec<int> L1mu_quality,
                float minPt, float maxEta, int minQuality, int sign, float deltaEtamax){
int sum = 0;
for (int i = 0; i < L1mu_pt.size(); i++) {
  if (!singleMuCheck(L1mu_pt[i], L1mu_eta[i], L1mu_quality[i], minPt, maxEta, minQuality)) continue;
  for (int j = i+1; j < L1mu_pt.size(); j++) {
    if (L1mu_charge[i]*L1mu_charge[j] != sign) continue;
    if (!singleMuCheck(L1mu_pt[j], L1mu_eta[j], L1mu_quality[j], minPt, maxEta, minQuality)) continue;
    if (L1deltaEta(L1mu_eta[i], L1mu_eta[j]) > deltaEtamax) continue;
    sum++;
  }
}
return (sum > 0);
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
#data = data.Define("J_eta1", "atanh(J_pz1/J_p1)")
#data = data.Define("J_eta2", "atanh(J_pz2/J_p2)")

data = data.Define("J_deltaEta", "abs(J_eta1 - J_eta2)")
data = data.Define("J_deltaPhi", "Double_t angle = abs(J_phi1 - J_phi2); if (angle > ROOT::Math::Pi()) angle = 2*ROOT::Math::Pi() - angle; return angle;")
data = data.Define("J_deltaR", "sqrt(J_deltaEta*J_deltaEta + J_deltaPhi*J_deltaPhi)")


data = data.Define("L10er1p4SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.4, 12, -1, 1.4)")
data = data.Define("L10er1p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 1.4)")
data = data.Define("L10er1p6SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.6, 12, -1, 1.4)")
data = data.Define("L10er1p7SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.7, 12, -1, 1.4)")
data = data.Define("L10er1p8SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.8, 12, -1, 1.4)")
data = data.Define("L10er1p9SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.9, 12, -1, 1.4)")
data = data.Define("L10er2p0SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 2.0, 12, -1, 1.4)")

data = data.Define("L10er2p0SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 2.0, 12, -1, 1.2)")
data = data.Define("L10er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 2.0, 12, -1, 1.6)")
data = data.Define("L10er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 2.0, 12, -1, 2.0)")
data = data.Define("L10er2p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 2.0, 12, -1, 2.4)")
data = data.Define("L10er2p0SQOSdR2p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 2.0, 12, -1, 2.8)")

data = data.Define("L11er2p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 1., 2.0, 12, -1, 2.4)")
data = data.Define("L12er2p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 2.4)")
data = data.Define("L13er2p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 2.4)")
data = data.Define("L14er2p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4., 2.0, 12, -1, 2.4)")
data = data.Define("L15er2p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 5., 2.0, 12, -1, 2.4)")

data = data.Define("L13er1p8SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.8, 12, -1, 2.4)")
data = data.Define("L13er1p6SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.6, 12, -1, 2.4)")
data = data.Define("L13er1p4SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.4, 12, -1, 2.4)")
data = data.Define("L13er1p2SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.2, 12, -1, 2.4)")
data = data.Define("L13er1p0SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.0, 12, -1, 2.4)")

data = data.Define("L11er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 1., 2.0, 12, -1, 2.0)")
data = data.Define("L12er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 2.0)")
#data = data.Define("L13er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 2.0)")
data = data.Define("L14er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4., 2.0, 12, -1, 2.0)")
data = data.Define("L15er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 5., 2.0, 12, -1, 2.0)")

data = data.Define("L13er1p8SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.8, 12, -1, 2.0)")
data = data.Define("L13er1p6SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.6, 12, -1, 2.0)")
#data = data.Define("L13er1p4SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.4, 12, -1, 2.0)")
data = data.Define("L13er1p2SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.2, 12, -1, 2.0)")
data = data.Define("L13er1p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.0, 12, -1, 2.0)")

data = data.Define("L10er1p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 1.6)")
data = data.Define("L10er1p5SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 1.8)")
data = data.Define("L10er1p5SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 2.0)")
data = data.Define("L10er1p5SQOSdR2p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 2.2)")
data = data.Define("L10er1p5SQOSdR2p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 2.4)")
data = data.Define("L10er1p5SQOSdR2p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 2.6)")
data = data.Define("L10er1p5SQOSdR2p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 2.8)")
data = data.Define("L10er1p5SQOSdR3p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 3.0)")

data = data.Define("L10er1p2SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.2, 12, -1, 1.2)")
data = data.Define("L10er1p2SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.2, 12, -1, 1.4)")
data = data.Define("L10er1p2SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.2, 12, -1, 1.6)")
data = data.Define("L10er1p2SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.2, 12, -1, 1.8)")
data = data.Define("L10er1p2SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.2, 12, -1, 2.0)")
data = data.Define("L10er1p2SQOSdR2p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.2, 12, -1, 2.2)")

data = data.Define("L12er1p5SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.2)")
data = data.Define("L12er1p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.4)")
data = data.Define("L12er1p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.6)")
data = data.Define("L12er1p5SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.8)")
data = data.Define("L12er1p5SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 2.0)")
data = data.Define("L12er1p5SQOSdR2p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.2, 12, -1, 2.2)")

data = data.Define("L13er1p5SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.2)")
data = data.Define("L13er1p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.4)")
data = data.Define("L13er1p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.6)")
data = data.Define("L13er1p5SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 1.8)")
data = data.Define("L13er1p5SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 1.5, 12, -1, 2.0)")

data = data.Define("L13er2p0SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 1.2)")
data = data.Define("L13er2p0SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 1.4)")
#data = data.Define("L13er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 1.6)")
data = data.Define("L13er2p0SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 1.8)")
#data = data.Define("L13er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 2.0)")

data = data.Define("L14er2p5SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.2)")
data = data.Define("L14er2p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.4)")
data = data.Define("L14er2p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.6)")
data = data.Define("L14er2p5SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.8)")
data = data.Define("L14er2p5SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 2.0)")

data = data.Define("L11er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 1., 2.0, 12, -1, 1.6)")
data = data.Define("L12er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 1.6)")
data = data.Define("L13er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 1.6)")
data = data.Define("L14er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4., 2.0, 12, -1, 1.6)")
data = data.Define("L15er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 5., 2.0, 12, -1, 1.6)")

data = data.Define("L13er1p8SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.8, 12, -1, 1.6)")
data = data.Define("L13er1p6SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.6, 12, -1, 1.6)")
data = data.Define("L13er1p4SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.4, 12, -1, 1.6)")
data = data.Define("L13er1p2SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.2, 12, -1, 1.6)")
data = data.Define("L13er1p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.0, 12, -1, 1.6)")

data = data.Define("L10er2p0SQOSdEta1p6", "flagEta(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0.0, 2.0, 12, -1, 1.6)")
data = data.Define("L10er2p0SQOSdEta1p5", "flagEta(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0.0, 2.0, 12, -1, 1.5)")

data = data.Define("L15er2p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 5., 2.5, 12, -1, 1.6)")
data = data.Define("L10er1p5SQOSdE1p2", "flagEta(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0.0, 1.5, 12, -1, 1.2)")
data = data.Define("L10er1p4SQOSdE1p2", "flagEta(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0.0, 1.4, 12, -1, 1.2)")

data = data.Define("L13er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 2.0)")
data = data.Define("L13er1p4SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.4, 12, -1, 2.0)")
data = data.Define("L10er1p4SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.4, 12, -1, 2.0)")
data = data.Define("L10er1p4SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.4, 12, -1, 1.6)")

#

testseeds = ["L10er1p5SQOSdR1p4", "L10er1p7SQOSdR1p4", "L10er2p0SQOSdR1p4", \
"L10er1p5SQOSdR1p8", "L10er1p5SQOSdR2p2", "L10er1p5SQOSdR2p6", \
"L10er2p0SQOSdR1p2", "L10er2p0SQOSdR1p6", "L10er2p0SQOSdR2p0", "L10er2p0SQOSdR2p4", "L10er2p0SQOSdR2p8", \
"L11er2p0SQOSdR2p4", "L12er2p0SQOSdR2p4", "L13er2p0SQOSdR2p4", "L14er2p0SQOSdR2p4", \
"L13er2p0SQOSdR2p4", "L13er1p8SQOSdR2p4", "L13er1p6SQOSdR2p4", "L13er1p4SQOSdR2p4", "L13er1p2SQOSdR2p4", \
"L11er2p0SQOSdR2p0", "L12er2p0SQOSdR2p0", "L13er2p0SQOSdR2p0", "L14er2p0SQOSdR2p0", \
"L13er2p0SQOSdR2p0", "L13er1p8SQOSdR2p0", "L13er1p6SQOSdR2p0", "L13er1p4SQOSdR2p0", "L13er1p2SQOSdR2p0", \
]

testseeds2 = ["L10er2p0SQOSdEta1p6", "L10er2p0SQOSdEta1p5", "L13er2p0SQOSdR1p4", "L14er2p5SQOSdR1p2", "L10er1p4SQOSdR1p4", \
"L11er2p0SQOSdR1p6", "L12er2p0SQOSdR1p6", "L13er2p0SQOSdR1p6", "L14er2p0SQOSdR1p6", \
"L13er2p0SQOSdR1p6", "L13er1p8SQOSdR1p6", "L13er1p6SQOSdR1p6", "L13er1p4SQOSdR1p6", "L13er1p2SQOSdR1p6", \
]

testseeds3 = ["L10er1p5SQOSdR1p4", "L10er2p0SQOSdR1p6", "L10er2p0SQOSdR2p0", "L10er2p0SQOSdR2p4", \
"L13er2p0SQOSdR2p4", "L13er1p8SQOSdR2p4", "L13er1p6SQOSdR2p4", \
"L13er2p0SQOSdR2p0", "L13er1p8SQOSdR2p0", "L13er1p6SQOSdR2p0", \
"L13er2p0SQOSdR1p6", "L13er1p8SQOSdR1p6", "L13er1p6SQOSdR1p6"]

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05") 
#denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
#denom2 = denom.Filter("J_mass > 4.5")

#dd1 = denom.Count().GetValue()
#dd2 = denom2.Count().GetValue()

denomA = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
denomB = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass > 4.5 & J_mass < 8.5")
denomC = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 5.3")
denomD = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass > 4.8 & J_mass < 5.8")

ddA = denomA.Count().GetValue()
ddB = denomB.Count().GetValue()
ddC = denomC.Count().GetValue()
ddD = denomD.Count().GetValue()

#optseeds = ["L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0", "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 + L13er2p0SQOSdR1p6 + L10er1p5SQOSdE1p2 > 0", \
#"L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0", "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0"]

#seeds0 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0"
#seeds1 = "L13er2p0SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds2 = "L10er1p4SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds3 = "L13er1p4SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds4 = "L10er1p4SQOSdR1p6 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds5 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0"

seeds0 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0"
seeds1 = "L10er1p4SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
seeds2 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0"
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L14er2p0SQOSdR1p6 > 0"
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0"

#optseeds = [seeds0, seeds1, seeds2, seeds3, seeds4, seeds5]
optseeds = [seeds3]

for seed in optseeds:
  numA = denomA.Filter(seed).Count().GetValue()
  numB = denomB.Filter(seed).Count().GetValue()
  numC = denomC.Filter(seed).Count().GetValue()
  numD = denomD.Filter(seed).Count().GetValue()
  print(numA/ddA, numB/ddB, numC/ddC, numD/ddD)

exit()

with open("efficiency3.txt", "r+") as fo:

#  print("Seed name", "18+X offline eff", "in [4.5, 8.5]", "X only offline eff", "in [4.5, 8.5]", file=fo)
  print("Seed name", "18+X offline eff in [0., 5.3]", "18+X only offline eff in [4.8, 5.8]", file=fo)

  for seed in testseeds3:
    #num1a = denom.Filter(seed+" + L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0").Count().GetValue()
    #num2a = denom2.Filter(seed+" + L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0").Count().GetValue()
    #num1b = denom.Filter(seed + " > 0").Count().GetValue()
    #num2b = denom2.Filter(seed + " > 0").Count().GetValue()
    #print(seed, num1a/dd1, num2a/dd2, num1b/dd1, num2b/dd2, file=fo)
    numA = denomA.Filter(seed+" + L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0").Count().GetValue()
    numB = denomB.Filter(seed+" + L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0").Count().GetValue()
    print(seed, numA/ddA, numB/ddB)

  print("\n", file=fo)

'''
print(data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1").Count().GetValue())
print(data.Filter("L10er1p5SQOSdR1p2 == 1").Count().GetValue())
print(data.Filter("L10er1p5dR == 1").Count().GetValue())
print(data.Filter("L14er2p5SQOSdR1p4 == 1").Count().GetValue())
print(data.Filter("L14dR == 1").Count().GetValue())
'''
# In[22]:

# In[35]:
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05") 
#.Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5") 

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
num0 = denom.Filter("L10er1p5SQOSdR1p2 + L14er2p5SQOSdR1p4 > 0")
num1 = denom.Filter("L10er1p5SQOSdR1p2 + L14er2p5SQOSdR1p4 + L10er1p6SQOSdR1p2 > 0")
num2 = denom.Filter("L10er1p5SQOSdR1p2 + L14er2p5SQOSdR1p4 + L10er1p7SQOSdR1p2 > 0")
num3 = denom.Filter("L10er1p5SQOSdR1p2 + L14er2p5SQOSdR1p4 + L10er1p8SQOSdR1p2 > 0")
num4 = denom.Filter("L10er1p5SQOSdR1p2 + L14er2p5SQOSdR1p4 + L10er1p9SQOSdR1p2 > 0")
num5 = denom.Filter("L10er1p5SQOSdR1p2 + L14er2p5SQOSdR1p4 + L10er2p0SQOSdR1p2 > 0")

xmin = 0.
xmax = 8.5
nbins = int((xmax-xmin)/0.02)

denhist = denom.Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
num0hist = num0.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num1hist = num1.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num2hist = num2.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num3hist = num3.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num4hist = num4.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num5hist = num4.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")

ratio0 = ROOT.TEfficiency(num0hist.GetPtr(), denhist.GetPtr())
ratio0.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio0.SetConfidenceLevel(0.68)

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

ratio2 = ROOT.TEfficiency(num2hist.GetPtr(), denhist.GetPtr())
ratio2.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio2.SetConfidenceLevel(0.68)

ratio3 = ROOT.TEfficiency(num3hist.GetPtr(), denhist.GetPtr())
ratio3.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio3.SetConfidenceLevel(0.68)

ratio4 = ROOT.TEfficiency(num4hist.GetPtr(), denhist.GetPtr())
ratio4.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio4.SetConfidenceLevel(0.68)

ratio5 = ROOT.TEfficiency(num5hist.GetPtr(), denhist.GetPtr())
ratio5.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio5.SetConfidenceLevel(0.68)

# In[36]:

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num0hist.SetFillColor(2)
num1hist.SetFillColor(3)
num2hist.SetFillColor(4)
num3hist.SetFillColor(5)
num4hist.SetFillColor(6)
num5hist.SetFillColor(7)

legend = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(), "All L1 seeds", "l")
legend.AddEntry(num0hist.GetPtr(), "Only 2018 seeds", "f")
legend.AddEntry(num1hist.GetPtr(), "2018 + DMu0_er1p6_SQ_OS_dR_Max1p2", "f")
legend.AddEntry(num2hist.GetPtr(), "2018 + DMu0_er1p7_SQ_OS_dR_Max1p2", "f")
legend.AddEntry(num3hist.GetPtr(), "2018 + DMu0_er1p8_SQ_OS_dR_Max1p2", "f")
legend.AddEntry(num4hist.GetPtr(), "2018 + DMu0_er1p9_SQ_OS_dR_Max1p2", "f")
legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max1p2", "f")

denhist.GetYaxis().SetRangeUser(1e4, 3e7)

denhist.Draw()
num5hist.Draw("same")
num4hist.Draw("same")
num3hist.Draw("same")
num2hist.Draw("same")
num1hist.Draw("same")
num0hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";m(#mu#mu) [GeV];Ratio 2018+X/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)

ratio0.SetLineColor(2)
ratio1.SetLineColor(3)
ratio2.SetLineColor(4)
ratio3.SetLineColor(5)
ratio4.SetLineColor(6)
ratio5.SetLineColor(7)

emptyHist.Draw("")
ratio0.Draw("esame")
ratio1.Draw("esame")
ratio2.Draw("esame")
ratio3.Draw("esame")
ratio4.Draw("esame")
ratio5.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("L1check_looserEta.png")
canvas1.Clear("")
canvas1.Close()

'''

