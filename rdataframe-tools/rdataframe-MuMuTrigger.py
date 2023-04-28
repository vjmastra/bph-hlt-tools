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

float JdeltaRL1(ROOT::VecOps::RVec<float> L1mu_pt, ROOT::VecOps::RVec<float> L1mu_eta, ROOT::VecOps::RVec<float> L1mu_phi, ROOT::VecOps::RVec<int> L1mu_charge){
float dRmin = 999;
for (int i = 0; i < L1mu_pt.size(); i++) {
  for (int j = i+1; j < L1mu_pt.size(); j++) {
    if (L1mu_charge[i]*L1mu_charge[j] > 0) continue;
    float dRval = L1deltaR(L1mu_pt[i], L1mu_eta[i], L1mu_phi[i], L1mu_pt[j], L1mu_eta[j], L1mu_phi[j]);
    if (dRval < dRmin) dRmin = dRval;
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
#data = data.Define("J_eta1", "atanh(J_pz1/J_p1)")
#data = data.Define("J_eta2", "atanh(J_pz2/J_p2)")

data = data.Define("J_deltaEta", "abs(J_eta1 - J_eta2)")
data = data.Define("J_deltaPhi", "Double_t angle = abs(J_phi1 - J_phi2); if (angle > ROOT::Math::Pi()) angle = 2*ROOT::Math::Pi() - angle; return angle;")
data = data.Define("J_deltaR", "sqrt(J_deltaEta*J_deltaEta + J_deltaPhi*J_deltaPhi)")

data = data.Define("J_deltaR_L1", "JdeltaRL1(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge)")

#data = data.Define("L10er1p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 1.4)")
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

data = data.Define("L11er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 1., 2.0, 12, -1, 1.6)")
data = data.Define("L12er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 2., 2.0, 12, -1, 1.6)")
#data = data.Define("L13er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 1.6)")
#data = data.Define("L14er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4., 2.0, 12, -1, 1.6)")
data = data.Define("L15er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 5., 2.0, 12, -1, 1.6)")

data = data.Define("L13er1p8SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.8, 12, -1, 1.6)")
data = data.Define("L13er1p6SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.6, 12, -1, 1.6)")
data = data.Define("L13er1p4SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.4, 12, -1, 1.6)")
data = data.Define("L13er1p2SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.2, 12, -1, 1.6)")
data = data.Define("L13er1p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.0, 12, -1, 1.6)")

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

data = data.Define("L10er1p5SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 1.2)")
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

#data = data.Define("L14er2p5SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.2)")
data = data.Define("L14er2p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.4)")
data = data.Define("L14er2p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.6)")
data = data.Define("L14er2p5SQOSdR1p8", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.8)")
data = data.Define("L14er2p5SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 2.0)")

#Dima seeds
data = data.Define("L10er1p5SQOSdE1p2", "flagEta(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0.0, 1.5, 12, -1, 1.2)")
data = data.Define("L10er1p4SQOSdE1p2", "flagEta(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0.0, 1.4, 12, -1, 1.2)")
data = data.Define("L14er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4., 2.0, 12, -1, 1.6)")
data = data.Define("L15er2p5SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 5., 2.5, 12, -1, 1.6)")

#Vincenzo seeds
data = data.Define("L10er1p5SQOSdR1p4", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.5, 12, -1, 1.4)") #2018
data = data.Define("L14er2p5SQOSdR1p2", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4.0, 2.5, 12, -1, 1.2)") #2018
data = data.Define("L13er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 1.6)") #new
#data = data.Define("L14er2p0SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 4., 2.0, 12, -1, 1.6)") #backup - shared with Dima

#Riccardo seeds1
data = data.Define("L13er2p0SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 2.0, 12, -1, 2.0)")
data = data.Define("L13er1p4SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 3., 1.4, 12, -1, 2.0)")
data = data.Define("L10er1p4SQOSdR2p0", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.4, 12, -1, 2.0)")
data = data.Define("L10er1p4SQOSdR1p6", "flag(L1mu_pt, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality, 0., 1.4, 12, -1, 1.6)")

#seeds0 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0"
#seeds1 = "L13er2p0SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds2 = "L10er1p4SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds3 = "L13er1p4SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
#seeds4 = "L10er1p4SQOSdR1p6 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
seeds5 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0"

seeds0 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0"
seeds1 = "L10er1p4SQOSdR2p0 + L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0"
seeds2 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0"
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L14er2p0SQOSdR1p6 > 0"
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0"

#

#print(data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1").Count().GetValue())
#print(data.Filter("L10er1p5SQOSdR1p4 == 1").Count().GetValue())
#print(data.Filter("L10er1p5dR == 1").Count().GetValue())
#print(data.Filter("L14er2p5SQOSdR1p2 == 1").Count().GetValue())
#print(data.Filter("L14dR == 1").Count().GetValue())

# In[22]:

'''
data.Filter("mu1HLTmatched + mu2HLTmatched > 0").Count().GetValue()


# In[23]:

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
inclusivecut = mumudata.Filter("HLTLowMassInclusive == 1")
displacedcut = mumudata.Filter("HLTLowMassDisplaced == 1")

mumuhist = mumudata.Histo1D(("", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", 575, 0., 11.5), "J_mass")
incluhist = inclusivecut.Histo1D(("J_mass", "dimuon mass", 575, 0., 11.5), "J_mass")
displhist = displacedcut.Histo1D(("J_mass", "dimuon mass", 575, 0., 11.5), "J_mass")

# In[24]:

mumuhist.SetStats(0)

incluhist.SetFillColor(2)
displhist.SetFillColor(3)

string = ROOT.TLatex()
string.SetTextSize(0.04)

line = ROOT.TLine(8.5, 1, 8.5, 3e8)
line.SetLineColor(ROOT.kBlue)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)

legend = ROOT.TLegend(0.4, 0.75, 0.68, 0.89)
legend.AddEntry(mumuhist.GetPtr(), "Dimuon candidates", "l")
legend.AddEntry(incluhist.GetPtr(), "w. inclusive #mu#mu trigger fired", "f")
legend.AddEntry(displhist.GetPtr(), "w. displaced #mu#mu trigger fired", "f")

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)
canvas1.SetLogy()

mumuhist.GetYaxis().SetRangeUser(1, 3e8)
mumuhist.Draw()
incluhist.Draw("same")
displhist.Draw("same")

line.Draw("same")
legend.Draw("same")

canvas1.Print("./DMLM_HLT.png")
canvas1.Clear()
canvas1.Close()
'''
# In[35]:
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05") 
#.Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5") 

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
num0 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0")
num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p6SQOSdR1p4 > 0")
num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p7SQOSdR1p4 > 0")
num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p8SQOSdR1p4 > 0")
num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p9SQOSdR1p4 > 0")
num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR1p4 > 0")

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
legend.AddEntry(num1hist.GetPtr(), "2018 + DMu0_er1p6_SQ_OS_dR_Max1p4", "f")
legend.AddEntry(num2hist.GetPtr(), "2018 + DMu0_er1p7_SQ_OS_dR_Max1p4", "f")
legend.AddEntry(num3hist.GetPtr(), "2018 + DMu0_er1p8_SQ_OS_dR_Max1p4", "f")
legend.AddEntry(num4hist.GetPtr(), "2018 + DMu0_er1p9_SQ_OS_dR_Max1p4", "f")
legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max1p4", "f")

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


# In[37]:
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
#.Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5") 

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
num0 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0")
num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR1p4 > 0")
num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR1p6 > 0")
num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR1p8 > 0")
num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR2p0 > 0")
num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR2p2 > 0")

#num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR1p6 > 0")
#num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR1p8 > 0")
#num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR2p2 > 0")
#num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR2p6 > 0")
#num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er1p5SQOSdR3p0 > 0")

#num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR1p2 > 0")
#num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR1p6 > 0")
#num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR2p0 > 0")
#num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR2p4 > 0")
#num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR2p8 > 0")

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
legend.AddEntry(num1hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max1p4", "f")
legend.AddEntry(num2hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num3hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max1p8", "f")
legend.AddEntry(num4hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max2p0", "f")
legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max2p2", "f")

#legend.AddEntry(num1hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max1p6", "f")
#legend.AddEntry(num2hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max1p8", "f")
#legend.AddEntry(num3hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max2p2", "f")
#legend.AddEntry(num4hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max2p6", "f")
#legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er1p5_SQ_OS_dR_Max3p0", "f")

#legend.AddEntry(num1hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max1p2", "f")
#legend.AddEntry(num2hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max1p6", "f")
#legend.AddEntry(num3hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num4hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max2p4", "f")
#legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max2p8", "f")

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
ratio5.Draw("esame")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio2.Draw("esame")
ratio1.Draw("esame")
ratio0.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("L1check_0er1p5_dRscan.png")
#canvas1.Print("L1check_0er1p5_dRlargeScan.png")
#canvas1.Print("L1check_0er2p0_dRlargeScan.png")
canvas1.Clear("")
canvas1.Close()

'''
# In[37]:
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
num0 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0")
#num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L14er2p0SQOSdR2p0 > 0")
#num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR2p0 > 0")
#num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L12er2p0SQOSdR2p0 > 0")
#num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L11er2p0SQOSdR2p0 > 0")
#num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR2p0 > 0")

num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L14er2p0SQOSdR1p6 > 0")
num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0")
num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L12er2p0SQOSdR1p6 > 0")
num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L11er2p0SQOSdR1p6 > 0")
num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0SQOSdR1p6 > 0")

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
#legend.AddEntry(num1hist.GetPtr(), "2018 + DMu4_er2p0_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num2hist.GetPtr(), "2018 + DMu3_er2p0_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num3hist.GetPtr(), "2018 + DMu2_er2p0_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num4hist.GetPtr(), "2018 + DMu1_er2p0_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max2p0", "f")

legend.AddEntry(num1hist.GetPtr(), "2018 + DMu4_er2p0_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num2hist.GetPtr(), "2018 + DMu3_er2p0_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num3hist.GetPtr(), "2018 + DMu2_er2p0_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num4hist.GetPtr(), "2018 + DMu1_er2p0_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num5hist.GetPtr(), "2018 + DMu0_er2p0_SQ_OS_dR_Max1p6", "f")

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
ratio5.Draw("esame")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio2.Draw("esame")
ratio1.Draw("esame")
ratio0.Draw("esame")

pad2.SetLeftMargin(0.15)
#canvas1.Print("L1check_ptscan_er2p0_dR2p0.png")
canvas1.Print("L1check_ptscan_er2p0_dR1p6.png")
canvas1.Clear("")
canvas1.Close()


#mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
#mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
num0 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0")
#num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p2SQOSdR2p0 > 0")
#num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p4SQOSdR2p0 > 0")
#num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p6SQOSdR2p0 > 0")
#num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p8SQOSdR2p0 > 0")
#num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR2p0 > 0")

num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p2SQOSdR1p6 > 0")
num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p4SQOSdR1p6 > 0")
num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p6SQOSdR1p6 > 0")
num4 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er1p8SQOSdR1p6 > 0")
num5 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0")

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
#legend.AddEntry(num1hist.GetPtr(), "2018 + DMu3_er1p2_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num2hist.GetPtr(), "2018 + DMu3_er1p4_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num3hist.GetPtr(), "2018 + DMu3_er1p6_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num4hist.GetPtr(), "2018 + DMu3_er1p8_SQ_OS_dR_Max2p0", "f")
#legend.AddEntry(num5hist.GetPtr(), "2018 + DMu3_er2p0_SQ_OS_dR_Max2p0", "f")

legend.AddEntry(num1hist.GetPtr(), "2018 + DMu3_er1p2_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num2hist.GetPtr(), "2018 + DMu3_er1p4_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num3hist.GetPtr(), "2018 + DMu3_er1p6_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num4hist.GetPtr(), "2018 + DMu3_er1p8_SQ_OS_dR_Max1p6", "f")
legend.AddEntry(num5hist.GetPtr(), "2018 + DMu3_er2p0_SQ_OS_dR_Max1p6", "f")

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
ratio5.Draw("esame")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio2.Draw("esame")
ratio1.Draw("esame")
ratio0.Draw("esame")

pad2.SetLeftMargin(0.15)
#canvas1.Print("L1check_3_erScan_dR2p0.png")
canvas1.Print("L1check_3_erScan_dR1p6.png")
canvas1.Clear("")
canvas1.Close()

exit()
'''
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
#num0 = denom.Filter("L15er2p5SQOSdR1p6 > 0")
#num1 = denom.Filter("L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 > 0")
#num2 = denom.Filter("L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p5SQOSdE1p2 > 0")
#num2 = denom.Filter("L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0")

num0 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0")
num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0")
#num1 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L14er2p0SQOSdR1p6 > 0")

xmin = 0.
xmax = 8.5
nbins = int((xmax-xmin)/0.02)

denhist = denom.Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
num0hist = num0.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num1hist = num1.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
#num2hist = num2.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")

ratio0 = ROOT.TEfficiency(num0hist.GetPtr(), denhist.GetPtr())
ratio0.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio0.SetConfidenceLevel(0.68)

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

#ratio2 = ROOT.TEfficiency(num2hist.GetPtr(), denhist.GetPtr())
#ratio2.SetStatisticOption(ROOT.TEfficiency.kFCP)
#ratio2.SetConfidenceLevel(0.68)

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
#num2hist.SetFillColor(4)

legend = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
#legend.AddEntry(num0hist.GetPtr(), "Seed A", "f")
#legend.AddEntry(num1hist.GetPtr(), "Seed A+B", "f")
#legend.AddEntry(num2hist.GetPtr(), "Seed A+B+C2", "f")

legend.AddEntry(num0hist.GetPtr(), "Seeds 2018", "f")
legend.AddEntry(num1hist.GetPtr(), "Seeds 2018+D1", "f")
#legend.AddEntry(num1hist.GetPtr(), "Seeds 2018+D2", "f")

denhist.GetYaxis().SetRangeUser(1e4, 3e7)

denhist.Draw()
#num2hist.Draw("same")
num1hist.Draw("same")
num0hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";m(#mu#mu) [GeV];Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)

ratio0.SetLineColor(2)
ratio1.SetLineColor(3)
#ratio2.SetLineColor(4)

emptyHist.Draw("")
#ratio2.Draw("esame")
ratio1.Draw("esame")
ratio0.Draw("esame")

pad2.SetLeftMargin(0.15)
#canvas1.Print("L1check_Dima_backup.png")
canvas1.Print("L1check_Vince_main.png")
canvas1.Clear("")
canvas1.Close()
exit()
denomA = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1")
denomB = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass > 4.5 & J_mass < 8.5")
denomC = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 5.3")
denomD = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass > 4.8 & J_mass < 5.8")

ddA = denomA.Count().GetValue()
ddB = denomB.Count().GetValue()
ddC = denomC.Count().GetValue()
ddD = denomD.Count().GetValue()

print("Option", "Full", "[4.5,8.5]", "[0,5.3]", "[4.8,5.8]")

#seeds = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0"
seeds = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L14er2p0SQOSdR1p6 > 0"

numA = denomA.Filter(seeds).Count().GetValue()
numB = denomB.Filter(seeds).Count().GetValue()
numC = denomC.Filter(seeds).Count().GetValue()
numD = denomD.Filter(seeds).Count().GetValue()
print(seeds, numA/ddA, numB/ddB, numC/ddC, numD/ddD)

'''

'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")

denom = mumudata.Filter("HLTLowMassInclusive == 1").Filter("J_mass < 8.5")

num0 = denom.Filter(seeds0)
num1 = denom.Filter(seeds1)
num2 = denom.Filter(seeds2)
num3 = denom.Filter(seeds3)
num4 = denom.Filter(seeds4)
num5 = denom.Filter(seeds5)

#num0 = denom.Filter("L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0")
#num1 = denom.Filter("L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 + L13er2p0SQOSdR1p6 + L10er1p5SQOSdE1p2 > 0")
#num2 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0")
#num3 = denom.Filter("L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0")

xmin = 0.
xmax = 8.5
nbins = int((xmax-xmin)/0.02)

denhist = denom.Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
num0hist = num0.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num1hist = num1.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num2hist = num2.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num3hist = num3.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num4hist = num4.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")
num5hist = num5.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")

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

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num0hist.SetLineColor(2)
num1hist.SetLineColor(3)
num2hist.SetLineColor(4)
num3hist.SetLineColor(7)
num4hist.SetLineColor(6)
num5hist.SetLineColor(8)

legend = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num0hist.GetPtr(), "Dima Option 1", "l")
legend.AddEntry(num1hist.GetPtr(), "dEta -> DMu3er2p0SQOSdR2p0", "l")
legend.AddEntry(num2hist.GetPtr(), "dEta -> DMu0er1p4SQOSdR2p0", "l")
legend.AddEntry(num3hist.GetPtr(), "dEta -> DMu3er1p4SQOSdR2p0", "l")
legend.AddEntry(num4hist.GetPtr(), "dEta -> DMu0er1p4SQOSdR1p6", "l")
legend.AddEntry(num5hist.GetPtr(), "2018 seeds", "l")

legend = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num0hist.GetPtr(), "Opt 1", "l")
legend.AddEntry(num1hist.GetPtr(), "deltaR Opt 2", "l")
legend.AddEntry(num2hist.GetPtr(), "Opt 3 - main", "l")
legend.AddEntry(num3hist.GetPtr(), "Opt 3 - backup", "l")
legend.AddEntry(num4hist.GetPtr(), "2018 seeds", "l")

#legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
#legend.AddEntry(num0hist.GetPtr(), "Option 1", "l")
#legend.AddEntry(num1hist.GetPtr(), "Option 2", "l")
#legend.AddEntry(num2hist.GetPtr(), "Option 3", "l")
#legend.AddEntry(num3hist.GetPtr(), "2018 seeds", "l")

denhist.GetYaxis().SetRangeUser(1e4, 3e7)

denhist.Draw()
#num5hist.Draw("same")
num4hist.Draw("same")
num3hist.Draw("same")
num2hist.Draw("same")
num1hist.Draw("same")
num0hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";m(#mu#mu) [GeV];Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

ratio0.SetLineColor(2)
ratio1.SetLineColor(3)
ratio2.SetLineColor(4)
ratio3.SetLineColor(7)
ratio4.SetLineColor(6)
ratio5.SetLineColor(8)

emptyHist.Draw("")
#ratio5.Draw("esame")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio2.Draw("esame")
ratio1.Draw("esame")
ratio0.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("L1check_main_VS_summary.png")
canvas1.Clear("")
canvas1.Close()

'''
'''
#pt efficiency

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1").Filter("J_mass > 4.8 & J_mass < 5.8")#.Filter("abs(J_eta1) < 1.4 & abs(J_eta2) < 1.4")

seeds0 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0dEta1p5 > 0" #2022
seeds1 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0" #opt1
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0" #opt3
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0" #2018

denom = mumudata.Filter(seeds0)
num1 = denom.Filter(seeds1)
num3 = denom.Filter(seeds3)
num4 = denom.Filter(seeds4)

xmin = 0.
xmax = 30.
nbins = int((xmax-xmin)/0.5)

denhist = denom.Histo1D(("J_pt", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.5 GeV", nbins, xmin, xmax), "J_pt")
num1hist = num1.Histo1D(("J_pt", "dimuon pt", nbins, xmin, xmax), "J_pt")
num3hist = num3.Histo1D(("J_pt", "dimuon pt", nbins, xmin, xmax), "J_pt")
num4hist = num4.Histo1D(("J_pt", "dimuon pt", nbins, xmin, xmax), "J_pt")

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

ratio3 = ROOT.TEfficiency(num3hist.GetPtr(), denhist.GetPtr())
ratio3.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio3.SetConfidenceLevel(0.68)

ratio4 = ROOT.TEfficiency(num4hist.GetPtr(), denhist.GetPtr())
ratio4.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio4.SetConfidenceLevel(0.68)

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num1hist.SetLineColor(3)
num3hist.SetLineColor(7)
num4hist.SetLineColor(6)

legend = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num1hist.GetPtr(), "Option 1", "l")
legend.AddEntry(num3hist.GetPtr(), "Option 3", "l")
legend.AddEntry(num4hist.GetPtr(), "2018 L1 seeds", "l")

denhist.Draw()
num4hist.Draw("same")
num3hist.Draw("same")
num1hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";p_{T}(#mu#mu) [GeV];Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

ratio1.SetLineColor(3)
ratio3.SetLineColor(7)
ratio4.SetLineColor(6)

emptyHist.Draw("")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio1.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("ptEfficiency_4p8_5p8.png")
canvas1.Clear("")
canvas1.Close()
'''

#deltaEta efficiency
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1").Filter("J_mass > 4.8 & J_mass < 5.8")#.Filter("abs(J_eta1) < 1.4 & abs(J_eta2) < 1.4")

seeds0 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0dEta1p5 > 0" #2022
seeds1 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0" #opt1
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0" #opt3
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0" #2018

denom = mumudata.Filter(seeds0)
num1 = denom.Filter(seeds1)
num3 = denom.Filter(seeds3)
num4 = denom.Filter(seeds4)

xmin = 0.
xmax = 2.0
nbins = int((xmax-xmin)/0.05)

denhist = denom.Histo1D(("J_dEta",";#Delta#eta(#mu^{+}#mu^{-});Candidates / 0.05", nbins, xmin, xmax), "J_deltaEta")
num1hist = num1.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaEta")
num3hist = num3.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaEta")
num4hist = num4.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaEta")

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

ratio3 = ROOT.TEfficiency(num3hist.GetPtr(), denhist.GetPtr())
ratio3.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio3.SetConfidenceLevel(0.68)

ratio4 = ROOT.TEfficiency(num4hist.GetPtr(), denhist.GetPtr())
ratio4.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio4.SetConfidenceLevel(0.68)

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num1hist.SetLineColor(1)
num3hist.SetLineColor(2)
num4hist.SetLineColor(3)

legend = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num1hist.GetPtr(), "Option 1", "l")
legend.AddEntry(num3hist.GetPtr(), "Option 3", "l")
legend.AddEntry(num4hist.GetPtr(), "2018 L1 seeds", "l")

denhist.Draw()
num4hist.Draw("same")
num3hist.Draw("same")
num1hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";#Delta#eta(#mu#mu);Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

ratio1.SetLineColor(1)
ratio3.SetLineColor(2)
ratio4.SetLineColor(3)

emptyHist.Draw("")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio1.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("dEtaEfficiency_4p8_5p8.png")
canvas1.Clear("")
canvas1.Close()

#deltaEta efficiency

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1").Filter("J_mass > 3.0 & J_mass < 3.2")#.Filter("abs(J_eta1) < 1.4 & abs(J_eta2) < 1.4")

seeds0 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0dEta1p5 > 0" #2022
seeds1 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0" #opt1
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0" #opt3
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0" #2018

denom = mumudata.Filter(seeds0)
num1 = denom.Filter(seeds1)
num3 = denom.Filter(seeds3)
num4 = denom.Filter(seeds4)

xmin = 0.
xmax = 1.5
nbins = int((xmax-xmin)/0.05)

denhist = denom.Histo1D(("J_dEta",";#Delta#eta(#mu^{+}#mu^{-});Candidates / 0.05", nbins, xmin, xmax), "J_deltaEta")
num1hist = num1.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaEta")
num3hist = num3.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaEta")
num4hist = num4.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaEta")

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

ratio3 = ROOT.TEfficiency(num3hist.GetPtr(), denhist.GetPtr())
ratio3.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio3.SetConfidenceLevel(0.68)

ratio4 = ROOT.TEfficiency(num4hist.GetPtr(), denhist.GetPtr())
ratio4.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio4.SetConfidenceLevel(0.68)

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num1hist.SetLineColor(1)
num3hist.SetLineColor(2)
num4hist.SetLineColor(3)

legend = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num1hist.GetPtr(), "Option 1", "l")
legend.AddEntry(num3hist.GetPtr(), "Option 3", "l")
legend.AddEntry(num4hist.GetPtr(), "2018 L1 seeds", "l")

denhist.Draw()
num4hist.Draw("same")
num3hist.Draw("same")
num1hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";#Delta#eta(#mu#mu);Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

ratio1.SetLineColor(1)
ratio3.SetLineColor(2)
ratio4.SetLineColor(3)

emptyHist.Draw("")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio1.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("dEtaEfficiency_3p0_3p2.png")
canvas1.Clear("")
canvas1.Close()

#deltaR efficiency

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1").Filter("J_mass > 3.0 & J_mass < 3.2")#.Filter("abs(J_eta1) < 1.4 & abs(J_eta2) < 1.4")

seeds0 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0dEta1p5 > 0" #2022
seeds1 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0" #opt1
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0" #opt3
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0" #2018

denom = mumudata.Filter(seeds0)
num1 = denom.Filter(seeds1)
num3 = denom.Filter(seeds3)
num4 = denom.Filter(seeds4)

xmin = 0.
xmax = 1.5
nbins = int((xmax-xmin)/0.05)

denhist = denom.Histo1D(("J_dEta",";#DeltaR(#mu^{+}#mu^{-});Candidates / 0.05", nbins, xmin, xmax), "J_deltaR")
num1hist = num1.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaR")
num3hist = num3.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaR")
num4hist = num4.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaR")

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

ratio3 = ROOT.TEfficiency(num3hist.GetPtr(), denhist.GetPtr())
ratio3.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio3.SetConfidenceLevel(0.68)

ratio4 = ROOT.TEfficiency(num4hist.GetPtr(), denhist.GetPtr())
ratio4.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio4.SetConfidenceLevel(0.68)

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num1hist.SetLineColor(1)
num3hist.SetLineColor(2)
num4hist.SetLineColor(3)

legend = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num1hist.GetPtr(), "Option 1", "l")
legend.AddEntry(num3hist.GetPtr(), "Option 3", "l")
legend.AddEntry(num4hist.GetPtr(), "2018 L1 seeds", "l")

denhist.Draw()
num4hist.Draw("same")
num3hist.Draw("same")
num1hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";#DeltaR(#mu#mu);Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

ratio1.SetLineColor(1)
ratio3.SetLineColor(2)
ratio4.SetLineColor(3)

emptyHist.Draw("")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio1.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("dREfficiency_3p0_3p2.png")
canvas1.Clear("")
canvas1.Close()

#deltaR efficiency

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1").Filter("J_mass > 4.8 & J_mass < 5.8")#.Filter("abs(J_eta1) < 1.4 & abs(J_eta2) < 1.4")

seeds0 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L10er2p0dEta1p5 > 0" #2022
seeds1 = "L15er2p5SQOSdR1p6 + L14er2p0SQOSdR1p6 + L10er1p4SQOSdE1p2 > 0" #opt1
seeds3 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 + L13er2p0SQOSdR1p6 > 0" #opt3
seeds4 = "L10er1p5SQOSdR1p4 + L14er2p5SQOSdR1p2 > 0" #2018

denom = mumudata.Filter(seeds0)
num1 = denom.Filter(seeds1)
num3 = denom.Filter(seeds3)
num4 = denom.Filter(seeds4)

xmin = 0.
xmax = 2.5
nbins = int((xmax-xmin)/0.05)

denhist = denom.Histo1D(("J_dEta",";#DeltaR(#mu^{+}#mu^{-});Candidates / 0.1 GeV", nbins, xmin, xmax), "J_deltaR")
num1hist = num1.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaR")
num3hist = num3.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaR")
num4hist = num4.Histo1D(("J_dEta", "dimuon pt", nbins, xmin, xmax), "J_deltaR")

ratio1 = ROOT.TEfficiency(num1hist.GetPtr(), denhist.GetPtr())
ratio1.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio1.SetConfidenceLevel(0.68)

ratio3 = ROOT.TEfficiency(num3hist.GetPtr(), denhist.GetPtr())
ratio3.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio3.SetConfidenceLevel(0.68)

ratio4 = ROOT.TEfficiency(num4hist.GetPtr(), denhist.GetPtr())
ratio4.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio4.SetConfidenceLevel(0.68)

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
num1hist.SetLineColor(1)
num3hist.SetLineColor(2)
num4hist.SetLineColor(3)

legend = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(),  "2022 L1 seeds", "l")
legend.AddEntry(num1hist.GetPtr(), "Option 1", "l")
legend.AddEntry(num3hist.GetPtr(), "Option 3", "l")
legend.AddEntry(num4hist.GetPtr(), "2018 L1 seeds", "l")

denhist.Draw()
num4hist.Draw("same")
num3hist.Draw("same")
num1hist.Draw("same")

legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";#DeltaR(#mu#mu);Ratio option/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

ratio1.SetLineColor(1)
ratio3.SetLineColor(2)
ratio4.SetLineColor(3)

emptyHist.Draw("")
ratio4.Draw("esame")
ratio3.Draw("esame")
ratio1.Draw("esame")

pad2.SetLeftMargin(0.15)
canvas1.Print("dREfficiency_4p8_5p8.png")
canvas1.Clear("")
canvas1.Close()

#2D hist
'''
mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1")
'''
#

histmassdr = mumudata.Histo2D(("", ";m(#mu#mu)[GeV];#DeltaR", 200, 0, 8.5, 200, 0, 3.0), "J_mass", "J_deltaR")

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histmassdr.SetStats(0)
histmassdr.Draw("colz")

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.Print("mass_dR.png")
c1.Clear()
c1.Close()

#

histmassdeta = mumudata.Histo2D(("", ";m(#mu#mu)[GeV];#Delta#eta", 200, 0, 8.5, 200, 0, 2.0), "J_mass", "J_deltaEta")

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histmassdeta.SetStats(0)
histmassdeta.Draw("colz")

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.Print("mass_dEta.png")
c1.Clear()
c1.Close()

#

histptdetaJ = mumudata.Filter("J_mass > 3.0 & J_mass < 3.2").Histo2D(("", ";#Delta #eta;p_{T}(#mu#mu) [GeV]", 200, 0, 3.0, 200, 5., 50), "J_deltaEta", "J_pt")

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histptdetaJ.SetStats(0)
histptdetaJ.Draw("colz")

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.Print("dEta_pt_3p0_3p2.png")
c1.Clear()
c1.Close()

#
'''
histptdrJ = mumudata.Filter("J_mass > 3.0 & J_mass < 3.2").Histo2D(("", ";#DeltaR;p_{T}(#mu#mu) [GeV]", 200, 0, 1.5, 200, 0., 50), "J_deltaR", "J_pt")

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histptdrJ.SetStats(0)
histptdrJ.Draw("colz")

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.Print("dR_pt_3p0_3p2.png")
c1.Clear()
c1.Close()
#
exit()
histptdetaB = mumudata.Filter("J_mass > 4.8 & J_mass < 5.8").Histo2D(("", ";#Delta #eta;p_{T}(#mu#mu) [GeV]", 200, 0, 3.0, 200, 5., 50), "J_deltaEta", "J_pt")
'''
c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histptdetaB.SetStats(0)
histptdetaB.Draw("colz")

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.Print("dEta_pt_4p8_5p8.png")
c1.Clear()
c1.Close()

#

histptdrB = mumudata.Filter("J_mass > 4.8 & J_mass < 5.8").Histo2D(("", ";#DeltaR;p_{T}(#mu#mu) [GeV]", 200, 0, 3.0, 200, 5., 50), "J_deltaR", "J_pt")

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histptdrB.SetStats(0)
histptdrB.Draw("colz")

c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.Print("dR_pt_4p8_5p8.png")
c1.Clear()
c1.Close()
'''

#Check of dima's plot

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass > 3.0 & J_mass < 3.2")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 3.0 & J_Prob > 0.05")
mumudata = mumudata.Filter("HLTLowMassInclusive == 1")

denom = mumudata.Filter("L10er2p0dEta1p5 > 0").Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5")
num = denom.Filter("L10er1p5dR > 0")

xmin = 0.
xmax = 2.0
nbins = int((xmax-xmin)/0.05)

denhist = denom.Histo1D(("dR_reco", ";#DeltaR^{off}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.05", nbins, xmin, xmax), "J_deltaR")
numhist = num.Histo1D(("dR_reco", "dimuon mass", nbins, xmin, xmax), "J_deltaR")

ratio = ROOT.TEfficiency(numhist.GetPtr(), denhist.GetPtr())
ratio.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio.SetConfidenceLevel(0.68)

emptyHist = ROOT.TH1D("emptyHist", ";#DeltaR(#mu#mu);Ratio", nbins, xmin, xmax)
emptyHist.GetXaxis().SetTitle(denhist.GetXaxis().GetTitle())
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

emptyHist.Draw("")
ratio.Draw("same")

c1.Print("dReff_dR1p4_dE2p0.png")
c1.Clear()
c1.Close()

xmin = 0.
xmax = 2.0
nbins = int((xmax-xmin)/0.05)

denhist = denom.Histo1D(("dR_reco", ";#DeltaR^{off}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.05", nbins, xmin, xmax), "J_deltaR_L1")
numhist = num.Histo1D(("dR_reco", "dimuon mass", nbins, xmin, xmax), "J_deltaR_L1")

ratio = ROOT.TEfficiency(numhist.GetPtr(), denhist.GetPtr())
ratio.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio.SetConfidenceLevel(0.68)

emptyHist = ROOT.TH1D("emptyHist", ";#DeltaR(#mu#mu);Ratio", nbins, xmin, xmax)
emptyHist.GetXaxis().SetTitle(denhist.GetXaxis().GetTitle())
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)
emptyHist.GetYaxis().SetLabelSize(0.06)

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

emptyHist.Draw("")
ratio.Draw("same")

c1.Print("dRL1eff_dR1p4_dE2p0.png")
c1.Clear()
c1.Close()

exit()




























'''
mumudata = data.Filter("J_mass > 0. & J_mass < 2.9")
mumudata = mumudata.Filter("J_pt1 > 2.0 & J_pt2 > 2.0 & J_Prob > 0.05")

histpteta = mumudata.Histo2D(("", "2.9 < m < 3.3, prob > 5%;#eta(#mu);p_{T}(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta1", "J_pt1")
histpteta2 = mumudata.Histo2D(("", ";#eta(#mu);p_T(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta2", "J_pt2")

histpteta.Add(histpteta2.GetPtr())

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histpteta.SetStats(0)
histpteta.Draw("colz")

c1.Print("pt_eta_0_2p9.png")
c1.Clear()
c1.Close()

mumudata = data.Filter("J_mass > 2.9 & J_mass < 3.3")
mumudata = mumudata.Filter("J_pt1 > 2.0 & J_pt2 > 2.0 & J_Prob > 0.05")

histpteta = mumudata.Histo2D(("", "2.9 < m < 3.3, prob > 5%;#eta(#mu);p_{T}(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta1", "J_pt1")
histpteta2 = mumudata.Histo2D(("", ";#eta(#mu);p_T(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta2", "J_pt2")

histpteta.Add(histpteta2.GetPtr())

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histpteta.SetStats(0)
histpteta.Draw("colz")

c1.Print("pt_eta_2p9_3p3.png")
c1.Clear()
c1.Close()

mumudata = data.Filter("J_mass > 0 & J_mass < 4.0")
mumudata = mumudata.Filter("J_pt1 > 2.0 & J_pt2 > 2.0 & J_Prob > 0.05")

histpteta = mumudata.Histo2D(("", "0. < m < 4.0, prob > 5%;#eta(#mu);p_{T}(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta1", "J_pt1")
histpteta2 = mumudata.Histo2D(("", ";#eta(#mu);p_T(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta2", "J_pt2")

histpteta.Add(histpteta2.GetPtr())

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histpteta.SetStats(0)
histpteta.Draw("colz")

c1.Print("pt_eta_0_4p0.png")
c1.Clear()
c1.Close()

mumudata = data.Filter("J_mass > 4.0 & J_mass < 8.5")
mumudata = mumudata.Filter("J_pt1 > 2.0 & J_pt2 > 2.0 & J_Prob > 0.05")

histpteta = mumudata.Histo2D(("", "4.0 < m < 8.5, prob > 5%;#eta(#mu);p_{T}(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta1", "J_pt1")
histpteta2 = mumudata.Histo2D(("", ";#eta(#mu);p_T(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta2", "J_pt2")

histpteta.Add(histpteta2.GetPtr())

c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histpteta.SetStats(0)
histpteta.Draw("colz")

c1.Print("pt_eta_4p0_8p5.png")
c1.Clear()
c1.Close()
'''

mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 4.0 & J_pt2 > 4.0 & J_Prob > 0.1") \
.Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5") 

denom = mumudata.Filter("HLTLowMassDisplaced == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5")
num = denom.Filter("L10er1p5dR + L14dR > 0")

xmin = 0.
xmax = 8.5
nbins = int((xmax-xmin)/0.2)

denhist = denom.Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
numhist = num.Histo1D(("J_mass", "dimuon mass", nbins, xmin, xmax), "J_mass")

ratio = ROOT.TEfficiency(numhist.GetPtr(), denhist.GetPtr())
ratio.SetStatisticOption(ROOT.TEfficiency.kFCP)
ratio.SetConfidenceLevel(0.68)


# In[38]:


canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)

pad1.Draw("")
pad2.Draw("")

pad1.cd()
pad1.SetLogy()

denhist.SetStats(0)
numhist.SetFillColor(3)

legend = ROOT.TLegend(0.7, 0.8, 0.89, 0.89)
legend.AddEntry(denhist.GetPtr(), "All L1 seeds", "l")
legend.AddEntry(numhist.GetPtr(), "Only 2018 seeds", "f")

#denhist.GetYaxis().SetRangeUser(0, 4e5)

denhist.Draw()
numhist.Draw("same")
legend.Draw("same")
pad1.SetLeftMargin(0.15)

pad2.cd()

emptyHist = ROOT.TH1D("emptyHist", ";m(#mu#mu) [GeV];Ratio 2018/2022", nbins, xmin, xmax)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)

emptyHist.Draw("")
ratio.Draw("esame")
pad2.SetLeftMargin(0.15)

canvas1.Draw()


# In[39]:


mumudata = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")
mumudata = mumudata.Filter("J_pt1 > 2.0 & J_pt2 > 2.0 & J_Prob > 0.01")
#Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5")

diff = mumudata.Filter("HLTLowMassInclusive == 1").Filter("L10er2p0dEta1p5 == 1").Filter("J_mass < 8.5").Filter("L10er1p5dR + L14dR == 0")

xmin = 2.9
xmax = 3.3
nbins = int((xmax-xmin)/0.02)

diff = diff.Filter("J_mass > " + str(xmin) + " & J_mass < " + str(xmax))
diff = diff.Filter("J_pt1 > 4. & J_pt2 > 3 & J_Prob > 0.1")
diff = diff.Filter("J_deltaEta < 1.6")

masshist = diff.Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
massBarHist0 = diff.Filter("abs(J_eta1) < 2.0 & abs(J_eta2) < 2.0").Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
massBarHist1 = diff.Filter("abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5").Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
massBarHist2 = diff.Filter("abs(J_eta1) < 1.4 & abs(J_eta2) < 1.4").Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")
massBarHist3 = diff.Filter("abs(J_eta1) < 1.2 & abs(J_eta2) < 1.2").Histo1D(("J_mass", ";m(#mu^{+}#mu^{-}) [GeV];Candidates / 20 MeV", nbins, xmin, xmax), "J_mass")


# In[40]:


c0 = ROOT.TCanvas("c0", "c0", 1200, 1200)

masshist.SetStats(0)
massBarHist0.SetFillColor(2)
massBarHist1.SetFillColor(3)
massBarHist2.SetFillColor(4)
massBarHist3.SetFillColor(5)

masshist.GetYaxis().SetRangeUser(0, 50000)

masshist.Draw("")
massBarHist0.Draw("same")
massBarHist1.Draw("same")
massBarHist2.Draw("same")
massBarHist3.Draw("same")

c0.Draw("")


# In[29]:


histpteta = diff.Histo2D(("", ";#eta(#mu);p_{T}(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta1", "J_pt1")
histpteta2 = diff.Histo2D(("", ";#eta(#mu);p_T(#mu) [GeV]", 200, -2.5, 2.5, 200, 0., 10), "J_eta2", "J_pt2")

histpteta.Add(histpteta2.GetPtr())


# In[30]:


c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histpteta.SetStats(0)
histpteta.Draw("colz")

c1.Draw("")


# In[185]:


histptdeta = diff.Histo2D(("", ";#Delta #eta;p_{T}(#mu#mu) [GeV]", 200, 0, 3.0, 200, 0., 10), "J_deltaEta", "J_pt")


# In[186]:


c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histptdeta.SetStats(0)
histptdeta.Draw("colz")

c1.Draw("")


# In[187]:


histptdr = diff.Histo2D(("", ";#Delta R;p_{T}(#mu#mu) [GeV]", 200, 0, 3.0, 200, 0., 10), "J_deltaR", "J_pt")


# In[188]:


c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histptdr.SetStats(0)
histptdr.Draw("colz")

c1.Draw("")


# In[130]:


histetadeta = diff.Histo2D(("", ";#eta(#mu#mu);#Delta #eta", 200, -2.5, 2.5, 200, 0., 3), "J_eta1", "J_deltaR")


# In[131]:


c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histetadeta.SetStats(0)
histetadeta.Draw("colz")

c1.Draw("")


# In[132]:


histeta1deta = diff.Histo2D(("", ";#eta(#mu);#Delta #eta", 200, -2.5, 2.5, 200, 0., 2), "J_eta1", "J_deltaR")
histeta2deta = diff.Histo2D(("", ";#eta(#mu);#Delta #eta", 200, -2.5, 2.5, 200, 0., 2), "J_eta2", "J_deltaR")
histeta1deta.Add(histeta2deta.GetPtr())


# In[133]:


c1 = ROOT.TCanvas("c1", "c1", 1200, 1200)

histeta1deta.SetStats(0)
histeta1deta.Draw("colz")

c1.Draw("")


# In[ ]:





# In[ ]:





# In[ ]:





# In[44]:


histptdr = denom.Histo2D(("", ";m(#mu^{+}#mu^{-}) [GeV];#Delta R", 400, 0., 8.5, 100, 0., 2.5), "J_mass", "J_deltaR")


# In[32]:


denpthist = denom.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 400, 0., 80), "J_pt")
numpthist = num.Histo1D(("", "dimuon mass", 400, 0., 80), "J_pt")


# In[38]:


denpthist.SetStats(0)

numpthist.SetFillColor(3)

string = ROOT.TLatex()
string.SetTextSize(0.04)

legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.85)
legend.AddEntry(denpthist.GetPtr(), "All L1 seeds", "l")
legend.AddEntry(numpthist.GetPtr(), "Only 2018 seeds", "f")

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)
canvas1.SetLogy()

denpthist.GetYaxis().SetRangeUser(1, 1e7)

denpthist.Draw()
numpthist.Draw("same")

#string.DrawLatex(8.6, 5.2e4, "DoubleMuonLowMass - Run 355558-355559")
#string.DrawLatex(8.5, 6.2e4, "ParkingDoubleMuonLowMass0 - Run 355872")
string.DrawLatex(5.0, 3.2e8, "ParkingDoubleMuonLowMassX - Run2022G")

legend.Draw("same")

canvas1.SetLeftMargin(0.15)
canvas1.Draw()


# In[46]:


denetahist = denom.Histo1D(("", ";#eta(#mu^{+}#mu^{-}) [GeV];Candidates / 0.01", 500, -2.5, 2.5), "J_eta")
numetahist = num.Histo1D(("", "dimuon mass", 500, -2.5, 2.5), "J_eta")


# In[47]:


denetahist.SetStats(0)

numetahist.SetFillColor(3)

string = ROOT.TLatex()
string.SetTextSize(0.04)

legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.85)
legend.AddEntry(denetahist.GetPtr(), "All L1 seeds", "l")
legend.AddEntry(numetahist.GetPtr(), "Only 2018 seeds", "f")

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)
#canvas1.SetLogy()

denetahist.GetYaxis().SetRangeUser(1, 1e6)

denetahist.Draw()
numetahist.Draw("same")

#string.DrawLatex(8.6, 5.2e4, "DoubleMuonLowMass - Run 355558-355559")
#string.DrawLatex(8.5, 6.2e4, "ParkingDoubleMuonLowMass0 - Run 355872")
string.DrawLatex(5.0, 3.2e8, "ParkingDoubleMuonLowMassX - Run2022G")

legend.Draw("same")

canvas1.SetLeftMargin(0.15)
canvas1.Draw()


# In[51]:


denpt1hist = denom.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 400, 0., 10), "J_pt1")
numpt1hist = num.Histo1D(("", "dimuon mass", 400, 0., 10), "J_pt1")

denpt2hist = denom.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 400, 0., 10), "J_pt2")
numpt2hist = num.Histo1D(("", "dimuon mass", 400, 0., 10), "J_pt2")

denpt1hist.Add(denpt2hist.GetPtr())
numpt1hist.Add(numpt2hist.GetPtr())


# In[52]:


denpt1hist.SetStats(0)

numpt1hist.SetFillColor(3)

string = ROOT.TLatex()
string.SetTextSize(0.04)

legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.85)
legend.AddEntry(denpt1hist.GetPtr(), "All L1 seeds", "l")
legend.AddEntry(numpt1hist.GetPtr(), "Only 2018 seeds", "f")

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)
canvas1.SetLogy()

denpt1hist.GetYaxis().SetRangeUser(1, 1e8)

denpt1hist.Draw()
numpt1hist.Draw("same")

#string.DrawLatex(8.6, 5.2e4, "DoubleMuonLowMass - Run 355558-355559")
#string.DrawLatex(8.5, 6.2e4, "ParkingDoubleMuonLowMass0 - Run 355872")
string.DrawLatex(5.0, 3.2e8, "ParkingDoubleMuonLowMassX - Run2022G")

legend.Draw("same")

canvas1.SetLeftMargin(0.15)
canvas1.Draw()


# In[55]:


deneta1hist = denom.Histo1D(("", ";#eta(#mu) [GeV];Candidates / 0.01", 500, -2.5, 2.5), "J_eta1")
numeta1hist = num.Histo1D(("", "dimuon mass", 500, -2.5, 2.5), "J_eta1")

deneta2hist = denom.Histo1D(("", ";#eta(#mu) [GeV];Candidates / 0.01", 500, -2.5, 2.5), "J_eta2")
numeta2hist = num.Histo1D(("", "dimuon mass", 500, -2.5, 2.5), "J_eta2")

deneta1hist.Add(deneta2hist.GetPtr())
numeta1hist.Add(numeta2hist.GetPtr())


# In[56]:


deneta1hist.SetStats(0)

numeta1hist.SetFillColor(3)

string = ROOT.TLatex()
string.SetTextSize(0.04)

legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.85)
legend.AddEntry(deneta1hist.GetPtr(), "All L1 seeds", "l")
legend.AddEntry(numeta1hist.GetPtr(), "Only 2018 seeds", "f")

canvas1 = ROOT.TCanvas()
canvas1.SetCanvasSize(1200, 1200)
canvas1.SetLogy()

deneta1hist.GetYaxis().SetRangeUser(1, 1e8)

deneta1hist.Draw()
numeta1hist.Draw("same")

#string.DrawLatex(8.6, 5.2e4, "DoubleMuonLowMass - Run 355558-355559")
#string.DrawLatex(8.5, 6.2e4, "ParkingDoubleMuonLowMass0 - Run 355872")
string.DrawLatex(5.0, 3.2e8, "ParkingDoubleMuonLowMassX - Run2022G")

legend.Draw("same")

canvas1.SetLeftMargin(0.15)
canvas1.Draw()


# In[258]:


#canvas1.Print("MuMu_MuMuTrigger_2022G.pdf")


# In[49]:


#efficiency "tnp"

data0 = data.Filter("HLTLowMassInclusive == 1 & J_mass < 8.5 & J_Prob > 0.01 & J_pt1 > 2. & J_pt2 > 2 & \
abs(J_eta1) < 1.5 & abs(J_eta2) < 1.5")

numtnp = data0.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1")
num2tnp = numtnp.Filter("L10er1p5dR + L14dR > 0")
dentnp = data0.Filter("mu1HLTmatched + mu2HLTmatched > 0")


# In[50]:


pt1TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu) [GeV];Candidates / 0.2 GeV", 50, 0., 20), "J_pt1")
pt2TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu) [GeV];Candidates / 0.2 GeV", 50, 0., 20), "J_pt2")

pt1TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu) [GeV];Candidates / 0.2 GeV", 50, 0., 20), "J_pt1")
pt2TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu) [GeV];Candidates / 0.2 GeV", 50, 0., 20), "J_pt2")

pt1TnpDenHist = dentnp.Filter("mu1HLTmatched == 1").Histo1D(("", ";p_{T}(#mu) [GeV];Candidates / 0.2 GeV", 50, 0., 20), "J_pt1")
pt2TnpDenHist = dentnp.Filter("mu2HLTmatched == 1").Histo1D(("", ";p_{T}(#mu) [GeV];Candidates / 0.2 GeV", 50, 0., 20), "J_pt2")

pt1TnpNumHist.Add(pt2TnpNumHist.GetPtr())
pt1TnpNumHist2.Add(pt2TnpNumHist2.GetPtr())
pt1TnpDenHist.Add(pt2TnpDenHist.GetPtr())

effPtHist = ROOT.TEfficiency(pt1TnpNumHist.GetPtr(), pt1TnpDenHist.GetPtr())
effPtHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effPtHist.SetConfidenceLevel(0.68)

effPtHist2 = ROOT.TEfficiency(pt1TnpNumHist2.GetPtr(), pt1TnpDenHist.GetPtr())
effPtHist2.SetStatisticOption(ROOT.TEfficiency.kFCP)
effPtHist2.SetConfidenceLevel(0.68)


# In[51]:


c0 = ROOT.TCanvas()

emptyHist = ROOT.TH1D("emptyHist", ";p_{T}(#mu);Efficiency", 10, 0, 20)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)

legend = ROOT.TLegend(0.65, 0.2, 0.89, 0.4)
legend.AddEntry(effPtHist, "All L1 seeds", "l")
legend.AddEntry(effPtHist2, "Only 2018 seeds", "l")

line = ROOT.TLine(4., 0, 4., 1.05)
line.SetLineColor(ROOT.kBlue)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)

effPtHist2.SetLineColor(2)

emptyHist.Draw("")
effPtHist.Draw("esame")
effPtHist2.Draw("esame")
line.Draw("same")
legend.Draw("same")

c0.Draw("")


# In[52]:


eta1TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, -2.5, 2.5), "J_eta1")
eta2TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, -2.5, 2.5), "J_eta2")

eta1TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, -2.5, 2.5), "J_eta1")
eta2TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, -2.5, 2.5), "J_eta2")

eta1TnpDenHist = dentnp.Filter("mu1HLTmatched == 1").Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, -2.5, 2.5), "J_eta1")
eta2TnpDenHist = dentnp.Filter("mu2HLTmatched == 1").Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, -2.5, 2.5), "J_eta2")

eta1TnpNumHist.Add(eta2TnpNumHist.GetPtr())
eta1TnpNumHist2.Add(eta2TnpNumHist2.GetPtr())
eta1TnpDenHist.Add(eta2TnpDenHist.GetPtr())

effEtaHist = ROOT.TEfficiency(eta1TnpNumHist.GetPtr(), eta1TnpDenHist.GetPtr())
effEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effEtaHist.SetConfidenceLevel(0.68)

effEtaHist2 = ROOT.TEfficiency(eta1TnpNumHist2.GetPtr(), eta1TnpDenHist.GetPtr())
effEtaHist2.SetStatisticOption(ROOT.TEfficiency.kFCP)
effEtaHist2.SetConfidenceLevel(0.68)


# In[53]:


c0 = ROOT.TCanvas()

emptyHist = ROOT.TH1D("emptyHist", ";#eta(#mu);Efficiency", 10, -2.5, 2.5)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)

legend = ROOT.TLegend(0.65, 0.2, 0.89, 0.4)
legend.AddEntry(effEtaHist, "All L1 seeds", "l")
legend.AddEntry(effEtaHist2, "Only 2018 seeds", "l")

effEtaHist2.SetLineColor(2)
emptyHist.Draw("")
effEtaHist.Draw("same")
effEtaHist2.Draw("same")
legend.Draw("same")

c0.Draw("")


# In[80]:


dR1TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")
dR2TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")

dR1TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")
dR2TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")

dR1TnpDenHist = dentnp.Filter("mu1HLTmatched == 1").Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")
dR2TnpDenHist = dentnp.Filter("mu2HLTmatched == 1").Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")

dR1TnpNumHist.Add(dR2TnpNumHist.GetPtr())
dR1TnpNumHist2.Add(dR2TnpNumHist2.GetPtr())
dR1TnpDenHist.Add(dR2TnpDenHist.GetPtr())

effDrHist = ROOT.TEfficiency(dR1TnpNumHist.GetPtr(), dR1TnpDenHist.GetPtr())
effDrHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effDrHist.SetConfidenceLevel(0.68)

effDrHist2 = ROOT.TEfficiency(dR1TnpNumHist2.GetPtr(), dR1TnpDenHist.GetPtr())
effDrHist2.SetStatisticOption(ROOT.TEfficiency.kFCP)
effDrHist2.SetConfidenceLevel(0.68)


# In[81]:


c0 = ROOT.TCanvas()

emptyHist = ROOT.TH1D("emptyHist", ";#eta(#mu);Efficiency", 10, 0., 2.5)
emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyHist.SetStats(0)

legend = ROOT.TLegend(0.15, 0.2, 0.39, 0.4)
legend.AddEntry(effEtaHist, "All L1 seeds", "l")
legend.AddEntry(effEtaHist2, "Only 2018 seeds", "l")

effDrHist2.SetLineColor(2)
emptyHist.Draw("")
effDrHist.Draw("same")
effDrHist2.Draw("same")
legend.Draw("same")

c0.Draw("")


# In[85]:


hist = numtnp.Histo1D(("", "", 40, 0, 80), "nVtx")
#dR2TnpNumHist = numtnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")

#dR1TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")
#dR2TnpNumHist2 = num2tnp.Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")

#dR1TnpDenHist = dentnp.Filter("mu1HLTmatched == 1").Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")
#dR2TnpDenHist = dentnp.Filter("mu2HLTmatched == 1").Histo1D(("", ";p_{T}(#mu^{+}#mu^{-}) [GeV];Candidates / 0.2 GeV", 50, 0, 2.5), "J_deltaR")

#dR1TnpNumHist.Add(dR2TnpNumHist.GetPtr())
#dR1TnpNumHist2.Add(dR2TnpNumHist2.GetPtr())
#dR1TnpDenHist.Add(dR2TnpDenHist.GetPtr())

#effDrHist = ROOT.TEfficiency(dR1TnpNumHist.GetPtr(), dR1TnpDenHist.GetPtr())
#effDrHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
#effDrHist.SetConfidenceLevel(0.68)

#effDrHist2 = ROOT.TEfficiency(dR1TnpNumHist2.GetPtr(), dR1TnpDenHist.GetPtr())
#effDrHist2.SetStatisticOption(ROOT.TEfficiency.kFCP)
#effDrHist2.SetConfidenceLevel(0.68)


# In[86]:


c0 = ROOT.TCanvas()

hist.Draw("")

c0.Draw("")


# In[33]:


mumudata2 = mumudata.Filter("J_pt1 > 4. & J_pt2 > 4. & J_pt > 8 & J_Prob > 0.1 & abs(J_eta) < 2.4")


# In[ ]:





# In[233]:


dfcands = pd.DataFrame(mumudata2.Filter("J_mass > 2.9 & J_mass < 3.3").AsNumpy(columns=["J_mass"]))


# In[234]:


len(dfcands)


# In[98]:


jpsimass = ROOT.RooRealVar("J_mass", "m(#mu#mu) [GeV]", 2.9, 3.3)
mumuroodata = ROOT.RooDataSet.from_numpy({"J_mass": dfcands["J_mass"]}, [jpsimass])
mumuroohist = mumuroodata.binnedClone()


# In[99]:


mumuroohist.Print()


# In[100]:


jpsiframe = jpsimass.frame(Title="Dimuon Candidate Mass")


# In[101]:


mean = ROOT.RooRealVar("#mu", "mean of gaussians", 3.095, 3.09, 3.10) #mean value, min value, max value
sigma1 = ROOT.RooRealVar("#sigma_{barrel}", "width1", 0.02, 0.01, 0.05) # most probably this are not that good!
sigma2 = ROOT.RooRealVar("#sigma_{endcap}", "width2", 0.04, 0.02, 0.1)

sig1 = ROOT.RooGaussian("signal1", "signal1", jpsimass, mean, sigma1)
sig2 = ROOT.RooGaussian("signal2", "signal2", jpsimass, mean, sigma2)

a0 = ROOT.RooRealVar("a0", "a0", -0.5, -1, 1)
a1 = ROOT.RooRealVar("a1", "a1", -0.1, -1, +1)
bkg = ROOT.RooChebychev("bkg", "Background", jpsimass, [a0]) ## to increase the degree, just increase the coefficients

bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.1, 0.01, 0.2)
barfrac = ROOT.RooRealVar("f_{barrel}", "fraction of higher res", 0.4, 0.2, 0.6)
model = ROOT.RooAddPdf("model", "g+cheb", [bkg, sig1, sig2], [bkgfrac, barfrac])
#model = ROOT.RooAddPdf("model", "g+g", [sig1, sig2], [barfrac])


# In[102]:


r_full = model.fitTo(mumuroohist, Save=True)


# In[103]:


c0 = ROOT.TCanvas("canvas0", "canvas0", 1200, 600)

mumuroohist.plotOn(jpsiframe)
model.plotOn(jpsiframe) # By default only fitted range is shown
model.paramOn(jpsiframe, ROOT.RooFit.Parameters([mean, barfrac, bkgfrac, sigma1, sigma2]), ROOT.RooFit.Layout(0.6, 0.9, 0.9))

jpsiframe.Draw()
c0.Draw()


# In[104]:


num = mumuroodata.numEntries()*(1 - bkgfrac.getVal())
print(num)


# In[105]:


print("nevt x fb:", num/13.4)


# In[28]:


c0.Print("Jpsi_fit_Run2022F.png")


# In[106]:


dfcands = pd.DataFrame(mumudata2.Filter("J_mass > 3.5 & J_mass < 3.9").AsNumpy(columns=["J_mass"]))


# In[107]:


len(dfcands)


# In[108]:


jpsimass = ROOT.RooRealVar("J_mass", "m(#mu#mu) [GeV]", 3.5, 3.9)
mumuroodata = ROOT.RooDataSet.from_numpy({"J_mass": dfcands["J_mass"]}, [jpsimass])
mumuroohist = mumuroodata.binnedClone()


# In[109]:


jpsiframe = jpsimass.frame(Title="Dimuon Candidate Mass")


# In[110]:


mean = ROOT.RooRealVar("#mu", "mean of gaussians", 3.67, 3.6, 3.75) #mean value, min value, max value
sigma1 = ROOT.RooRealVar("#sigma_{barrel}", "width1", 0.02, 0.01, 0.05) # most probably this are not that good!
sigma2 = ROOT.RooRealVar("#sigma_{endcap}", "width2", 0.045, 0.02, 0.1)

sig1 = ROOT.RooGaussian("signal1", "signal1", jpsimass, mean, sigma1)
sig2 = ROOT.RooGaussian("signal2", "signal2", jpsimass, mean, sigma2)

a0 = ROOT.RooRealVar("a0", "a0", -0.3, -1, 1)
a1 = ROOT.RooRealVar("a1", "a1", 0.03, -1, +1)
bkg = ROOT.RooChebychev("bkg", "Background", jpsimass, [a0, a1]) ## to increase the degree, just increase the coefficients

bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.5, 0.4, 0.6)
barfrac = ROOT.RooRealVar("f_{barrel}", "fraction of higher res", 0.2, 0.1, 0.3)
model = ROOT.RooAddPdf("model", "g+cheb", [bkg, sig1, sig2], [bkgfrac, barfrac])
#model = ROOT.RooAddPdf("model", "g+g", [sig1, sig2], [barfrac])


# In[111]:


r_full = model.fitTo(mumuroohist, Save=True)


# In[112]:


c0 = ROOT.TCanvas("canvas0", "canvas0", 1200, 600)

mumuroodata.plotOn(jpsiframe)
model.plotOn(jpsiframe) # By default only fitted range is shown
model.paramOn(jpsiframe, ROOT.RooFit.Parameters([mean, barfrac, sigma1, sigma2]), ROOT.RooFit.Layout(0.6, 0.9, 0.9))

jpsiframe.Draw()
c0.Draw()


# In[113]:


print(mumuroodata.numEntries()*(1 - bkgfrac.getVal()))


# In[114]:


c0.Print("psi2s_fit_2022F.png")


# In[115]:


dfcands = pd.DataFrame(mumudata2.Filter("J_mass > 0.95 & J_mass < 1.1").AsNumpy(columns=["J_mass"]))


# In[116]:


len(dfcands)


# In[117]:


jpsimass = ROOT.RooRealVar("J_mass", "m(#mu#mu) [GeV]", 0.95, 1.1)
mumuroodata = ROOT.RooDataSet.from_numpy({"J_mass": dfcands["J_mass"]}, [jpsimass])
mumuroohist = mumuroodata.binnedClone()


# In[118]:


jpsiframe = jpsimass.frame(Title="Dimuon Candidate Mass")


# In[119]:


mean = ROOT.RooRealVar("#mu", "mean of gaussians", 1.01, 0.95, 1.05) #mean value, min value, max value
sigma1 = ROOT.RooRealVar("#sigma_{barrel}", "width1", 0.01, 0.005, 0.05) # most probably this are not that good!
sigma2 = ROOT.RooRealVar("#sigma_{endcap}", "width2", 0.03, 0.005, 0.05)

sig1 = ROOT.RooGaussian("signal1", "signal1", jpsimass, mean, sigma1)
sig2 = ROOT.RooGaussian("signal2", "signal2", jpsimass, mean, sigma2)

a0 = ROOT.RooRealVar("a0", "a0", -0.02, -1, 1)
a1 = ROOT.RooRealVar("a1", "a1", -0.1, -1, +1)
bkg = ROOT.RooChebychev("bkg", "Background", jpsimass, [a0, a1]) ## to increase the degree, just increase the coefficients

bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.7, 0.6, 0.8)
barfrac = ROOT.RooRealVar("f_{barrel}", "fraction of higher res", 0.2, 0.1, 0.3)
model = ROOT.RooAddPdf("model", "g+cheb", [bkg, sig1, sig2], [bkgfrac, barfrac])
#model = ROOT.RooAddPdf("model", "g+g", [sig1, sig2], [barfrac])


# In[120]:


r_full = model.fitTo(mumuroohist, Save=True)


# In[121]:


c0 = ROOT.TCanvas("canvas0", "canvas0", 1200, 600)

mumuroodata.plotOn(jpsiframe)
model.plotOn(jpsiframe) # By default only fitted range is shown
model.paramOn(jpsiframe, ROOT.RooFit.Parameters([mean, sigma1, sigma2, barfrac, bkgfrac]), ROOT.RooFit.Layout(0.6, 0.9, 0.9))

jpsiframe.Draw()
c0.Draw()


# In[122]:


c0.Print("phi_fit_2022F.png")


# In[300]:


dfcands = pd.DataFrame(mumudata2.Filter("J_mass > 9.0 & J_mass < 11. & J_pt > 10. & abs(J_eta) < 2.0").AsNumpy(columns=["J_mass"]))


# In[301]:


len(dfcands)


# In[302]:


upsmass = ROOT.RooRealVar("J_mass", "m(#mu#mu) [GeV]", 9.0, 11.)
mumuroodata = ROOT.RooDataSet.from_numpy({"J_mass": dfcands["J_mass"]}, [upsmass])
mumuroohist = mumuroodata.binnedClone()


# In[303]:


upsframe = upsmass.frame(Title="Dimuon Candidate Mass")


# In[304]:


meanu1 = ROOT.RooRealVar("#mu_{Y(1S)}", "mean of gaussians", 9.447, 9.44, 9.46) #mean value, min value, max value
meanu2 = ROOT.RooRealVar("#mu_{Y(2S)}", "mean of gaussians", 10.015, 9.9, 10.1) #mean value, min value, max value
meanu3 = ROOT.RooRealVar("#mu_{Y(3S)}", "mean of gaussians", 10.338, 10.2, 10.4) #mean value, min value, max value

meanu1.setConstant()
meanu2.setConstant()
meanu3.setConstant()

sigmau1 = ROOT.RooRealVar("#sigma_{Y(1S)}", "width1", 0.08, 0.05, 0.1) # most probably this are not that good!
sigmau2 = ROOT.RooRealVar("#sigma_{Y(2S)}", "width2", 0.08, 0.05, 0.1) # most probably this are not that good!
sigmau3 = ROOT.RooRealVar("#sigma_{Y(3S)}", "width3", 0.08, 0.05, 0.1) # most probably this are not that good!

#sigmau1.setConstant()
#sigmau2.setConstant()
#sigmau3.setConstant()

ncb = ROOT.RooRealVar("nu1", "n cb y1s", 1.5, 0., 3.0)
alphacb = ROOT.RooRealVar("alphau1", "alpha cb y1s", 5.0, -10., 100.0)

alphacb.setConstant()

#sigu1 = ROOT.RooGaussian("signal1", "signal1", upsmass, meanu1, sigmau1)
sigu1 = ROOT.RooCrystalBall("signal1", "signal1", upsmass, meanu1, sigmau1, ncb, alphacb)
sigu2 = ROOT.RooGaussian("signal2", "signal2", upsmass, meanu2, sigmau2)
sigu3 = ROOT.RooGaussian("signal3", "signal3", upsmass, meanu3, sigmau3)

f0 = ROOT.RooRealVar("f0", "f0", -0.3, -1, +0)
f1 = ROOT.RooRealVar("f1", "f1", -0.1, -1., +0)
bkgu = ROOT.RooChebychev("bkgu", "Background", upsmass, [f0, f1]) ## to increase the degree, just increase the coefficients

bkgfracu = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.38, 0.2, 0.6)
sigfrac1 = ROOT.RooRealVar("sigfrac1", "fraction of Y(1S)", 0.40, 0.2, 0.5)
sigfrac2 = ROOT.RooRealVar("sigfrac2", "fraction of Y(2S)", 0.14, 0.05, 0.2)

#sigfrac1.setConstant()
#sigfrac2.setConstant()

modelu = ROOT.RooAddPdf("modelu", "gs+cheb", [bkgu, sigu1, sigu2, sigu3], [bkgfracu, sigfrac1, sigfrac2])


# In[305]:


upsfit = modelu.fitTo(mumuroohist, Save=True)


# In[306]:


c2 = ROOT.TCanvas("canvas2", "canvas2", 1200, 600)

mumuroodata.plotOn(upsframe)
modelu.plotOn(upsframe) # By default only fitted range is shown
modelu.paramOn(upsframe, ROOT.RooFit.Parameters([meanu1, meanu2, meanu3, sigmau1, sigmau2, sigmau3, bkgfracu]), ROOT.RooFit.Layout(0.65, 0.9, 0.9))

upsframe.Draw()
c2.Draw()


# In[307]:


len(dfcands)*(1-bkgfracu.getVal())


# In[195]:


c2.Print("ups_fit_2022F.png")


# In[158]:


pthist0 = data.Histo1D(("J_pt", "dimuon pt", 250, 0, 30.), "J_pt")
pthist1 = mumudata.Histo1D(("J_pt", "dimuon pt", 250, 0, 30.), "J_pt")
pthist2 = mumudata.Filter("J_mass > 2.9 & J_mass < 3.3 & J_Prob > 0.1").Histo1D(("J_pt", "dimuon pt", 250, 0, 30.), "J_pt")


# In[114]:


etahist0 = data.Histo1D(("J_eta", "dimuon eta", 250, -3, 3.), "J_eta")
etahist1 = mumudata.Histo1D(("J_eta", "dimuon eta", 250, -3, 3.), "J_eta")
etahist2 = mumudata.Filter("J_mass > 2.9 & J_mass < 3.3 & J_Prob > 0.1").Histo1D(("J_eta", "dimuon pt", 250, 0, 30.), "J_eta")


# In[115]:


c0 = ROOT.TCanvas()
etahist0.Draw()
etahist1.Draw("same")
etahist2.Draw("same")
c0.Draw()


# In[56]:


jpsicut = data.Filter("J_mass > 2.8")\
                .Filter("J_mass < 4.0")\

jpsihist = jpsicut.Histo1D(("J_mass", "Jpsi mass", 250, 2.8, 4.0), "J_mass")


# In[57]:


canvas2 = ROOT.TCanvas()
canvas2.SetCanvasSize(1500, 500)
canvas2.SetLogy()
jpsihist.Draw()
canvas2.Draw()


# In[58]:


jpsicut = data.Filter("J_mass > 2.8")\
                .Filter("J_mass < 4.0")\

jpsihist = jpsicut.Histo1D(("J_mass", "Jpsi mass", 250, 2.8, 4.0), "J_mass")


# In[59]:


from ROOT import TF1

gaus1 = ROOT.TF1("gaus1", "gaus", 2.8, 4.0)
gaus2 = ROOT.TF1("gaus2", "gaus", 2.8, 4.0)
exp = ROOT.TF1("bkg", "[0]*exp(x/[1])", 2.8, 4.0)

gaus1.SetParameters(5e5, 3.1, 0.05)
gaus2.SetParameters(5e4, 3.7, 0.05)
exp.SetParameters(1e4, 1)

total = ROOT.TF1("total", "gaus1+gaus2+expo")
jpsihist.Fit("total")


# In[60]:


canvas2 = ROOT.TCanvas()
canvas2.SetCanvasSize(1500, 500)
canvas2.SetLogy()
jpsihist.Draw()
canvas2.Draw()


# In[61]:


upscut = data.Filter("J_mass > 8.5")\
                .Filter("J_mass < 11.5")
#upsDimacut = upscut.Filter("tri_LowMassInclusive == 1")

upshist = upscut.Histo1D(("J_mass", "Ups mass", 200, 8.5, 11.5), "J_mass")
#upsDimahist = upsDimacut.Histo1D(("J_mass", "Jpsi mass", 50, 8.5, 11.5), "J_mass")


# In[62]:


canvas3 = ROOT.TCanvas()
canvas3.SetCanvasSize(1500, 500)
upshist.Draw()
#upsDimahist.Draw("same")
canvas3.Draw()


# In[63]:


lowmasscut = data.Filter("J_mass > 0.15").Filter("J_mass < 1.15")
lowmassdispcut = lowmasscut.Filter("tri_LowMassDisplaced == 1")

lowmasshist = lowmasscut.Histo1D(("J_mass", "mumu mass", 100, 0.5, 0.6), "J_mass")
lowmassdisphist = lowmassdispcut.Histo1D(("J_mass", "mumu_mass", 100, 0.5, 0.6), "J_mass")


# In[64]:


lowmasshist.GetYaxis().SetRangeUser(0, 2000)
canvas4 = ROOT.TCanvas()
canvas4.SetCanvasSize(1500, 500)
lowmasshist.Draw()
#lowmassdisphist.Draw("same")
canvas4.Draw()


# In[59]:


mumudataC = dataC.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")\
.Filter("J_pt1 > 4. & J_pt2 > 4. & J_pt > 8 & J_Prob > 0.1 & abs(J_eta) < 2.4")\
.Filter("J_mass > 2.9 & J_mass < 3.3")

mumudataD = dataD.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")\
.Filter("J_pt1 > 4. & J_pt2 > 4. & J_pt > 8 & J_Prob > 0.1 & abs(J_eta) < 2.4")\
.Filter("J_mass > 2.9 & J_mass < 3.3")

mumudataE = dataE.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")\
.Filter("J_pt1 > 4. & J_pt2 > 4. & J_pt > 8 & J_Prob > 0.1 & abs(J_eta) < 2.4")\
.Filter("J_mass > 2.9 & J_mass < 3.3")

mumudataF = data.Filter("mu1HLTmatched == 1 & mu2HLTmatched == 1 & J_mass < 11.5")\
.Filter("J_pt1 > 4. & J_pt2 > 4. & J_pt > 8 & J_Prob > 0.1 & abs(J_eta) < 2.4")\
.Filter("J_mass > 2.9 & J_mass < 3.3")


# In[60]:


jpsihistC = mumudataC.Histo1D(("J_mass", "J/#psi mass distribution (norm.)", 100, 2.9, 3.3), "J_mass")
jpsihistD = mumudataD.Histo1D(("J_mass", "Jpsi mass", 100, 2.9, 3.3), "J_mass")
jpsihistE = mumudataE.Histo1D(("J_mass", "Jpsi mass", 100, 2.9, 3.3), "J_mass")
jpsihistF = mumudataF.Histo1D(("J_mass", "Jpsi mass", 100, 2.9, 3.3), "J_mass")


# In[61]:


c1 = ROOT.TCanvas()

legend = ROOT.TLegend(0.78, 0.6, 0.89, 0.75)

jpsihistC.Scale(1./jpsihistC.Integral())
jpsihistD.Scale(1./jpsihistD.Integral())
jpsihistE.Scale(1./jpsihistE.Integral())
jpsihistF.Scale(1./jpsihistF.Integral())

jpsihistC.SetLineWidth(2)
jpsihistD.SetLineWidth(2)
jpsihistE.SetLineWidth(2)
jpsihistF.SetLineWidth(2)

jpsihistC.SetLineColor(2)
jpsihistD.SetLineColor(3)
jpsihistE.SetLineColor(4)
jpsihistF.SetLineColor(5)

legend.AddEntry(jpsihistC.GetPtr(), "2022C", "l")
legend.AddEntry(jpsihistD.GetPtr(), "2022D", "l")
legend.AddEntry(jpsihistE.GetPtr(), "2022E", "l")
legend.AddEntry(jpsihistF.GetPtr(), "2022F", "l")

jpsihistC.Draw("")
jpsihistD.Draw("same")
jpsihistE.Draw("same")
jpsihistF.Draw("same")

legend.Draw("")

c1.Draw("")


# In[62]:


jpsihistC = mumudataC.Histo1D(("J_mass", "J/#psi p_{T} distribution (norm.)", 100, 8., 50.), "J_pt")
jpsihistD = mumudataD.Histo1D(("J_mass", "Jpsi mass", 100, 8., 50.), "J_pt")
jpsihistE = mumudataE.Histo1D(("J_mass", "Jpsi mass", 100, 8., 50.), "J_pt")
jpsihistF = mumudataF.Histo1D(("J_mass", "Jpsi mass", 100, 8., 50.), "J_pt")


# In[63]:


c1 = ROOT.TCanvas()

legend = ROOT.TLegend(0.78, 0.6, 0.89, 0.75)

jpsihistC.Scale(1./jpsihistC.Integral())
jpsihistD.Scale(1./jpsihistD.Integral())
jpsihistE.Scale(1./jpsihistE.Integral())
jpsihistF.Scale(1./jpsihistF.Integral())

jpsihistC.SetLineWidth(2)
jpsihistD.SetLineWidth(2)
jpsihistE.SetLineWidth(2)
jpsihistF.SetLineWidth(2)

jpsihistC.SetLineColor(2)
jpsihistD.SetLineColor(3)
jpsihistE.SetLineColor(4)
jpsihistF.SetLineColor(5)

legend.AddEntry(jpsihistC.GetPtr(), "2022C", "l")
legend.AddEntry(jpsihistD.GetPtr(), "2022D", "l")
legend.AddEntry(jpsihistE.GetPtr(), "2022E", "l")
legend.AddEntry(jpsihistF.GetPtr(), "2022F", "l")

jpsihistC.Draw("")
jpsihistD.Draw("same")
jpsihistE.Draw("same")
jpsihistF.Draw("same")

legend.Draw("")

c1.Draw("")


# In[64]:


jpsihistC = mumudataC.Histo1D(("J_mass", "leading #mu p_{T} distribution (norm.)", 100, 3., 30.), "J_pt1")
jpsihistD = mumudataD.Histo1D(("J_mass", "Jpsi mass", 100, 3., 30.), "J_pt1")
jpsihistE = mumudataE.Histo1D(("J_mass", "Jpsi mass", 100, 3., 30.), "J_pt1")
jpsihistF = mumudataF.Histo1D(("J_mass", "Jpsi mass", 100, 3., 30.), "J_pt1")


# In[65]:


c1 = ROOT.TCanvas()

legend = ROOT.TLegend(0.78, 0.6, 0.89, 0.75)

jpsihistC.Scale(1./jpsihistC.Integral())
jpsihistD.Scale(1./jpsihistD.Integral())
jpsihistE.Scale(1./jpsihistE.Integral())
jpsihistF.Scale(1./jpsihistF.Integral())

jpsihistC.SetLineWidth(2)
jpsihistD.SetLineWidth(2)
jpsihistE.SetLineWidth(2)
jpsihistF.SetLineWidth(2)

jpsihistC.SetLineColor(2)
jpsihistD.SetLineColor(3)
jpsihistE.SetLineColor(4)
jpsihistF.SetLineColor(5)

legend.AddEntry(jpsihistC.GetPtr(), "2022C", "l")
legend.AddEntry(jpsihistD.GetPtr(), "2022D", "l")
legend.AddEntry(jpsihistE.GetPtr(), "2022E", "l")
legend.AddEntry(jpsihistF.GetPtr(), "2022F", "l")

jpsihistC.Draw("")
jpsihistD.Draw("same")
jpsihistE.Draw("same")
jpsihistF.Draw("same")

legend.Draw("")

c1.Draw("")


# In[66]:


jpsihistC = mumudataC.Histo1D(("J_mass", "J/#psi #eta distribution (norm.)", 100, -2.5, 2.5), "J_eta")
jpsihistD = mumudataD.Histo1D(("J_mass", "J/#psi #eta distribution (norm.)", 100, -2.5, 2.5), "J_eta")
jpsihistE = mumudataE.Histo1D(("J_mass", "Jpsi mass", 100, -2.5, 2.5), "J_eta")
jpsihistF = mumudataF.Histo1D(("J_mass", "Jpsi mass", 100, -2.5, 2.5), "J_eta")


# In[70]:


c1 = ROOT.TCanvas()

legend = ROOT.TLegend(0.78, 0.6, 0.89, 0.75)

jpsihistC.Scale(1./jpsihistC.Integral())
jpsihistD.Scale(1./jpsihistD.Integral())
jpsihistE.Scale(1./jpsihistE.Integral())
jpsihistF.Scale(1./jpsihistF.Integral())

jpsihistC.SetLineWidth(2)
jpsihistD.SetLineWidth(2)
jpsihistE.SetLineWidth(2)
jpsihistF.SetLineWidth(2)

jpsihistC.SetLineColor(2)
jpsihistD.SetLineColor(3)
jpsihistE.SetLineColor(4)
jpsihistF.SetLineColor(5)

legend.AddEntry(jpsihistC.GetPtr(), "2022C", "l")
legend.AddEntry(jpsihistD.GetPtr(), "2022D", "l")
legend.AddEntry(jpsihistE.GetPtr(), "2022E", "l")
legend.AddEntry(jpsihistF.GetPtr(), "2022F", "l")

jpsihistD.SetTitle("J/#psi #eta distribution (norm.)")

jpsihistD.Draw("")
jpsihistC.Draw("same")
jpsihistE.Draw("same")
jpsihistF.Draw("same")

legend.Draw("")

c1.Draw("")


# In[ ]:




