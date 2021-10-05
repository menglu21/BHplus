import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

class BHProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("lhe_nlepton", "I")
    self.out.branch("n_tight_muon", "I")
    self.out.branch("n_loose_muon", "I")
    self.out.branch("n_tight_ele", "I")
    self.out.branch("n_loose_ele", "I")
    self.out.branch("n_tight_jet_in24", "I")
    self.out.branch("n_tight_jet_in47", "I")
    self.out.branch("n_bjet_DeepB", "I")
    self.out.branch("btag_SFall", "F")
    self.out.branch("HT_24", "F")
    self.out.branch("HT_47", "F")
    self.out.branch("j1_pt", "F")
    self.out.branch("j1_eta", "F")
    self.out.branch("j1_phi", "F")
    self.out.branch("j1_mass", "F")
    self.out.branch("j2_pt", "F")
    self.out.branch("j2_eta", "F")
    self.out.branch("j2_phi", "F")
    self.out.branch("j2_mass", "F")
    self.out.branch("j3_pt", "F")
    self.out.branch("j3_eta", "F")
    self.out.branch("j3_phi", "F")
    self.out.branch("j3_mass", "F")
    self.out.branch("j4_pt", "F")
    self.out.branch("j4_eta", "F")
    self.out.branch("j4_phi", "F")
    self.out.branch("j4_mass", "F")
    self.out.branch("DeepB_j1_pt", "F")
    self.out.branch("DeepB_j1_eta", "F")
    self.out.branch("DeepB_j1_phi", "F")
    self.out.branch("DeepB_j1_mass", "F")
    self.out.branch("DeepB_j2_pt", "F")
    self.out.branch("DeepB_j2_eta", "F")
    self.out.branch("DeepB_j2_phi", "F")
    self.out.branch("DeepB_j2_mass", "F")
    self.out.branch("DeepB_j3_pt", "F")
    self.out.branch("DeepB_j3_eta", "F")
    self.out.branch("DeepB_j3_phi", "F")
    self.out.branch("DeepB_j3_mass", "F")
    self.out.branch("bh_nl", "B")
    self.out.branch("bh_jets", "B")
    self.out.branch("bh_region", "I")
    self.out.branch("bh_l1_id", "I")
    self.out.branch("bh_l1_pdgid", "I")
    self.out.branch("bh_l1_pt", "F")
    self.out.branch("bh_l1_eta", "F")
    self.out.branch("bh_l1_phi", "F")
    self.out.branch("bh_l1_mass", "F")
    self.out.branch("bh_met", "F")
    self.out.branch("bh_met_phi", "F")
    self.out.branch("bh_dr_l1j1", "F")
    self.out.branch("bh_dr_l1j2", "F")
    self.out.branch("bh_dr_l1j3", "F")
    self.out.branch("bh_dr_l1j4", "F")
    self.out.branch("bh_mlj1", "F")
    self.out.branch("bh_mlj2", "F")
    self.out.branch("bh_mlj3", "F")
    self.out.branch("bh_mlj4", "F")
    self.out.branch("bh_mlj1j2", "F")
    self.out.branch("bh_mlj1j3", "F")
    self.out.branch("bh_mlj1j4", "F")
    self.out.branch("bh_mlj2j3", "F")
    self.out.branch("bh_mlj2j4", "F")
    self.out.branch("bh_mlj3j4", "F")
    self.out.branch("bh_mlb1", "F")
    self.out.branch("bh_mlb2", "F")
    self.out.branch("bh_mlb3", "F")
    self.out.branch("bh_mlb1b2", "F")
    self.out.branch("bh_mlb1b3", "F")
    self.out.branch("bh_mlb2b3", "F")
    self.out.branch("WZ_region", "I")
    self.out.branch("WZ_zl1_id", "I")
    self.out.branch("WZ_zl2_id", "I")
    self.out.branch("WZ_wl_id", "I")
    self.out.branch("WZ_zl1_pdgid", "I")
    self.out.branch("WZ_zl2_pdgid", "I")
    self.out.branch("WZ_wl_pdgid", "I")
    self.out.branch("WZ_zl1_pt", "F")
    self.out.branch("WZ_zl1_eta", "F")
    self.out.branch("WZ_zl1_phi", "F")
    self.out.branch("WZ_zl1_mass", "F")
    self.out.branch("WZ_zl2_pt", "F")
    self.out.branch("WZ_zl2_eta", "F")
    self.out.branch("WZ_zl2_phi", "F")
    self.out.branch("WZ_zl2_mass", "F")
    self.out.branch("WZ_l3_pt", "F")
    self.out.branch("WZ_l3_eta", "F")
    self.out.branch("WZ_l3_phi", "F")
    self.out.branch("WZ_l3_mass", "F")
    self.out.branch("WZ_Z_mass", "F")
    self.out.branch("WZ_Z_pt", "F")
    self.out.branch("WZ_Z_eta", "F")
    self.out.branch("WZ_Z_phi", "F")
    self.out.branch("WZ_met", "F")
    self.out.branch("DY_region", "I")
    self.out.branch("DY_l1_id", "I")
    self.out.branch("DY_l2_id", "I")
    self.out.branch("DY_l1_pdgid", "I")
    self.out.branch("DY_l2_pdgid", "I")
    self.out.branch("DY_l1_pt", "F")
    self.out.branch("DY_l1_eta", "F")
    self.out.branch("DY_l1_phi", "F")
    self.out.branch("DY_l1_mass", "F")
    self.out.branch("DY_l2_pt", "F")
    self.out.branch("DY_l2_eta", "F")
    self.out.branch("DY_l2_phi", "F")
    self.out.branch("DY_l2_mass", "F")
    self.out.branch("DY_z_mass", "F")
    self.out.branch("DY_z_pt", "F")
    self.out.branch("DY_z_eta", "F")
    self.out.branch("DY_z_phi", "F")
    self.out.branch("DY_drll", "F")
    self.out.branch("tightJets_id_in24","I",lenVar="nJet")
    self.out.branch("tightJets_id_in47","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_c_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightElectrons_id","I",lenVar="nElectron")
    self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    self.out.branch("tightMuons_id","I",lenVar="nMuon")
    self.out.branch("additional_looseMuons_id","I",lenVar="nMuon")
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
	for iobj in range(0,event.nTrigObj):
	  if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
	    HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    lhe_nlepton=0
    if self.is_lhe:
      lheparticle = Collection(event, 'LHEPart')
      for ilhe in range(0, event.nLHEPart):
	if lheparticle[ilhe].status==1 and (abs(lheparticle[ilhe].pdgId)==11 or abs(lheparticle[ilhe].pdgId)==13 or abs(lheparticle[ilhe].pdgId)==15):
	  lhe_nlepton=lhe_nlepton+1

    self.out.fillBranch("lhe_nlepton", lhe_nlepton)

    # total number of ele+muon, currently require at least 1 leptons
    if ((event.nMuon + event.nElectron) < 1): return False

    # lepton pt threshold according to the HLT
    if self.year=="2016":
      ele_pt=30
      muon_pt=26
    if self.year=="2017":
      ele_pt=35
      muon_pt=30
    if self.year=="2018":
      ele_pt=35
      muon_pt=26

    # Muon selection: tight cut-based ID + tight PF iso, or loose cut-based ID + loose PF iso, with pt > 20 GeV
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    tightMuons = []
    tightMuons_pdgid = []
    tightMuons_id = []
    additional_looseMuons = []
    additional_looseMuons_pdgid = []
    additional_looseMuons_id = []
    for imu in range(0, event.nMuon):
      if (muons[imu].tightId):
        if (muons[imu].pfRelIso04_all<0.15 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>15):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
      elif (muons[imu].looseId):
        if (muons[imu].pfRelIso04_all<0.25 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>15):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          additional_looseMuons.append(muon_v4_temp.Clone())
          additional_looseMuons_pdgid.append(muons[imu].pdgId)
          additional_looseMuons_id.append(imu)

    n_tight_muon = len(tightMuons)
    n_loose_muon = len(additional_looseMuons)

    self.out.fillBranch("n_tight_muon", n_tight_muon)
    self.out.fillBranch("n_loose_muon", n_loose_muon)
    tightMuons_id.extend(np.zeros(event.nMuon-len(tightMuons_id),int)-1)
    additional_looseMuons_id.extend(np.zeros(event.nMuon-len(additional_looseMuons_id),int)-1)
    self.out.fillBranch("tightMuons_id", tightMuons_id)
    self.out.fillBranch("additional_looseMuons_id", additional_looseMuons_id)

    # electron selection: tight (veto) cut-based ID + impact parameter cut, with pt > 15 GeV
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    tightElectrons = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []
    for iele in range(0, event.nElectron):
      if (electrons[iele].cutBased==4):
        if (((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz)<0.1) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>15):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          tightElectrons.append(electron_v4_temp.Clone())
          tightElectrons_pdgid.append(electrons[iele].pdgId)
          tightElectrons_id.append(iele)
      elif (electrons[iele].cutBased>0):
        if (((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz)<0.1) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>15):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          additional_vetoElectrons.append(electron_v4_temp.Clone())
          additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
          additional_vetoElectrons_id.append(iele)

    n_tight_ele = len(tightElectrons)
    n_loose_ele = len(additional_vetoElectrons)
    self.out.fillBranch("n_tight_ele", n_tight_ele)
    self.out.fillBranch("n_loose_ele", n_loose_ele)
    tightElectrons_id.extend(np.zeros(event.nElectron-len(tightElectrons_id),int)-1)
    additional_vetoElectrons_id.extend(np.zeros(event.nElectron-len(additional_vetoElectrons_id),int)-1)
    self.out.fillBranch("tightElectrons_id", tightElectrons_id)
    self.out.fillBranch("additional_vetoElectrons_id", additional_vetoElectrons_id)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets
    # tightLepVeto PF jets (ak4), UL 2016/2017/2018 (jetId 110=6), medium B-tag WP
    # UL17 DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # UL18 DeepCSV=(nanoaod btagDeepB) loose: 0.1208, medium: 0.4168, tight: 0.7665
    # UL17 DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476
    # UL18 DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0490, medium: 0.2783, tight: 0.7100

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # UL17 CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
    # UL18 CvsL=btagDeepFlavCvL: 0.099, CvsB=btagDeepFlavCvB: 0.325
    # c-tag not available in NANOAOD yet

    jets = Collection(event, 'Jet')

    j1_pt=-99
    j1_eta=-99
    j1_phi=-99
    j1_mass=-99
    j2_pt=-99
    j2_eta=-99
    j2_phi=-99
    j2_mass=-99
    j3_pt=-99
    j3_eta=-99
    j3_phi=-99
    j3_mass=-99
    j4_pt=-99
    j4_eta=-99
    j4_phi=-99
    j4_mass=-99
    DeepB_j1_pt=-99
    DeepB_j1_eta=-99
    DeepB_j1_phi=-99
    DeepB_j1_mass=-99
    DeepB_j2_pt=-99
    DeepB_j2_eta=-99
    DeepB_j2_phi=-99
    DeepB_j2_mass=-99
    DeepB_j3_pt=-99
    DeepB_j3_eta=-99
    DeepB_j3_phi=-99
    DeepB_j3_mass=-99

    tightJets_id_in24 = []
    tightJets_id_in47 = []

    tightJets_b_DeepCSVmedium_id = []

    # require DeltaR between Jets and tight leptons greater than 0.4
    jet_v4_temp=TLorentzVector()
    for ijet in range(0, event.nJet):
      pass_jet_lep_Dr=1
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass[ijet])
      for ilep in range(0,len(tightLeptons)):
	if jet_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_jet_lep_Dr=0

      if not (pass_jet_lep_Dr>0):continue

      if self.year=="2016":
        if jets[ijet].jetId==6 and event.Jet_pt[ijet]>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4: 
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
	    tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepB > 0.4941):
              tightJets_b_DeepCSVmedium_id.append(ijet)

      if (self.year=="2017"):
	if jets[ijet].jetId==6 and event.Jet_pt[ijet]>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4:
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
            tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepFlavB > 0.3040):
              tightJets_b_DeepCSVmedium_id.append(ijet)

      if (self.year=="2018"):
	if jets[ijet].jetId==6 and event.Jet_pt[ijet]>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4:
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
            tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepFlavB > 0.2783):
              tightJets_b_DeepCSVmedium_id.append(ijet)

    HT_24=0
    HT_47=0
    for ijet in range(0,len(tightJets_id_in24)):
      HT_24=HT_24+event.Jet_pt[tightJets_id_in24[ijet]]

    HT_47=HT_24
    for ijet in range(0,len(tightJets_id_in47)):
      HT_47=HT_47+event.Jet_pt[tightJets_id_in47[ijet]]

    self.out.fillBranch("HT_24",HT_24)
    self.out.fillBranch("HT_47",HT_47)

    n_tight_jet_in24 = len(tightJets_id_in24)
    n_tight_jet_in47 = len(tightJets_id_in24)+len(tightJets_id_in47)
    n_bjet_DeepB = len(tightJets_b_DeepCSVmedium_id)

    btag_SFall=1.
    if n_bjet_DeepB>0 and self.is_mc:
      for ib in range(0,n_bjet_DeepB):
        btag_SFall=btag_SFall*event.Jet_btagSF_deepjet_M[tightJets_b_DeepCSVmedium_id[ib]]

    self.out.fillBranch("n_tight_jet_in24",n_tight_jet_in24)
    self.out.fillBranch("n_tight_jet_in47",n_tight_jet_in47)
    self.out.fillBranch("n_bjet_DeepB",n_bjet_DeepB)
    self.out.fillBranch("btag_SFall",btag_SFall)

    if n_tight_jet_in24>3:
      j4_pt=event.Jet_pt[tightJets_id_in24[3]]
      j4_eta=event.Jet_eta[tightJets_id_in24[3]]
      j4_phi=event.Jet_phi[tightJets_id_in24[3]]
      j4_mass=event.Jet_mass[tightJets_id_in24[3]]
      j3_pt=event.Jet_pt[tightJets_id_in24[2]]
      j3_eta=event.Jet_eta[tightJets_id_in24[2]]
      j3_phi=event.Jet_phi[tightJets_id_in24[2]]
      j3_mass=event.Jet_mass[tightJets_id_in24[2]]
      j2_pt=event.Jet_pt[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass[tightJets_id_in24[0]]
    if n_tight_jet_in24==3:
      j3_pt=event.Jet_pt[tightJets_id_in24[2]]
      j3_eta=event.Jet_eta[tightJets_id_in24[2]]
      j3_phi=event.Jet_phi[tightJets_id_in24[2]]
      j3_mass=event.Jet_mass[tightJets_id_in24[2]]
      j2_pt=event.Jet_pt[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass[tightJets_id_in24[0]]
    if n_tight_jet_in24==2:
      j2_pt=event.Jet_pt[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass[tightJets_id_in24[0]]
    if n_tight_jet_in24==1:
      j1_pt=event.Jet_pt[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass[tightJets_id_in24[0]]

    if n_bjet_DeepB>2:
      DeepB_j1_pt=event.Jet_pt[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_mass=event.Jet_mass[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j2_pt=event.Jet_pt[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_mass=event.Jet_mass[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j3_pt=event.Jet_pt[tightJets_b_DeepCSVmedium_id[2]]
      DeepB_j3_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[2]]
      DeepB_j3_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[2]]
      DeepB_j3_mass=event.Jet_mass[tightJets_b_DeepCSVmedium_id[2]]
    if n_bjet_DeepB==2:
      DeepB_j1_pt=event.Jet_pt[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_mass=event.Jet_mass[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j2_pt=event.Jet_pt[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_mass=event.Jet_mass[tightJets_b_DeepCSVmedium_id[1]]
    if n_bjet_DeepB==1:
      DeepB_j1_pt=event.Jet_pt[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_mass=event.Jet_mass[tightJets_b_DeepCSVmedium_id[0]]

    self.out.fillBranch("j1_pt",j1_pt)
    self.out.fillBranch("j1_eta",j1_eta)
    self.out.fillBranch("j1_phi",j1_phi)
    self.out.fillBranch("j1_mass",j1_mass)
    self.out.fillBranch("j2_pt",j2_pt)
    self.out.fillBranch("j2_eta",j2_eta)
    self.out.fillBranch("j2_phi",j2_phi)
    self.out.fillBranch("j2_mass",j2_mass)
    self.out.fillBranch("j3_pt",j3_pt)
    self.out.fillBranch("j3_eta",j3_eta)
    self.out.fillBranch("j3_phi",j3_phi)
    self.out.fillBranch("j3_mass",j3_mass)
    self.out.fillBranch("j4_pt",j4_pt)
    self.out.fillBranch("j4_eta",j4_eta)
    self.out.fillBranch("j4_phi",j4_phi)
    self.out.fillBranch("j4_mass",j4_mass)
    self.out.fillBranch("DeepB_j1_pt",DeepB_j1_pt)
    self.out.fillBranch("DeepB_j1_eta",DeepB_j1_eta)
    self.out.fillBranch("DeepB_j1_phi",DeepB_j1_phi)
    self.out.fillBranch("DeepB_j1_mass",DeepB_j1_mass)
    self.out.fillBranch("DeepB_j2_pt",DeepB_j2_pt)
    self.out.fillBranch("DeepB_j2_eta",DeepB_j2_eta)
    self.out.fillBranch("DeepB_j2_phi",DeepB_j2_phi)
    self.out.fillBranch("DeepB_j2_mass",DeepB_j2_mass)
    self.out.fillBranch("DeepB_j3_pt",DeepB_j3_pt)
    self.out.fillBranch("DeepB_j3_eta",DeepB_j3_eta)
    self.out.fillBranch("DeepB_j3_phi",DeepB_j3_phi)
    self.out.fillBranch("DeepB_j3_mass",DeepB_j3_mass)

    tightJets_id_in24.extend(np.zeros(event.nJet-len(tightJets_id_in24),int)-1)
    tightJets_id_in47.extend(np.zeros(event.nJet-len(tightJets_id_in47),int)-1)
    tightJets_b_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVmedium_id),int)-1)

    self.out.fillBranch("tightJets_id_in24",tightJets_id_in24)
    self.out.fillBranch("tightJets_id_in47",tightJets_id_in47)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)

    if len(tightLeptons)<1:return False

    #  BH region (the requirement on bjet will separate this region to signal region and w-region)
    #region: only 1 tight leptons, at least three jets
    #bh region lepton number selections
    bh_nl=False
    #bh region jet and bjet selection
    bh_jets=False
    #bh region tag, 1: muon, 2:ele
    bh_region=0
    bh_l1_id=-1
    bh_l1_pdgid=-99
    bh_l1_pt=-99
    bh_l1_eta=-99
    bh_l1_phi=-99
    bh_l1_mass=-99
    bh_met=-99
    bh_met_phi=-99
    bh_dr_l1j1=-99
    bh_dr_l1j2=-99
    bh_dr_l1j3=-99
    bh_dr_l1j4=-99
    bh_mlj1=-99
    bh_mlj2=-99
    bh_mlj3=-99
    bh_mlj4=-99
    bh_mlj1j2=-99
    bh_mlj1j3=-99
    bh_mlj1j4=-99
    bh_mlj2j3=-99
    bh_mlj2j4=-99
    bh_mlj3j4=-99
    bh_mlb1=-99
    bh_mlb2=-99
    bh_mlb3=-99
    bh_mlb1b2=-99
    bh_mlb1b3=-99
    bh_mlb2b3=-99
    
    # 2th lepton veto
    if len(tightLeptons)==1 and len(looseLeptons)==0:
      if (n_tight_muon==1 and tightMuons[0]>muon_pt) or (n_tight_ele==1 and tightElectrons[0]>ele_pt):
        bh_nl=True
    # at least three jets (bjet requirement will applied at plot level)
    if bh_nl and n_tight_jet_in24>2:
      bh_jets=True

    if bh_nl:
      if self.is_mc:
        bh_met=event.MET_T1Smear_pt
	bh_met_phi=event.MET_T1Smear_phi
      else:
        bh_met=event.MET_T1_pt
        bh_met_phi=event.MET_T1_phi
      if len(tightElectrons)==0:
	bh_region=1
	bh_l1_id=tightMuons_id[0]
	bh_l1_pdgid=tightMuons_pdgid[0]
	bh_l1_pt=tightMuons[0].Pt()
	bh_l1_eta=tightMuons[0].Eta()
	bh_l1_phi=tightMuons[0].Phi()
	bh_l1_mass=tightMuons[0].M()
      if len(tightElectrons)==1:
	bh_region=2
	bh_l1_id=tightElectrons_id[0]
	bh_l1_pdgid=tightElectrons_pdgid[0]
        bh_l1_pt=tightElectrons[0].Pt()
        bh_l1_eta=tightElectrons[0].Eta()
        bh_l1_phi=tightElectrons[0].Phi()
        bh_l1_mass=tightElectrons[0].M()

    if bh_region>0:
      l1_v4_temp=TLorentzVector()
      l1_v4_temp.SetPtEtaPhiM(bh_l1_pt,bh_l1_eta,bh_l1_phi,bh_l1_mass)
      j1_v4_temp=TLorentzVector()
      j2_v4_temp=TLorentzVector()
      j3_v4_temp=TLorentzVector()
      j4_v4_temp=TLorentzVector()
      b1_v4_temp=TLorentzVector()
      b2_v4_temp=TLorentzVector()
      b3_v4_temp=TLorentzVector()
      if n_tight_jet_in24>3:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
	j2_v4_temp.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
	j3_v4_temp.SetPtEtaPhiM(j3_pt,j3_eta,j3_phi,j3_mass)
	j4_v4_temp.SetPtEtaPhiM(j4_pt,j4_eta,j4_phi,j4_mass)
	bh_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
	bh_dr_l1j2=l1_v4_temp.DeltaR(j2_v4_temp)
	bh_dr_l1j3=l1_v4_temp.DeltaR(j3_v4_temp)
	bh_dr_l1j4=l1_v4_temp.DeltaR(j4_v4_temp)
	bh_mlj1=(l1_v4_temp+j1_v4_temp).M()
	bh_mlj2=(l1_v4_temp+j2_v4_temp).M()
	bh_mlj3=(l1_v4_temp+j3_v4_temp).M()
	bh_mlj4=(l1_v4_temp+j4_v4_temp).M()
	bh_mlj1j2=(l1_v4_temp+j1_v4_temp+j2_v4_temp).M()
	bh_mlj1j3=(l1_v4_temp+j1_v4_temp+j3_v4_temp).M()
	bh_mlj1j4=(l1_v4_temp+j1_v4_temp+j4_v4_temp).M()
	bh_mlj2j3=(l1_v4_temp+j2_v4_temp+j3_v4_temp).M()
	bh_mlj2j4=(l1_v4_temp+j2_v4_temp+j4_v4_temp).M()
	bh_mlj3j4=(l1_v4_temp+j3_v4_temp+j4_v4_temp).M()
      if n_tight_jet_in24==3:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
        j2_v4_temp.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
        j3_v4_temp.SetPtEtaPhiM(j3_pt,j3_eta,j3_phi,j3_mass)
        bh_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
        bh_dr_l1j2=l1_v4_temp.DeltaR(j2_v4_temp)
        bh_dr_l1j3=l1_v4_temp.DeltaR(j3_v4_temp)
        bh_mlj1=(l1_v4_temp+j1_v4_temp).M()
        bh_mlj2=(l1_v4_temp+j2_v4_temp).M()
        bh_mlj3=(l1_v4_temp+j3_v4_temp).M()
	bh_mlj1j2=(l1_v4_temp+j1_v4_temp+j2_v4_temp).M()
	bh_mlj1j3=(l1_v4_temp+j1_v4_temp+j3_v4_temp).M()
	bh_mlj2j3=(l1_v4_temp+j2_v4_temp+j3_v4_temp).M()
      if n_tight_jet_in24==2:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
        j2_v4_temp.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
        bh_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
        bh_dr_l1j2=l1_v4_temp.DeltaR(j2_v4_temp)
        bh_mlj1=(l1_v4_temp+j1_v4_temp).M()
        bh_mlj2=(l1_v4_temp+j2_v4_temp).M()
	bh_mlj1j2=(l1_v4_temp+j1_v4_temp+j2_v4_temp).M()
      if n_tight_jet_in24==1:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
        bh_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
        bh_mlj1=(l1_v4_temp+j1_v4_temp).M()

      if n_bjet_DeepB>2:
	b1_v4_temp.SetPtEtaPhiM(DeepB_j1_pt, DeepB_j1_eta,DeepB_j1_phi,DeepB_j1_mass)
	b2_v4_temp.SetPtEtaPhiM(DeepB_j2_pt, DeepB_j2_eta,DeepB_j2_phi,DeepB_j2_mass)
	b3_v4_temp.SetPtEtaPhiM(DeepB_j3_pt, DeepB_j3_eta,DeepB_j3_phi,DeepB_j3_mass)
	bh_mlb1=(l1_v4_temp+b1_v4_temp).M()
	bh_mlb2=(l1_v4_temp+b2_v4_temp).M()
	bh_mlb3=(l1_v4_temp+b3_v4_temp).M()
	bh_mlb1b2=(l1_v4_temp+b1_v4_temp+b2_v4_temp).M()
	bh_mlb1b3=(l1_v4_temp+b1_v4_temp+b3_v4_temp).M()
	bh_mlb2b3=(l1_v4_temp+b2_v4_temp+b3_v4_temp).M()
      if n_bjet_DeepB==2:
	b1_v4_temp.SetPtEtaPhiM(DeepB_j1_pt, DeepB_j1_eta,DeepB_j1_phi,DeepB_j1_mass)
	b2_v4_temp.SetPtEtaPhiM(DeepB_j2_pt, DeepB_j2_eta,DeepB_j2_phi,DeepB_j2_mass)
	bh_mlb1=(l1_v4_temp+b1_v4_temp).M()
	bh_mlb2=(l1_v4_temp+b2_v4_temp).M()
	bh_mlb1b2=(l1_v4_temp+b1_v4_temp+b2_v4_temp).M()
      if n_bjet_DeepB==1:
	b1_v4_temp.SetPtEtaPhiM(DeepB_j1_pt, DeepB_j1_eta,DeepB_j1_phi,DeepB_j1_mass)
	bh_mlb1=(l1_v4_temp+b1_v4_temp).M()

    self.out.fillBranch("bh_nl", bh_nl)
    self.out.fillBranch("bh_jets", bh_jets)
    self.out.fillBranch("bh_region", bh_region)
    self.out.fillBranch("bh_l1_id", bh_l1_id)
    self.out.fillBranch("bh_l1_pdgid", bh_l1_pdgid)
    self.out.fillBranch("bh_l1_pt", bh_l1_pt)
    self.out.fillBranch("bh_l1_eta", bh_l1_eta)
    self.out.fillBranch("bh_l1_phi", bh_l1_phi)
    self.out.fillBranch("bh_l1_mass", bh_l1_mass)
    self.out.fillBranch("bh_met", bh_met)
    self.out.fillBranch("bh_met_phi", bh_met_phi)
    self.out.fillBranch("bh_dr_l1j1", bh_dr_l1j1)
    self.out.fillBranch("bh_dr_l1j2", bh_dr_l1j2)
    self.out.fillBranch("bh_dr_l1j3", bh_dr_l1j3)
    self.out.fillBranch("bh_dr_l1j4", bh_dr_l1j4)
    self.out.fillBranch("bh_mlj1", bh_mlj1)
    self.out.fillBranch("bh_mlj2", bh_mlj2)
    self.out.fillBranch("bh_mlj3", bh_mlj3)
    self.out.fillBranch("bh_mlj4", bh_mlj4)
    self.out.fillBranch("bh_mlj1j2", bh_mlj1j2)
    self.out.fillBranch("bh_mlj1j3", bh_mlj1j3)
    self.out.fillBranch("bh_mlj1j4", bh_mlj1j4)
    self.out.fillBranch("bh_mlj2j3", bh_mlj2j3)
    self.out.fillBranch("bh_mlj2j4", bh_mlj2j4)
    self.out.fillBranch("bh_mlj3j4", bh_mlj3j4)
    self.out.fillBranch("bh_mlb1", bh_mlb1)
    self.out.fillBranch("bh_mlb2", bh_mlb2)
    self.out.fillBranch("bh_mlb3", bh_mlb3)
    self.out.fillBranch("bh_mlb1b2", bh_mlb1b2)
    self.out.fillBranch("bh_mlb1b3", bh_mlb1b3)
    self.out.fillBranch("bh_mlb2b3", bh_mlb2b3)

    ###################
    # WZ region
    ##################

    #WZ region lepton number selections
    WZ_nl=False
    #WZ region b-jet selection
    WZ_nb=False
    #WZ region lepton kinematics selctions
    WZ_leptons=False
    #WZ region MET selection
    WZ_MET=False
    #WZ region tag, 0: fail to pass the WZ selection, 1:3 muon, 2:2muon, 3:1muon, 4:0 muon
    WZ_region=0
    WZ_zl1_id=-1
    WZ_zl2_id=-1
    WZ_wl_id=-1
    WZ_zl1_pdgid=-99
    WZ_zl2_pdgid=-99
    WZ_wl_pdgid=-99
    WZ_zl1_pt=-99
    WZ_zl1_eta=-99
    WZ_zl1_phi=-99
    WZ_zl1_mass=-99
    WZ_zl2_pt=-99
    WZ_zl2_eta=-99
    WZ_zl2_phi=-99
    WZ_zl2_mass=-99
    WZ_l3_pt=-99
    WZ_l3_eta=-99
    WZ_l3_phi=-99
    WZ_l3_mass=-99
    WZ_Z_mass=-99
    WZ_Z_pt=-99
    WZ_Z_eta=-99
    WZ_Z_phi=-99
    WZ_met=-99
    # the first two leading leptons with pt >20, 3rd lepton pt >15, 4th lepton veto
    if len(tightLeptons)==3 and tightLeptons[1].Pt()>20 and len(looseLeptons)==0:
      WZ_nl=True
    # no bjet
    if WZ_nl and tightJets_b_DeepCSVmedium_id[0]==-1:
      WZ_nb=True

    # mll>4 regardless the flavor and charge sign
    if WZ_nb and (tightLeptons[0]+tightLeptons[1]).M()>4 and (tightLeptons[2]+tightLeptons[1]).M()>4 and (tightLeptons[0]+tightLeptons[2]).M()>4:
      WZ_leptons=True
    
    if WZ_leptons and ((self.is_mc and event.MET_T1Smear_pt>30) or (event.MET_T1_pt>30 and (not (self.is_mc)))):
      WZ_MET=True

    if WZ_MET:
      if self.is_mc:
        WZ_met=event.MET_T1Smear_pt
      else:
        WZ_met=event.MET_T1_pt
      # 3 muons case
      if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1]+tightMuons_pdgid[2])==13:
	#two combination 0+2 or 1+2
        if (tightMuons_pdgid[0]-tightMuons_pdgid[1])==0:
          if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[1]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[1].Pt()
	    WZ_l3_eta=tightMuons[1].Eta()
	    WZ_l3_phi=tightMuons[1].Phi()
	    WZ_l3_mass=tightMuons[1].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[2]).Phi()

	  if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
            WZ_zl1_pdgid=tightMuons_pdgid[1]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[0]
	    WZ_zl1_pt=tightMuons[1].Pt()
	    WZ_zl1_eta=tightMuons[1].Eta()
	    WZ_zl1_phi=tightMuons[1].Phi()
	    WZ_zl1_mass=tightMuons[1].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[0].Pt()
	    WZ_l3_eta=tightMuons[0].Eta()
	    WZ_l3_phi=tightMuons[0].Phi()
	    WZ_l3_mass=tightMuons[0].M()
	    WZ_Z_mass=(tightMuons[1]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[1]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[1]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[1]+tightMuons[2]).Phi()
	#two combination 0+1 or 1+2
	elif (tightMuons_pdgid[0]-tightMuons_pdgid[2])==0:
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[1]
            WZ_wl_id=tightMuons_id[2]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[1]
            WZ_wl_pdgid=tightMuons_pdgid[2]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[1].Pt()
	    WZ_zl2_eta=tightMuons[1].Eta()
	    WZ_zl2_phi=tightMuons[1].Phi()
	    WZ_zl2_mass=tightMuons[1].M()
	    WZ_l3_pt=tightMuons[2].Pt()
	    WZ_l3_eta=tightMuons[2].Eta()
	    WZ_l3_phi=tightMuons[2].Phi()
	    WZ_l3_mass=tightMuons[2].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()

	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
            WZ_zl1_pdgid=tightMuons_pdgid[1]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[0]
	    WZ_zl1_pt=tightMuons[1].Pt()
	    WZ_zl1_eta=tightMuons[1].Eta()
	    WZ_zl1_phi=tightMuons[1].Phi()
	    WZ_zl1_mass=tightMuons[1].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[0].Pt()
	    WZ_l3_eta=tightMuons[0].Eta()
	    WZ_l3_phi=tightMuons[0].Phi()
	    WZ_l3_mass=tightMuons[0].M()
	    WZ_Z_mass=(tightMuons[1]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[1]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[1]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[1]+tightMuons[2]).Phi()
	#two combination 0+1 or 0+2
	else:
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[1]
            WZ_wl_id=tightMuons_id[2]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[1]
            WZ_wl_pdgid=tightMuons_pdgid[2]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[1].Pt()
	    WZ_zl2_eta=tightMuons[1].Eta()
	    WZ_zl2_phi=tightMuons[1].Phi()
	    WZ_zl2_mass=tightMuons[1].M()
	    WZ_l3_pt=tightMuons[2].Pt()
	    WZ_l3_eta=tightMuons[2].Eta()
	    WZ_l3_phi=tightMuons[2].Phi()
	    WZ_l3_mass=tightMuons[2].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[1]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[1].Pt()
	    WZ_l3_eta=tightMuons[1].Eta()
	    WZ_l3_phi=tightMuons[1].Phi()
	    WZ_l3_mass=tightMuons[1].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[2]).Phi()

      # 2 muons case
      if len(tightElectrons)==1 and (tightMuons_pdgid[0]-tightMuons_pdgid[1])==0:
	if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	  WZ_region=2
	  WZ_zl1_id=tightMuons_id[0]
	  WZ_zl2_id=tightMuons_id[1]
	  WZ_wl_id=tightElectrons_id[0]
	  WZ_zl1_pdgid=tightMuons_pdgid[0]
	  WZ_zl2_pdgid=tightMuons_pdgid[1]
	  WZ_wl_pdgid=tightElectrons_pdgid[0]
	  WZ_zl1_pt=tightMuons[0].Pt()
	  WZ_zl1_eta=tightMuons[0].Eta()
	  WZ_zl1_phi=tightMuons[0].Phi()
	  WZ_zl1_mass=tightMuons[0].M()
	  WZ_zl2_pt=tightMuons[1].Pt()
	  WZ_zl2_eta=tightMuons[1].Eta()
	  WZ_zl2_phi=tightMuons[1].Phi()
	  WZ_zl2_mass=tightMuons[1].M()
	  WZ_l3_pt=tightElectrons[0].Pt()
	  WZ_l3_eta=tightElectrons[0].Eta()
	  WZ_l3_phi=tightElectrons[0].Phi()
	  WZ_l3_mass=tightElectrons[0].M()
	  WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	  WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	  WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	  WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()

      # 1 muon case
      if len(tightElectrons)==2 and (tightElectrons_pdgid[0]-tightElectrons_pdgid[1])==0:
	if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	  WZ_region=3
	  WZ_zl1_id=tightElectrons_id[0]
	  WZ_zl2_id=tightElectrons_id[1]
	  WZ_wl_id=tightMuons_id[0]
	  WZ_zl1_pdgid=tightElectrons_pdgid[0]
	  WZ_zl2_pdgid=tightElectrons_pdgid[1]
	  WZ_wl_pdgid=tightMuons_pdgid[0]
	  WZ_zl1_pt=tightElectrons[0].Pt()
	  WZ_zl1_eta=tightElectrons[0].Eta()
	  WZ_zl1_phi=tightElectrons[0].Phi()
	  WZ_zl1_mass=tightElectrons[0].M()
	  WZ_zl2_pt=tightElectrons[1].Pt()
	  WZ_zl2_eta=tightElectrons[1].Eta()
	  WZ_zl2_phi=tightElectrons[1].Phi()
	  WZ_zl2_mass=tightElectrons[1].M()
	  WZ_l3_pt=tightMuons[0].Pt()
	  WZ_l3_eta=tightMuons[0].Eta()
	  WZ_l3_phi=tightMuons[0].Phi()
	  WZ_l3_mass=tightMuons[0].M()
	  WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	  WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	  WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	  WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()

      # 0 muon case
      if len(tightElectrons)==3 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1]+tightElectrons_pdgid[2])==11:
	#two combination 0+2 or 1+2
        if (tightElectrons_pdgid[0]-tightElectrons_pdgid[1])==0:
          if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[1]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[1].Pt()
	    WZ_l3_eta=tightElectrons[1].Eta()
	    WZ_l3_phi=tightElectrons[1].Phi()
	    WZ_l3_mass=tightElectrons[1].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[2]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
            WZ_zl1_pdgid=tightElectrons_pdgid[1]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[0]
	    WZ_zl1_pt=tightElectrons[1].Pt()
	    WZ_zl1_eta=tightElectrons[1].Eta()
	    WZ_zl1_phi=tightElectrons[1].Phi()
	    WZ_zl1_mass=tightElectrons[1].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[0].Pt()
	    WZ_l3_eta=tightElectrons[0].Eta()
	    WZ_l3_phi=tightElectrons[0].Phi()
	    WZ_l3_mass=tightElectrons[0].M()
	    WZ_Z_mass=(tightElectrons[1]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[1]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[1]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[1]+tightElectrons[2]).Phi()
	#two combination 0+1 or 1+2
	elif (tightElectrons_pdgid[0]-tightElectrons_pdgid[2])==0:
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[1]
            WZ_wl_id=tightElectrons_id[2]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[1]
            WZ_wl_pdgid=tightElectrons_pdgid[2]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[1].Pt()
	    WZ_zl2_eta=tightElectrons[1].Eta()
	    WZ_zl2_phi=tightElectrons[1].Phi()
	    WZ_zl2_mass=tightElectrons[1].M()
	    WZ_l3_pt=tightElectrons[2].Pt()
	    WZ_l3_eta=tightElectrons[2].Eta()
	    WZ_l3_phi=tightElectrons[2].Phi()
	    WZ_l3_mass=tightElectrons[2].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
            WZ_zl1_pdgid=tightElectrons_pdgid[1]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[0]
	    WZ_zl1_pt=tightElectrons[1].Pt()
	    WZ_zl1_eta=tightElectrons[1].Eta()
	    WZ_zl1_phi=tightElectrons[1].Phi()
	    WZ_zl1_mass=tightElectrons[1].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[0].Pt()
	    WZ_l3_eta=tightElectrons[0].Eta()
	    WZ_l3_phi=tightElectrons[0].Phi()
	    WZ_l3_mass=tightElectrons[0].M()
	    WZ_Z_mass=(tightElectrons[1]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[1]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[1]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[1]+tightElectrons[2]).Phi()
	#two combination 0+1 or 0+2
	else:
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[1]
            WZ_wl_id=tightElectrons_id[2]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[1]
            WZ_wl_pdgid=tightElectrons_pdgid[2]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[1].Pt()
	    WZ_zl2_eta=tightElectrons[1].Eta()
	    WZ_zl2_phi=tightElectrons[1].Phi()
	    WZ_zl2_mass=tightElectrons[1].M()
	    WZ_l3_pt=tightElectrons[2].Pt()
	    WZ_l3_eta=tightElectrons[2].Eta()
	    WZ_l3_phi=tightElectrons[2].Phi()
	    WZ_l3_mass=tightElectrons[2].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[1]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[1].Pt()
	    WZ_l3_eta=tightElectrons[1].Eta()
	    WZ_l3_phi=tightElectrons[1].Phi()
	    WZ_l3_mass=tightElectrons[1].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[2]).Phi()
	
    self.out.fillBranch("WZ_region", WZ_region)
    self.out.fillBranch("WZ_zl1_id", WZ_zl1_id)
    self.out.fillBranch("WZ_zl2_id", WZ_zl2_id)
    self.out.fillBranch("WZ_wl_id", WZ_wl_id)
    self.out.fillBranch("WZ_zl1_pdgid", WZ_zl1_pdgid)
    self.out.fillBranch("WZ_zl2_pdgid", WZ_zl2_pdgid)
    self.out.fillBranch("WZ_wl_pdgid", WZ_wl_pdgid)
    self.out.fillBranch("WZ_zl1_pt", WZ_zl1_pt)
    self.out.fillBranch("WZ_zl1_eta", WZ_zl1_eta)
    self.out.fillBranch("WZ_zl1_phi", WZ_zl1_phi)
    self.out.fillBranch("WZ_zl1_mass", WZ_zl1_mass)
    self.out.fillBranch("WZ_zl2_pt", WZ_zl2_pt)
    self.out.fillBranch("WZ_zl2_eta", WZ_zl2_eta)
    self.out.fillBranch("WZ_zl2_phi", WZ_zl2_phi)
    self.out.fillBranch("WZ_zl2_mass", WZ_zl2_mass)
    self.out.fillBranch("WZ_l3_pt", WZ_l3_pt)
    self.out.fillBranch("WZ_l3_eta", WZ_l3_eta)
    self.out.fillBranch("WZ_l3_phi", WZ_l3_phi)
    self.out.fillBranch("WZ_l3_mass", WZ_l3_mass)
    self.out.fillBranch("WZ_Z_mass", WZ_Z_mass)
    self.out.fillBranch("WZ_Z_pt", WZ_Z_pt)
    self.out.fillBranch("WZ_Z_eta", WZ_Z_eta)
    self.out.fillBranch("WZ_Z_phi", WZ_Z_phi)
    self.out.fillBranch("WZ_met", WZ_met)
    

    #  DDDD   YY      YY  (opposite sign) region: two opposite sign lepton, with |mll-91.1876|<15
    #  D   D    YY  YY   
    #  D    D     YY
    #  D   D      YY
    #  DDDD       YY

    #DY region lepton number selections (ttbar region are similar)
    DY_nl=False
    #DY region tag, 0: fail to pass the DY selection, 1:2 muon, 2:1 muon, 3:0 muon
    DY_region=0
    DY_l1_id=-1
    DY_l2_id=-1
    DY_l1_pdgid=-99
    DY_l2_pdgid=-99
    DY_l1_pt=-99
    DY_l1_eta=-99
    DY_l1_phi=-99
    DY_l1_mass=-99
    DY_l2_pt=-99
    DY_l2_eta=-99
    DY_l2_phi=-99
    DY_l2_mass=-99
    DY_z_mass=-99
    DY_z_pt=-99
    DY_z_eta=-99
    DY_z_phi=-99
    DY_drll=-99

    # the two leptons with pt >20, 3th lepton veto
    if len(tightLeptons)==2 and tightLeptons[1].Pt()>20 and len(looseLeptons)==0:
      DY_nl=True
    if DY_nl:
      # 2 muons case
      if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1])==0:
	DY_region=1
	DY_l1_id=tightMuons_id[0]
	DY_l2_id=tightMuons_id[1]
	DY_l1_pdgid=tightMuons_pdgid[0]
	DY_l2_pdgid=tightMuons_pdgid[1]
	DY_l1_pt=tightMuons[0].Pt()
	DY_l1_eta=tightMuons[0].Eta()
	DY_l1_phi=tightMuons[0].Phi()
	DY_l1_mass=tightMuons[0].M()
	DY_l2_pt=tightMuons[1].Pt()
	DY_l2_eta=tightMuons[1].Eta()
	DY_l2_phi=tightMuons[1].Phi()
	DY_l2_mass=tightMuons[1].M()
	DY_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	DY_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	DY_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	DY_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	DY_drll=tightLeptons[0].DeltaR(tightLeptons[1])
      # 2 eles case
      if len(tightElectrons)==2 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1])==0:
	DY_region=3
        DY_l1_id=tightElectrons_id[0]
        DY_l2_id=tightElectrons_id[1]
        DY_l1_pdgid=tightElectrons_pdgid[0]
        DY_l2_pdgid=tightElectrons_pdgid[1]
	DY_l1_pt=tightElectrons[0].Pt()
	DY_l1_eta=tightElectrons[0].Eta()
	DY_l1_phi=tightElectrons[0].Phi()
	DY_l1_mass=tightElectrons[0].M()
	DY_l2_pt=tightElectrons[1].Pt()
	DY_l2_eta=tightElectrons[1].Eta()
	DY_l2_phi=tightElectrons[1].Phi()
	DY_l2_mass=tightElectrons[1].M()
	DY_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	DY_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	DY_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	DY_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	DY_drll=tightLeptons[0].DeltaR(tightLeptons[1])
      # 1 ele case
      if len(tightElectrons)==1 and (sign(tightMuons_pdgid[0])+sign(tightElectrons_pdgid[0]))==0:
	DY_region=2
        DY_l1_id=tightMuons_id[0]
        DY_l2_id=tightElectrons_id[0]
        DY_l1_pdgid=tightMuons_pdgid[0]
        DY_l2_pdgid=tightElectrons_pdgid[0]
	DY_l1_pt=tightMuons[0].Pt()
	DY_l1_eta=tightMuons[0].Eta()
	DY_l1_phi=tightMuons[0].Phi()
	DY_l1_mass=tightMuons[0].M()
	DY_l2_pt=tightElectrons[0].Pt()
	DY_l2_eta=tightElectrons[0].Eta()
	DY_l2_phi=tightElectrons[0].Phi()
	DY_l2_mass=tightElectrons[0].M()
	DY_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	DY_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	DY_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	DY_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	DY_drll=tightLeptons[0].DeltaR(tightLeptons[1])

    self.out.fillBranch("DY_region", DY_region)
    self.out.fillBranch("DY_l1_id", DY_l1_id)
    self.out.fillBranch("DY_l2_id", DY_l2_id)
    self.out.fillBranch("DY_l1_pdgid", DY_l1_pdgid)
    self.out.fillBranch("DY_l2_pdgid", DY_l2_pdgid)
    self.out.fillBranch("DY_l1_pt", DY_l1_pt)
    self.out.fillBranch("DY_l1_eta", DY_l1_eta)
    self.out.fillBranch("DY_l1_phi", DY_l1_phi)
    self.out.fillBranch("DY_l1_mass", DY_l1_mass)
    self.out.fillBranch("DY_l2_pt", DY_l2_pt)
    self.out.fillBranch("DY_l2_eta", DY_l2_eta)
    self.out.fillBranch("DY_l2_phi", DY_l2_phi)
    self.out.fillBranch("DY_l2_mass", DY_l2_mass)
    self.out.fillBranch("DY_z_mass", DY_z_mass)
    self.out.fillBranch("DY_z_pt", DY_z_pt)
    self.out.fillBranch("DY_z_eta", DY_z_eta)
    self.out.fillBranch("DY_z_phi", DY_z_phi)
    self.out.fillBranch("DY_drll", DY_drll)

    if not (bh_nl or WZ_region >0 or DY_region>0):
      return False

    return True

BH2016 = lambda: BHProducer("2016")
BH2017 = lambda: BHProducer("2017")
BH2018 = lambda: BHProducer("2018")
