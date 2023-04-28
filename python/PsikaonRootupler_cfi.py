import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKaon',
                          #dimuons = cms.InputTag("selectedPatMuons"),
                          dimuons = cms.InputTag("slimmedMuons"),
                          #Trak = cms.InputTag("cleanPatTrackCands"),
                          #Trak = cms.InputTag("packedPFCandidates"),
                          #Trak_lowpt = cms.InputTag("lostTracks"),
                          #primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerInput = cms.InputTag("slimmedPatTrigger"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          algInputTag = cms.InputTag("gtStage2Digis", "", "RECO"),
                          l1Muons = cms.InputTag("gmtStage2Digis", "Muon", "RECO"), 
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          )
