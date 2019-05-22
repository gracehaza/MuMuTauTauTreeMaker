import FWCore.ParameterSet.Config as cms

lumiTree = cms.EDAnalyzer("LumiTree",
        genEventInfo = cms.InputTag("generator"),
        nevents = cms.InputTag('lumiSummary','numberOfEvents'),
        summedWeights = cms.InputTag('lumiSummary','sumOfWeightedEvents')
)

HLTEle = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
        HLTPaths = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"),
        eventSetupPathsKey = cms.string(''),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False) # throw exception on unknown path names
)

MuonPtEtaCut = cms.EDFilter("MuonPtEtaCut",
        muonTag = cms.InputTag("slimmedMuons"),
        Eta = cms.double(2.4),
        Pt = cms.double(3.0),
        minNumObjsToPassFilter = cms.uint32(2)
)

MuonID = cms.EDFilter("MuonID",
        muonTag = cms.InputTag('MuonPtEtaCut'),
        muonID = cms.string('loose'),
        minNumObjsToPassFilter = cms.int32(2)
)

LeadingMuonIso = cms.EDFilter("LeadingMuonIso",
        muonTag = cms.InputTag('MuonID'),
        relIsoCutVal = cms.double(0.25), # 0.25 for iso, -1 for ignoring iso
        passRelIso = cms.bool(True) #False = Non-Iso DiMu, True = Iso-DiMu
)

TrigMuMatcher = cms.EDFilter("TrigMuMatcher",
        muonsTag = cms.InputTag('LeadingMuonIso'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v"),
        dRCut = cms.double(0.1),
        mu1PtCut = cms.double(26.0)
)

SecondMuonSelector = cms.EDFilter("SecondMuonSelector",
        muonTag = cms.InputTag('MuonID'),
        mu1Tag = cms.InputTag('TrigMuMatcher'),
        relIsoCutVal = cms.double(0.25), # .25 for iso, -1 for ignoring iso
        passRelIso = cms.bool(True), #False = Non-Iso DiMu, True = Iso-DiMu
        dRCut = cms.double(-1), # -1 = no dR cut, >0 for dR low threshold
        oppositeSign = cms.bool(True), # False for SameSignDiMu, True regular
)

DiMuonMassSelector = cms.EDFilter("DiMuonMassSelector",
        mu1Tag = cms.InputTag('TrigMuMatcher'),
        mu2Tag = cms.InputTag('SecondMuonSelector'),
        minMass = cms.double(30),
        maxMass = cms.double(200)
)

DiMuonAnalyzer = cms.EDAnalyzer('DiMuonAnalyzer',
        Mu1Mu2 = cms.InputTag("DiMuonMassSelector"),
        Vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
        isMC = cms.bool(False),
        Generator = cms.InputTag("generator")
)
