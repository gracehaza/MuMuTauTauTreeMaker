import FWCore.ParameterSet.Config as cms

lumiTree = cms.EDAnalyzer("LumiTree",
        genEventInfo = cms.InputTag("generator"),
        nevents = cms.InputTag('lumiSummary','numberOfEvents'),
        summedWeights = cms.InputTag('lumiSummary','sumOfWeightedEvents'),
)

HLTEle = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
        HLTPaths = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*","HLT_IsoMu27_v*","HLT_IsoTkMu27_v*"),
        eventSetupPathsKey = cms.string(''),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False), # throw exception on unknown path names
)

MuonID = cms.EDFilter("MuonID",
        muonTag = cms.InputTag('slimmedMuons'),
        muonID = cms.string('medium'),
        minNumObjsToPassFilter = cms.int32(1),
)

MuonSelector = cms.EDFilter("MuonSelector",
        muonTag = cms.InputTag("MuonID"),
        relIsoCutVal = cms.double(-1), # positive number for iso threshold, -1 for ignoring iso
        normalRelIso = cms.bool(True), #True = Iso-mu; False = inverted Iso-mu
        Eta = cms.double(2.5),
        Pt = cms.double(3.0),
)

TrigMuMatcher = cms.EDFilter("TrigMuMatcher",
        muonsTag = cms.InputTag('MuonSelector'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v","HLT_IsoMu27_v","HLT_IsoTkMu27_v"),
        dRCut = cms.double(0.15),
        muPtCut = cms.double(26.0),
)

TauHadSelector = cms.EDFilter("TauHadSelector",
        tauTag = cms.InputTag('slimmedTausNewID'),
        tauDiscriminatorTag = cms.vstring('decayModeFindingNewDMs'),
        passDiscriminator = cms.bool(True),
        pTMin = cms.double(8.0),
        etaMax = cms.double(2.4),
)

JetSelector = cms.EDFilter("JetSelector",
        jetTag = cms.InputTag('slimmedJets'),
        jetIdName = cms.string("Tight"), # reference: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
        etaCut = cms.double(2.4),
        ptCut = cms.double(20),
)

ZTauMuTauHadAnalyzer = cms.EDAnalyzer('ZTauMuTauHadAnalyzer',
        MuTag = cms.InputTag("TrigMuMatcher"),
        TauTag = cms.InputTag("TauHadSelector"),
        JetTag = cms.InputTag("JetSelector"),
        MetTag = cms.InputTag("slimmedMETs"),
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        isMC = cms.bool(False),
)
