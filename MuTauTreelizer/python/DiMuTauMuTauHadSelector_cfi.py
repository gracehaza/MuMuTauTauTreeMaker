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

TrigMuMatcher = cms.EDFilter("TrigMuMatcher",
        muonsTag = cms.InputTag('slimmedMuons'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v","HLT_IsoMu27_v*","HLT_IsoTkMu27_v*"),
        dRCut = cms.double(0.15),
        muPtCut = cms.double(26.0),
)

MuonPtEtaCut = cms.EDFilter("MuonPtEtaCut",
        muonTag = cms.InputTag("TrigMuMatcher"),
        Eta = cms.double(2.4),
        Pt = cms.double(3.0),
        minNumObjsToPassFilter = cms.uint32(3),
)

MuonID = cms.EDFilter("MuonID",
        muonTag = cms.InputTag('MuonPtEtaCut'),
        muonID = cms.string('loose'),
        minNumObjsToPassFilter = cms.int32(3),
)

LeadingMuonIso = cms.EDFilter("LeadingMuonIso",
        muonTag = cms.InputTag('MuonID'),
        relIsoCutVal = cms.double(-1), # 0.25 for iso, -1 for ignoring iso
        passRelIso = cms.bool(False), #False = Non-Iso DiMu, True = Iso-DiMu
)

SecondMuonSelector = cms.EDFilter("SecondMuonSelector",
        muonTag = cms.InputTag('MuonID'),
        mu1Tag = cms.InputTag('LeadingMuonIso'),
        relIsoCutVal = cms.double(-1), # .25 for iso, -1 for ignoring iso
        passRelIso = cms.bool(False), #False = Non-Iso DiMu, True = Iso-DiMu
        dRCut = cms.double(-1), # -1 = no dR cut, >0 for dR low threshold
        oppositeSign = cms.bool(True), # False for SameSignDiMu, True regular
)

DiMuonMassSelector = cms.EDFilter("DiMuonMassSelector",
        mu1Tag = cms.InputTag('LeadingMuonIso'),
        mu2Tag = cms.InputTag('SecondMuonSelector'),
        minMass = cms.double(3),
        maxMass = cms.double(400),
)

ThirdMuonSelector = cms.EDFilter("ThirdMuonSelector",
        muonTag = cms.InputTag('MuonID'),
        mu1mu2Tag = cms.InputTag('DiMuonMassSelector'),
        dRCut = cms.double(-1), # -1 = no dR cut, >0 for dR(of mu3 from mu1 and mu2) low threshold
)

TauHadSelector = cms.EDFilter("TauHadSelector",
        tauTag = cms.InputTag('NewTauIDsEmbedded'), # output of configuration: "TauIdMVA.py"
        tauDiscriminatorTag = cms.vstring('decayModeFinding'),
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

MuMuTauMuTauHadAnalyzer = cms.EDAnalyzer('MuMuTauMuTauHadAnalyzer',
        Mu1Mu2Tag = cms.InputTag("DiMuonMassSelector"),
        Mu3Tag = cms.InputTag("ThirdMuonSelector"),
        TauTag = cms.InputTag("TauHadSelector"),
        JetTag = cms.InputTag("JetSelector"),
        MetTag = cms.InputTag("slimmedMETs"),
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        isMC = cms.bool(False),
)
