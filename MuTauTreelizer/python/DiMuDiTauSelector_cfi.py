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
        minNumObjsToPassFilter = cms.uint32(3)
)

MuonID = cms.EDFilter("MuonID",
        muonTag = cms.InputTag('MuonPtEtaCut'),
        muonID = cms.string('loose'),
        minNumObjsToPassFilter = cms.int32(3)
)

LeadingMuonIso = cms.EDFilter("LeadingMuonIso",
        muonTag = cms.InputTag('MuonID'),
        relIsoCutVal = cms.double(-1), # 0.25 for iso, -1 for ignoring iso
        passRelIso = cms.bool(False) #False = Non-Iso DiMu, True = Iso-DiMu
)

TrigMuMatcher = cms.EDFilter("TrigMuMatcher",
        muonsTag = cms.InputTag('LeadingMuonIso'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v"),
        dRCut = cms.double(0.1),
        mu1PtCut = cms.double(26.0)
)

SecondThirdMuonSelector = cms.EDFilter("SecondThirdMuonSelector",
        muonTag = cms.InputTag('MuonID'),
        mu1Tag = cms.InputTag('TrigMuMatcher'),
        relIsoCutVal = cms.double(-1), # .25 for iso, -1 for ignoring iso
        passRelIso = cms.bool(False), #False = Non-Iso DiMu, True = Iso-DiMu
        oppositeSign = cms.bool(True), # False for SameSignTriMuon, True for at least one pair of oppositely charged muons
)

TauHadSelector = cms.EDFilter("TauHadSelector",
        tauTag = cms.InputTag('NewTauIDsEmbedded'), # output of configuration: "TauIdMVA.py"
        #tauTag = cms.InputTag('selectedPatTausMuonCleaned'),
        tauDiscriminatorTag = cms.vstring('decayModeFinding'),
        passDiscriminator = cms.bool(True),
        pTMin = cms.double(8.0),
        etaMax = cms.double(2.4),
)

MuMuTauTauAnalyzer = cms.EDAnalyzer('MuMuTauTauAnalyzer',
        Mu1Tag = cms.InputTag("TrigMuMatcher"),
        Mu2Mu3Tag = cms.InputTag("SecondThirdMuonSelector"),
        TauHadTag = cms.InputTag("TauHadSelector"),
        Vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
        tauDiscriminatorTags = cms.vstring('byIsolationMVArun2017v2DBoldDMwLTraw2017',
            'byVVLooseIsolationMVArun2017v2DBoldDMwLT2017',
            'byVLooseIsolationMVArun2017v2DBoldDMwLT2017',
            'byLooseIsolationMVArun2017v2DBoldDMwLT2017',
            'byMediumIsolationMVArun2017v2DBoldDMwLT2017',
            'byTightIsolationMVArun2017v2DBoldDMwLT2017',
            'byVTightIsolationMVArun2017v2DBoldDMwLT2017',
            'byVVTightIsolationMVArun2017v2DBoldDMwLT2017',
        ),
        isMC = cms.bool(False),
        Generator = cms.InputTag("generator"),
)
