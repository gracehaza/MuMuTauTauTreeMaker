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

HLTEleDiMu = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
        HLTPaths = cms.vstring("HLT_Dimuon20_Jpsi_v*", "HLT_Dimuon25_Jpsi_v*","HLT_DoubleMu4_JpsiTrk_Displaced_v*"),
        eventSetupPathsKey = cms.string(''),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False), # throw exception on unknown path names
)


MuonID = cms.EDFilter("MuonID",
        muonTag = cms.InputTag('slimmedMuons'),
        muonID = cms.string('loose'),
        minNumObjsToPassFilter = cms.int32(2),
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

ElectronCandSelector = cms.EDFilter("ElectronCandSelector",
        electronTag = cms.InputTag('slimmedElectrons'),
        # --- need the two parameters below for electron isolation computation ---
        rhoTag = cms.InputTag("fixedGridRhoAll"),
        effAreasConfigFile = cms.FileInPath("MuMuTauTauTreeMaker/MuTauTreelizer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
        # ========================================================================
        relIdName = cms.string("Loose"), # customize electron ID options: Veto, Loose, Medium, Tight
        passRelIso = cms.bool(False),
        etaCut = cms.double(2.5),
        ptCut = cms.double(7.0),
)

TauCandSelector = cms.EDFilter("TauCandSelector",
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

#postfix=''
DeepDiTauProducer = cms.EDProducer("DeepDiTauProducer",
#       src = cms.InputTag('slimmedJets'),
        slimmedJetTag = cms.InputTag('slimmedJets'),
        DeepDiTauConfiguration = cms.PSet(
            memmapped = cms.bool(False),
            graphDefinitions = cms.VPSet(
                cms.PSet(
                    name = cms.string('ditau2017v1'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/03Dec2020_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/03Dec2020_mean_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('ditau2017MDv1'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/17Feb2021_relu_tanh_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/17Feb2021_relu_tanh_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton_charm')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton_charm_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton_charm_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton_charm')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_means_sigmas.txt'),
                ),
            ),
        ),
    )

JetIdEmbedder = cms.EDProducer("JetIdEmbedder",
        slimmedJetTag = cms.InputTag('slimmedJets'),
        discriminator = cms.string('pileupJetId:fullDiscriminant'),
        ditau2017v1 = cms.InputTag("DeepDiTauProducer","ditau2017v1"),
        ditau2017MDv1 = cms.InputTag("DeepDiTauProducer","ditau2017MDv1"),
        DeepDiTau_boosted = cms.InputTag("DeepDiTauProducer","DeepDiTau_boosted"),
        DeepDiTau_boosted_massdeco = cms.InputTag("DeepDiTauProducer","DeepDiTau_boosted_massdeco"),
        DeepDiTau_boosted_nolepton_charm = cms.InputTag("DeepDiTauProducer","DeepDiTau_boosted_nolepton_charm"),
        DeepDiTau_boosted_nolepton_charm_massdeco = cms.InputTag("DeepDiTauProducer","DeepDiTau_boosted_nolepton_charm_massdeco"),
        DeepDiTau_boosted_nolepton_massdeco= cms.InputTag("DeepDiTauProducer","DeepDiTau_boosted_nolepton_massdeco"),
        DeepDiTau_boosted_nolepton = cms.InputTag("DeepDiTauProducer","DeepDiTau_boosted_nolepton"),
        DeepDiTau_massdeco = cms.InputTag("DeepDiTauProducer","DeepDiTau_massdeco"),
        DeepDiTau = cms.InputTag("DeepDiTauProducer","DeepDiTau"),
        DeepDiTau_nolepton_charm_massdeco = cms.InputTag("DeepDiTauProducer","DeepDiTau_nolepton_charm_massdeco"),
        DeepDiTau_nolepton_charm = cms.InputTag("DeepDiTauProducer","DeepDiTau_nolepton_charm"),
        DeepDiTau_nolepton_massdeco = cms.InputTag("DeepDiTauProducer","DeepDiTau_nolepton_massdeco"),
        DeepDiTau_nolepton = cms.InputTag("DeepDiTauProducer","DeepDiTau_nolepton"),                       
   )

GenMuonCandSelector = cms.EDFilter("GenMuonCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenElectronCandSelector = cms.EDFilter("GenElectronCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenTauMuCandSelector = cms.EDFilter("GenTauMuCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenTauEleCandSelector = cms.EDFilter("GenTauEleCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenTauHadCandSelector = cms.EDFilter("GenTauHadCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

DiMuDiTauAnalyzer = cms.EDAnalyzer('DiMuDiTauAnalyzer',
        GenParticleTag = cms.InputTag('prunedGenParticles'),
        MuTag = cms.InputTag("TrigMuMatcher"),
        EleTag = cms.InputTag("ElectronCandSelector"),
        TauTag = cms.InputTag("TauCandSelector"),
        JetTag = cms.InputTag("JetSelector"),
        MetTag = cms.InputTag("slimmedMETs"),
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        rhoTag = cms.InputTag("fixedGridRhoAll"),
        effAreasConfigFile = cms.FileInPath("MuMuTauTauTreeMaker/MuTauTreelizer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
        isMC = cms.bool(True),
        slimmedJetTag = cms.InputTag('JetIdEmbedder'),        
        GenMuTag = cms.InputTag('GenMuonCandSelector'),
        GenEleTag = cms.InputTag('GenElectronCandSelector'),
        GenTauMuTag = cms.InputTag('GenTauMuCandSelector'),
        GenTauEleTag = cms.InputTag('GenTauEleCandSelector'),
        GenTauHadTag = cms.InputTag('GenTauHadCandSelector'),
        PileupTag = cms.InputTag("slimmedAddPileupInfo"),
        Generator = cms.InputTag("generator"),
)
