import FWCore.ParameterSet.Config as cms

process = cms.Process("DiMuonDiTauTreelizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuDiTauSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_6131.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataB_17Nov2017_v2/190327_111844/0000/mumutautau_12.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataB_17Nov2017_v2/190327_111844/0000/mumutautau_15.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataB_17Nov2017_v2/190327_111844/0000/mumutautau_16.root',
    )
)

######### embed 2017v2 tauID into the miniAOD ###############
# reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Rerunning_of_the_tau_ID_on_M_AN1

from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA import *
myTool = TauIDEmbedder(process, cms,
        debug = True,
        toKeep = ["2017v2"]
)
myTool.runTauID()

############################################################

process.treelizer = cms.Sequence(
        process.lumiTree*
        process.HLTEle*
        process.MuonPtEtaCut*
        process.MuonID*
        process.LeadingMuonIso*
        process.TrigMuMatcher*
        process.SecondThirdMuonSelector*
        process.rerunMvaIsolationSequence*
        process.NewTauIDsEmbedded*
        process.TauHadSelector*
        process.MuMuTauTauAnalyzer
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
        fileName =  cms.string('MuMuTauTauTreelization.root')
)

process.p = cms.Path(process.treelizer)
