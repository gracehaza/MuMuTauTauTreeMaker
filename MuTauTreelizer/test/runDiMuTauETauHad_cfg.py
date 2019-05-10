import FWCore.ParameterSet.Config as cms

process = cms.Process("DiMuonTauETauHadTreelizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuTauETauHadSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_6131.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0000/mumutautau_8.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0000/mumutautau_80.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0000/mumutautau_801.root',
    )
)

######### embed 2017v2 tauID into the miniAOD ###############
# reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Rerunning_of_the_tau_ID_on_M_AN1

from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVAForElectronCleaned import *
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
        process.SecondMuonSelector*
        process.DiMuonMassSelector*
        process.ElectronSelector*
        process.rerunMvaIsolationSequence*
        process.NewTauIDsEmbedded*
        process.TauHadSelector
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
        fileName =  cms.string('MuMuTauTauTreelization.root')
)

process.p = cms.Path(process.treelizer)
