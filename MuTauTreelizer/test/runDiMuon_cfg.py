import FWCore.ParameterSet.Config as cms

process = cms.Process("DiMuonTreelizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_zmutauskim_6131.root'
    )
)

process.treelizer = cms.Sequence(
        process.lumiTree*
        process.HLTEle*
        process.MuonPtEtaCut*
        process.MuonID*
        process.LeadingMuonIso*
        process.TrigMuMatcher*
        process.SecondMuonSelector*
        process.DiMuonMassSelector*
        process.DiMuonAnalyzer
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
        fileName =  cms.string('MuMuTreelization.root')
)

process.p = cms.Path(process.treelizer)
