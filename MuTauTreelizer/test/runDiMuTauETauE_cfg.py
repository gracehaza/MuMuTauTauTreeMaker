import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/MiniAOD_H125AA19_DiMuDiTau_Fall17DRPremix_v1/190515_140053/0000/mumutautau_1.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.parseArguments()

process = cms.Process("DiMuonTauETauETreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runDiMuTauETauE_cfg.py isMC=1 ###############
##########################################################################

if options.isMC == 1:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuTauETauESelectorMC_cfi")

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuTauETauESelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# --- please specify the sample that you need to run for local test ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

############################################################

############################################################
process.treelizer = cms.Sequence(
        process.lumiTree*
        process.HLTEle*
        process.TrigMuMatcher*
        process.MuonPtEtaCut*
        process.MuonID*
        process.LeadingMuonIso*
        process.SecondThirdMuonSelector*
        process.DiMuonMassSelector*
        process.ElectronSelector*
        process.JetSelector*
        process.MuMuTauETauEAnalyzer
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
        fileName =  cms.string('MuMuTauETauETreelization.root')
)

process.p = cms.Path(process.treelizer)