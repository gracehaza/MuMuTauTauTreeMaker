import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/MiniAOD_H125AA19_DiMuDiTau_Fall17DRPremix_v1/190515_140053/0000/mumutautau_1.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.register('tauCluster', 3, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "different tau clusters")
options.parseArguments()

process = cms.Process("DiMuonTauETauHadTreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runDiMuTauETauHad_cfg.py isMC=1 ###############
##########################################################################

if options.isMC == 1:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuTauETauHadSelectorMC_cfi")

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuTauETauHadSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# --- please specify the sample that you need to run for local test ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

######### embed 2017v2 tauID into the miniAOD ###############
# reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Rerunning_of_the_tau_ID_on_M_AN1
if options.tauCluster <= 0:
    print " ====== use slimmedTaus cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTaus import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()

elif options.tauCluster == 1:
    print " ====== use slimmedTausBoosted (lower Tau Pt) cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTausBoosted import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()

elif options.tauCluster == 2:
    print " ====== use slimmedTausMuonCleaned cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTausMuonCleaned import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()

else:
    print " ====== use slimmedTausElectronCleaned cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTausElectronCleaned import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()
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
        process.rerunMvaIsolationSequence*
        process.NewTauIDsEmbedded*
        process.TauHadSelector*
        process.JetSelector*
        process.MuMuTauETauHadAnalyzer
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
        fileName =  cms.string('MuMuTauETauHadTreelization.root')
)

process.p = cms.Path(process.treelizer)
