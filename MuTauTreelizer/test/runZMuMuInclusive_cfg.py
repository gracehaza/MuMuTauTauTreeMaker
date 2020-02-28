import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_100.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.register('tauCluster', 2, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "different tau clusters")
options.parseArguments()

process = cms.Process("ZMuMuInclusiveTreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runZMuMuInclusive_cfg.py isMC=1 ###############
##########################################################################

if options.isMC == 1:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.ZMuMuInclusiveSelectorMC_cfi")

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.ZMuMuInclusiveSelector_cfi")

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
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * process.slimmedTausNewID)

elif options.tauCluster == 1:
    print " ====== use slimmedTausBoosted (lower Tau Pt) cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTausBoosted import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * process.slimmedTausNewID)

elif options.tauCluster == 2:
    print " ====== use slimmedTausMuonCleaned cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTausMuonCleaned import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * process.slimmedTausNewID)

elif options.tauCluster == 3:
    print " ====== use slimmedTausElectronCleaned cluster ======"
    from MuMuTauTauTreeMaker.MuTauTreelizer.TauIdMVA_slimmedTausElectronCleaned import *
    myTool = TauIDEmbedder(process, cms,
            debug = True,
            toKeep = ["2017v2"]
    )
    myTool.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * process.slimmedTausNewID)

else:
    print " ====== use slimmedTausNewID cluster ======"
    updatedTauName = "slimmedTausNewID"
    import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = False,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))
############################################################

if options.isMC == 1:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonID*
            process.MuonSelector*
            process.TrigMuMatcher*
            process.ElectronCandSelector*
            process.rerunTauIDSequence*
            process.TauCandSelector*
            process.JetSelector*
            process.GenMuonCandSelector*
            process.GenElectronCandSelector*
            process.GenTauMuCandSelector*
            process.GenTauEleCandSelector*
            process.GenTauHadCandSelector*
            process.ZMuMuInclusiveAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('ZMuMuTreelization_mc.root')
    )

else:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonID*
            process.MuonSelector*
            process.TrigMuMatcher*
            process.ElectronCandSelector*
            process.rerunTauIDSequence*
            process.TauCandSelector*
            process.JetSelector*
            process.ZMuMuInclusiveAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('ZMuMuTreelization_data.root')
    )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.p = cms.Path(process.treelizer)
