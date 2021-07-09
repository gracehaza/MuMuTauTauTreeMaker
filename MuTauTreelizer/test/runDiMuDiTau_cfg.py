import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/MiniAOD_H125AA19_DiMuDiTau_Fall17DRPremix_v1/190515_140053/0000/mumutautau_1.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.register('tauCluster', 2, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "different tau clusters")
options.register('muTrigger', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "number of muons triggered")
options.register('numThreads', 8, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Set number of CPU cores")
options.parseArguments()

process = cms.Process("DiMuonDiTauTreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runDiMuDiTau_cfg.py isMC=1 ###############
##########################################################################
if options.isMC == 1:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuDiTauSelectorMC_cfi")

    if options.muTrigger == 2:
        process.HLTFilter = cms.Sequence(process.HLTEleDiMu)

        process.TrigRecoMuMatcher = cms.Sequence(process.TrigDiMuMatcher)
#        process.DiMuDiTauEDAnalyzer = cms.Sequence(process.DiMuDiTauAnalyzerDiMuTrig)

    else:
        process.HLTFilter = cms.Sequence(process.HLTEle)
        process.TrigRecoMuMatcher = cms.Sequence(process.TrigMuMatcher)
#        process.DiMuDiTauEDAnalyzer = cms.Sequence(process.DiMuDiTauAnalyzer)

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuDiTauSelector_cfi")

    if options.muTrigger == 2:
        process.HLTFilter = cms.Sequence(process.HLTEleDiMu)
        process.TrigRecoMuMatcher = cms.Sequence(process.TrigDiMuMatcher)
        process.DiMuDiTauEDAnalyzer = cms.Sequence(process.DiMuDiTauAnalyzerDiMuTrig)

    else:
        process.HLTFilter = cms.Sequence(process.HLTEle)
        process.TrigRecoMuMatcher = cms.Sequence(process.TrigMuMatcher)
        process.DiMuDiTauEDAnalyzer = cms.Sequence(process.DiMuDiTauAnalyzer)

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
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTaus as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 1:
    print " ====== use slimmedTausMuonCleaned cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausMuonCleaned as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 2:
    print " ====== use slimmedTausElectronCleaned cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausElectronCleaned as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 3:
    print " ====== use slimmedTausMuonCleanedMedium cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausMuonCleanedMedium as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 4:
    print " ====== use slimmedTausElectronCleanedMedium cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausElectronCleanedMedium as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 5:
    print " ====== use slimmedTausMuonCleanedTight cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausMuonCleanedTight as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 6:
    print " ====== use slimmedTausElectronCleanedTight cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausElectronCleanedTight as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 7:
    print " ====== use slimmedTausBoosted cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausBoosted as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = False,
            updatedTauName = updatedTauName,
            PATTauProducer = cms.InputTag('cleanedSlimmedTausBoosted'),
            srcChargedIsoPtSum = cms.string('chargedIsoPtSumNoOverLap'),
            srcNeutralIsoPtSum = cms.string('neutralIsoPtSumNoOverLap'),
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.ca8PFJetsCHSprunedForBoostedTausPAT * getattr(process, "cleanedSlimmedTausBoosted") * process.rerunMvaIsolationSequence * getattr(process,updatedTauName))
############################################################

if options.isMC == 1:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTFilter*
            process.MuonID*
            process.MuonSelector*
            process.TrigRecoMuMatcher*
            process.ElectronCandSelector*
            process.rerunTauIDSequence*
            process.TauCandSelector*
            process.DeepDiTauProducer*
            process.JetIdEmbedder*
            process.GenMuonCandSelector*
            process.GenElectronCandSelector*
            process.GenTauMuCandSelector*
            process.GenTauEleCandSelector*
            process.GenTauHadCandSelector*
            process.DiMuDiTauAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('MuMuTauTauTreelization_mc_07Jul2021.root')
    )

else:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTFilter*
            process.MuonID*
            process.MuonSelector*
            process.TrigRecoMuMatcher*
            process.ElectronCandSelector*
            process.rerunTauIDSequence*
            process.TauCandSelector*
            process.DeepDiTauProducer*
            process.JetIdEmbedder*
            process.DiMuDiTauAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('MuMuTauTauTreelization_data.root')
    )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.options.numberOfThreads = cms.untracked.uint32(options.numThreads)
process.p = cms.Path(process.treelizer)
