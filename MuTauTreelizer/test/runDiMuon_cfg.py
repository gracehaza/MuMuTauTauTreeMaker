import FWCore.ParameterSet.Config as cms

process = cms.Process("DiMuonTreelizer")

process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data or MC: ##############
isMC = True
##########################################################################

if isMC == True:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuSelectorMC_cfi")

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# --- please specify the sample that you need to run for local test ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_100.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_101.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_102.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_103.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_104.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_105.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_106.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_107.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_108.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_109.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_110.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_111.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_112.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_113.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_114.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_115.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_116.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_117.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_118.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_119.root',
        '/store/group/phys_higgs/HiggsExo/fengwang/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MiniAOD_DYJetsToLL_M50_Fall17DRPremix_v1/190806_145316/0000/mumutautau_zskim_120.root',
        #'/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_zskim_6131.root',
        #'/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_zskim_6130.root',
        #'/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_zskim_6132.root',
        #'/store/group/phys_higgs/HiggsExo/fengwang/SingleMuon/MiniAOD_SMu_DataF_17Nov2017_v1/190327_151854/0006/mumutautau_zskim_6133.root',
    )
)

############################################################

process.treelizer = cms.Sequence(
        process.lumiTree*
        process.HLTEle*
        process.TrigMuMatcher*
        process.MuonPtEtaCut*
        process.MuonID*
        process.LeadingMuonIso*
        process.SecondMuonIso*
        process.DiMuonMassSelector*
        process.ThirdMuonIso*
        process.JetSelector*
        process.DiMuonAnalyzer
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
        fileName =  cms.string('MuMuTreelization.root')
)

process.p = cms.Path(process.treelizer)
