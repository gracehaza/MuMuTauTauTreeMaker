# MuMuTauTauTreelizerThis tool is under developing...Introduction for setting up the environment:$ export SCRAM_ARCH=slc6_amd64_gcc630$ cmsrel CMSSW_9_4_9_cand2$ cd CMSSW_9_4_9_cand2/src/$ git clone https://github.com/Fengwangdong/MuMuTauTauTreeMaker.git$ scram b -j8Run the diMuon treelizer:$ cd MuMuTauTauTreeMaker/MuTauTreelizer/test$ cmsRun runDiMuon_cfg.pyThis script will load the EDM filters and analyzers defined in "MuTauTreelizer/python/DiMuSelector_cfi.py", in which one may need to modify the input parameter settings for the C++ modules defined in "MuTauTreelizer/plugins".Run the diMuon + ditau (tau_#mu + tau_h) treelizer:$ cd MuMuTauTauTreeMaker/MuTauTreelizer/test$ cmsRun runDiMuDiTau_cfg.pyThis script will load the EDM filters and analyzers defined in "MuTauTreelizer/python/DiMuDiTauSelector_cfi.py", in which one may need to modify the input parameter settings for the C++ modules defined in "MuTauTreelizer/plugins".Run the diMuon + ditau (tau_e + tau_h) treelizer:$ cd MuMuTauTauTreeMaker/MuTauTreelizer/test$ cmsRun runDiMuTauETauHad_cfg.pyThis script will load the EDM filters and analyzers defined in "MuTauTreelizer/python/DiMuTauETauHadSelector_cfi.py", in which one may need to modify the input parameter settings for the C++ modules defined in "MuTauTreelizer/plugins".