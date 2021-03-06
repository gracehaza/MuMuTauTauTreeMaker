// -*- C++ -*-
//
// Package:    MuMuChannel/DiMuDiTauAnalyzer
// Class:      DiMuDiTauAnalyzer
// 
/**\class DiMuDiTauAnalyzer DiMuDiTauAnalyzer.cc MuMuChannel/DiMuDiTauAnalyzer/plugins/DiMuDiTauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 09 Apr 2019 16:30:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <math.h>
#include <string>
#include <iostream>

using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DiMuDiTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DiMuDiTauAnalyzer(const edm::ParameterSet&);
      ~DiMuDiTauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      std::vector<const reco::Candidate*> findTauMuEleVisDaughters(const reco::Candidate*);
      std::vector<const reco::Candidate*> findTauHadVisDaughters(const reco::Candidate*);
      int findTauPiZeros(const reco::Candidate*);
      int findTauChargedHadrons(const reco::Candidate*);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> MuTag;
      edm::EDGetTokenT<edm::View<pat::Electron>> EleTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauTag;
      edm::EDGetTokenT<edm::View<pat::Jet>> JetTag;
      edm::EDGetTokenT<edm::View<pat::MET>> MetTag;
      edm::EDGetTokenT<edm::View<reco::Vertex>> VertexTag;
      edm::EDGetTokenT<edm::View<pat::Jet>> slimmedJetTag;
      edm::EDGetTokenT<double> rhoTag;
      EffectiveAreas effectiveAreas;
      bool isMC;
      edm::EDGetTokenT<edm::View<PileupSummaryInfo>> PileupTag;
      edm::EDGetTokenT<GenEventInfoProduct> generator;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenMuTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenEleTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenTauMuTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenTauEleTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> GenTauHadTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> GenParticleTag;
   
/*
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);
  */
      TTree *objectTree;
      // --- below is the vectors of object variables ---
      
      // --- reconstructed muons ---
      vector<float> recoMuonPt;
      vector<float> recoMuonEta;
      vector<float> recoMuonPhi;
      vector<float> recoMuonEnergy;
      vector<int> recoMuonPDGId;
      vector<float> recoMuonIsolation;
      vector<float> recoMuonDXY;
      vector<float> recoMuonDZ;
      vector<int> recoMuonNTrackerLayers;
      vector<int> recoMuonTriggerFlag;
      vector<int> recoMuonRefToElectron;
      vector<int> recoMuonRefToTau;

      // --- reconstructed electrons ---
      vector<float> recoElectronPt;
      vector<float> recoElectronEta;
      vector<float> recoElectronPhi;
      vector<float> recoElectronEnergy;
      vector<int> recoElectronPDGId;
      vector<float> recoElectronIsolation;
      vector<int> recoElectronIdLoose;
      vector<int> recoElectronIdMedium;
      vector<int> recoElectronIdTight;
      vector<float> recoElectronEcalTrkEnergyPostCorr;
      vector<float> recoElectronEcalTrkEnergyErrPostCorr;
      vector<int> recoElectronRefToMuon;
      vector<int> recoElectronRefToTau;

      // --- reconstructed taus ---
      vector<float> recoTauPt;
      vector<float> recoTauEta;
      vector<float> recoTauPhi;
      vector<float> recoTauEnergy;
      vector<int> recoTauPDGId;
      vector<float> recoTauDecayMode;
      vector<float> recoTauIsoMVArawValue;
      vector<float> recoTauIsoMVAVVLoose;
      vector<float> recoTauIsoMVAVLoose;
      vector<float> recoTauIsoMVALoose;
      vector<float> recoTauIsoMVAMedium;
      vector<float> recoTauIsoMVATight;
      vector<float> recoTauIsoMVAVTight;
      vector<float> recoTauIsoMVAVVTight;
      
      vector<float> recoTauAntiMuMVALoose;
      vector<float> recoTauAntiMuMVATight;

      vector<float> recoTauAntiEleMVArawValue;
      vector<float> recoTauAntiEleMVAVLoose;
      vector<float> recoTauAntiEleMVALoose;
      vector<float> recoTauAntiEleMVAMedium;
      vector<float> recoTauAntiEleMVATight;
      vector<float> recoTauAntiEleMVAVTight;

      // *** deep tau ID variables ***
      vector<float> recoTauDeepVSeraw;
      vector<float> recoTauDeepVSjetraw;
      vector<float> recoTauDeepVSmuraw;
 
      vector<float> recoTauDeepVSeLoose;
      vector<float> recoTauDeepVSjetLoose;
      vector<float> recoTauDeepVSmuLoose;

      vector<float> recoTauDeepVSeMedium;
      vector<float> recoTauDeepVSjetMedium;
      vector<float> recoTauDeepVSmuMedium;

      vector<float> recoTauDeepVSeTight;
      vector<float> recoTauDeepVSjetTight;
      vector<float> recoTauDeepVSmuTight;

      vector<float> recoTauDeepVSeVLoose;
      vector<float> recoTauDeepVSjetVLoose;
      vector<float> recoTauDeepVSmuVLoose;

      vector<float> recoTauDeepVSeVTight;
      vector<float> recoTauDeepVSjetVTight;

      vector<float> recoTauDeepVSeVVLoose;
      vector<float> recoTauDeepVSjetVVLoose;

      vector<float> recoTauDeepVSeVVTight;
     vector<float> recoTauDeepVSjetVVTight;

      vector<float> recoTauDeepVSeVVVLoose;
      vector<float> recoTauDeepVSjetVVVLoose;

      vector<int> recoTauRefToMuon;
      vector<int> recoTauRefToElectron;
  
      // --- reconstructed jets ---
      vector<float> recoJetPt;
      vector<float> recoJetEta;
      vector<float> recoJetPhi;
      vector<float> recoJetEnergy;
      vector<float> recoJetCSV;
        
      // --- reconstructed MET ---
      vector<float> recoMET;
      vector<float> recoMETPhi;
      vector<float> recoMETPx;
      vector<float> recoMETPy;

      // --- pileup and reconstructed vertices ---
      int recoNPrimaryVertex;
      int recoNPU;
      int trueNInteraction;
  int eventID;

  vector<float> oldDeepDiTauValue;
  vector<float> oldDeepDiTauValueMD;

  vector<float> DeepDiTauValue_genmatched;
  vector<float> DeepDiTauValueMD_genmatched;

  vector<float> DeepDiTau_boosted_nolepton_charmValue;
  vector<float> DeepDiTau_boosted_nolepton_charm_massdecoValue;
  vector<float> DeepDiTau_boosted_nolepton_massdecoValue;
  vector<float> DeepDiTau_boosted_noleptonValue;

  vector<float> DeepDiTau_nolepton_charm_massdecoValue;
  vector<float> DeepDiTau_nolepton_charmValue;
  vector<float> DeepDiTau_nolepton_massdecoValue;
  vector<float> DeepDiTau_noleptonValue;


  vector<float> jet_index;
  vector<float> jet_daughters_index;
  vector<float> jet_daughters_pdgid;
  int jetindex = 0;
  // int jetdaughterindex = 0;

  vector<float> jet_pt;
  vector<float> jet_eta;
  vector<float> jet_phi;
  vector<float> jet_energy;
  vector<float> jet_mass;
  vector<float> jet_jetCharge;

  vector<float> jet_chargedMultiplicity;
  vector<float> jet_neutralMultiplicity;

  vector<float> jet_chargedHadronMultiplicity;
  vector<float> jet_neutralHadronMultiplicity;
  vector<float> jet_electronMultiplicity;
  vector<float> jet_muonMultiplicity;
  vector<float> jet_photonMultiplicity;

  vector<float> jet_chargedEmEnergy;
  vector<float> jet_neutralEmEnergy;
  vector<float> jet_chargedHadronEnergy;
  vector<float> jet_neutralHadronEnergy;
  vector<float> jet_muonEnergy;
  vector<float> jet_electronEnergy;
  vector<float> jet_photonEnergy;

  vector<float> jet_chargedEmEnergyFraction;
  vector<float> jet_neutralEmEnergyFraction;
  vector<float> jet_chargedHadronEnergyFraction;
  vector<float> jet_neutralHadronEnergyFraction;
  vector<float> jet_muonEnergyFraction;
  vector<float> jet_electronEnergyFraction;
  vector<float> jet_photonEnergyFraction;

  vector<float> jet_n60;
  vector<float> jet_n90;

  vector<float> jet_hadronFlavour;
  vector<float> jet_partonFlavour;

  vector<float> jet_pfJetBProbabilityBJetTags;
  vector<float> jet_pfJetProbabilityBJetTags;
  vector<float> jet_pfTrackCountingHighEffBJetTags;
  vector<float> jet_pfSimpleSecondaryVertexHighEffBJetTags;
  vector<float> jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags;
  vector<float> jet_pfCombinedSecondaryVertexV2BJetTags;
  vector<float> jet_pfCombinedInclusiveSecondaryVertexV2BJetTags;
  vector<float> jet_softPFMuonBJetTags;
  vector<float> jet_softPFElectronBJetTags;
  vector<float> jet_pfCombinedMVAV2BJetTags;
  vector<float> jet_pfCombinedCvsLJetTags;
  vector<float> jet_pfCombinedCvsBJetTags;
  vector<float> jet_pfDeepCSVJetTags_probb;
  vector<float> jet_pfDeepCSVJetTags_probc;
  vector<float> jet_pfDeepCSVJetTags_probudsg;
  vector<float> jet_pfDeepCSVJetTags_probbb;

  vector<float> jet_numberOfDaughters;

  vector<float> jet_daughter_pt;
  vector<float> jet_daughter_eta;
  vector<float> jet_daughter_phi;
  vector<float> jet_daughter_energy;
  vector<float> jet_daughter_mass;
  vector<float> jet_daughter_charge;

  vector<float> jet_daughter_etaAtVtx;
  vector<float> jet_daughter_phiAtVtx;
  vector<float> jet_daughter_vx;
  vector<float> jet_daughter_vy;
  vector<float> jet_daughter_vz;
  vector<float> jet_daughter_dxy;
  vector<float> jet_daughter_dxyError;
  vector<float> jet_daughter_dz;
  vector<float> jet_daughter_dzError;

  vector<float> jet_daughter_pdgId;

  vector<float> jet_daughter_pixelLayersWithMeasurement;
  vector<float> jet_daughter_stripLayersWithMeasurement;
  vector<float> jet_daughter_trackerLayersWithMeasurement;
  vector<float> jet_daughter_trackHighPurity;

  vector<float> jet_daughter_puppiWeight;
  vector<float> jet_daughter_puppiWeightNoLep;

  vector<float> jet_daughter_isIsolatedChargedHadron;
  vector<float> jet_daughter_isStandAloneMuon;
  vector<float> jet_daughter_isTrackerMuon;
  vector<float> jet_daughter_isGlobalMuon;
  vector<float> jet_daughter_isGoodEgamma;


  // gen part truth
  vector<bool> jet_isTauHTauH;
  vector<bool> jet_isTauHTauM;
  vector<bool> jet_isTauHTauE;
  vector<bool> jet_isTauMTauM;
  vector<bool> jet_isTauMTauE;
  vector<bool> jet_isTauETauE;
  vector<bool> jet_isTauH;
  vector<bool> jet_isTauM;
  vector<bool> jet_isTauE;
  vector<bool> jet_isB;
  vector<bool> jet_isC;
  vector<bool> jet_isS;
  vector<bool> jet_isUD;
  vector<bool> jet_isG;

  /*  vector<float> DeepDiTaujetPt;
  vector<float> DeepDiTaujetEta;
  vector<float> DeepDiTaujetPhi;
  vector<float> DeepDiTaujetEnergy;
  */ 


      // --- gen muons ----
      vector<float> genMuonPt;
      vector<float> genMuonEta;
      vector<float> genMuonPhi;
      vector<float> genMuonMass;
      vector<int> genMuonPDGId;
      vector<int> genMuonMotherPDGId;

      // --- gen electrons ----
      vector<float> genElectronPt;
      vector<float> genElectronEta;
      vector<float> genElectronPhi;
      vector<float> genElectronMass;
      vector<int> genElectronPDGId;
      vector<int> genElectronMotherPDGId;

      // --- gen tau_mu ----
      vector<float> genTauMuPt;
      vector<float> genTauMuEta;
      vector<float> genTauMuPhi;
      vector<float> genTauMuMass;
      vector<int> genTauMuPDGId;
      vector<int> genTauMuMotherPDGId;
      vector<float> genTauMuVisPt;
      vector<float> genTauMuVisMass;

      // --- gen tau_e ----
      vector<float> genTauElePt;
      vector<float> genTauEleEta;
      vector<float> genTauElePhi;
      vector<float> genTauEleMass;
      vector<int> genTauElePDGId;
      vector<int> genTauEleMotherPDGId;
      vector<float> genTauEleVisPt;
      vector<float> genTauEleVisMass;

      // --- gen tau_h ----
      vector<float> genTauHadPt;
      vector<float> genTauHadEta;
      vector<float> genTauHadPhi;
      vector<float> genTauHadMass;
      vector<int> genTauHadPDGId;
      vector<int> genTauHadMotherPDGId;
      vector<float> genTauHadVisPt;
      vector<float> genTauHadVisMass;
      vector<int> genTauHadNPionZero;
      vector<int> genTauHadNChargedHadrons;

      // --- event weight for MC ---
      float genEventWeight; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiMuDiTauAnalyzer::DiMuDiTauAnalyzer(const edm::ParameterSet& iConfig):
    MuTag(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuTag"))),
    EleTag(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("EleTag"))),
    TauTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauTag"))),
    JetTag(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
    MetTag(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("MetTag"))),
    VertexTag(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexTag"))),
    slimmedJetTag(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("slimmedJetTag"))),
    rhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"))),
    effectiveAreas((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
    PileupTag(consumes<edm::View<PileupSummaryInfo>>(iConfig.existsAs<edm::InputTag>("PileupTag") ? iConfig.getParameter<edm::InputTag>("PileupTag") : edm::InputTag())),
    generator(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag())),
    GenMuTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenMuTag") ? iConfig.getParameter<edm::InputTag>("GenMuTag") : edm::InputTag())),
    GenEleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenEleTag") ? iConfig.getParameter<edm::InputTag>("GenEleTag") : edm::InputTag())),
  GenTauMuTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauMuTag") ? iConfig.getParameter<edm::InputTag>("GenTauMuTag") : edm::InputTag())),
  GenTauEleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauEleTag") ? iConfig.getParameter<edm::InputTag>("GenTauEleTag") : edm::InputTag())),
  GenTauHadTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauHadTag") ? iConfig.getParameter<edm::InputTag>("GenTauHadTag") : edm::InputTag())),
  GenParticleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenParticleTag") ? iConfig.getParameter<edm::InputTag>("GenParticleTag") : edm::InputTag()))
//GenParticleTag(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleTag")))  

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC = iConfig.getParameter<bool>("isMC");
}


DiMuDiTauAnalyzer::~DiMuDiTauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiMuDiTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu;
   iEvent.getByToken(MuTag, pMu);

   edm::Handle<edm::View<pat::Electron>> pElectron;
   iEvent.getByToken(EleTag, pElectron);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauTag, pTau);
   /*
   edm::Handle<edm::View<pat::Jet>> pJet;
   iEvent.getByToken(JetTag, pJet);
   */
   edm::Handle<edm::View<pat::MET>> pMet;
   iEvent.getByToken(MetTag, pMet);
  
   edm::Handle<edm::View<reco::Vertex>> pVertex;
   iEvent.getByToken(VertexTag, pVertex);

   edm::Handle<edm::View<pat::Jet>> pslimJet;
   iEvent.getByToken(slimmedJetTag, pslimJet);

   edm::Handle<double> pRho;
   iEvent.getByToken(rhoTag, pRho);


   edm::Handle<edm::View<reco::GenParticle>> pGenParticles;
   iEvent.getByToken(GenParticleTag, pGenParticles);


   if (isMC)
   {
       edm::Handle<edm::View<reco::GenParticle>> pGenMu;
       iEvent.getByToken(GenMuTag, pGenMu);

       edm::Handle<edm::View<reco::GenParticle>> pGenEle;
       iEvent.getByToken(GenEleTag, pGenEle);

       edm::Handle<edm::View<reco::GenParticle>> pGenTauMu;
       iEvent.getByToken(GenTauMuTag, pGenTauMu);

       edm::Handle<edm::View<reco::GenParticle>> pGenTauEle;
       iEvent.getByToken(GenTauEleTag, pGenTauEle);

       edm::Handle<edm::View<reco::GenParticle>> pGenTauHad;
       iEvent.getByToken(GenTauHadTag, pGenTauHad);

       edm::Handle<GenEventInfoProduct> gen_ev_info;
       iEvent.getByToken(generator, gen_ev_info);

       edm::Handle<edm::View<PileupSummaryInfo>> pileup_info;
       iEvent.getByToken(PileupTag, pileup_info);

       //       edm::Handle<edm::View<reco::GenParticle>> pGenParticles;
       //iEvent.getByToken(GenParticleTag, pGenParticles);

       if(pGenParticles->size() > 0){
	 for(edm::View<reco::GenParticle>::const_iterator iParticle=pGenParticles->begin(); iParticle!=pGenParticles->end(); ++iParticle){
	 }

       }

      if (pGenMu->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iMuon=pGenMu->begin(); iMuon!=pGenMu->end(); iMuon++)
           {
               genMuonPt.push_back(iMuon->pt());
               genMuonEta.push_back(iMuon->eta());
               genMuonPhi.push_back(iMuon->phi());
               genMuonMass.push_back(iMuon->mass());
               genMuonPDGId.push_back(iMuon->pdgId());
               genMuonMotherPDGId.push_back(iMuon->mother()->pdgId());
           } // end for loop on gen muons
       } // end if pGenMu->size() > 0

       if (pGenEle->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iElectron=pGenEle->begin(); iElectron!=pGenEle->end(); iElectron++)
           {
               genElectronPt.push_back(iElectron->pt());
               genElectronEta.push_back(iElectron->eta());
               genElectronPhi.push_back(iElectron->phi());
               genElectronMass.push_back(iElectron->mass());
               genElectronPDGId.push_back(iElectron->pdgId());
               genElectronMotherPDGId.push_back(iElectron->mother()->pdgId());
           } // end for loop on gen electrons
       } // end if pGenEle->size() > 0

       if (pGenTauMu->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iTau=pGenTauMu->begin(); iTau!=pGenTauMu->end(); iTau++)
           {
               genTauMuPt.push_back(iTau->pt());
               genTauMuEta.push_back(iTau->eta());
               genTauMuPhi.push_back(iTau->phi());
               genTauMuMass.push_back(iTau->mass());
               genTauMuPDGId.push_back(iTau->pdgId());
               genTauMuMotherPDGId.push_back(iTau->mother()->pdgId());

               TLorentzVector sumTauMuVis;
               std::vector <const reco::Candidate*> daughters;
               daughters.clear();

               for (unsigned int iDau = 0; iDau < iTau->numberOfDaughters(); iDau++)
               {
                   const reco::Candidate* directDaughter = iTau->daughter(iDau);
                   daughters = findTauMuEleVisDaughters(directDaughter); // collect all the current daughter (if status == 1) or together its daughters (if status != 1) 

                   for (unsigned int jDau = 0; jDau < daughters.size(); jDau++)
                   {
                       TLorentzVector p4Daughter;
                       p4Daughter.SetPtEtaPhiM(daughters[jDau]->pt(), daughters[jDau]->eta(), daughters[jDau]->phi(), daughters[jDau]->mass());
                       sumTauMuVis = sumTauMuVis + p4Daughter;
                   } // end for loop on all generations of visible daughter particles of tau_mu
               } // end for loop on tau_mu direct daughter particles

               genTauMuVisPt.push_back(sumTauMuVis.Pt());
               genTauMuVisMass.push_back(sumTauMuVis.M());
           } // end for loop on gen tau_mu
       } // end if pGenTauMu->size() > 0

       if (pGenTauEle->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iTau=pGenTauEle->begin(); iTau!=pGenTauEle->end(); iTau++)
           {
               genTauElePt.push_back(iTau->pt());
               genTauEleEta.push_back(iTau->eta());
               genTauElePhi.push_back(iTau->phi());
               genTauEleMass.push_back(iTau->mass());
               genTauElePDGId.push_back(iTau->pdgId());
               genTauEleMotherPDGId.push_back(iTau->mother()->pdgId());

               TLorentzVector sumTauEleVis;
               std::vector <const reco::Candidate*> daughters;
               daughters.clear();

               for (unsigned int iDau = 0; iDau < iTau->numberOfDaughters(); iDau++)
               {
                   const reco::Candidate* directDaughter = iTau->daughter(iDau);
                   daughters = findTauMuEleVisDaughters(directDaughter); // collect all the current daughter (if status == 1) or together its daughters (if status != 1) 

                   for (unsigned int jDau = 0; jDau < daughters.size(); jDau++)
                   {
                       TLorentzVector p4Daughter;
                       p4Daughter.SetPtEtaPhiM(daughters[jDau]->pt(), daughters[jDau]->eta(), daughters[jDau]->phi(), daughters[jDau]->mass());
                       sumTauEleVis = sumTauEleVis + p4Daughter;
                   } // end for loop on all generations of visible daughter particles of tau_e
               } // end for loop on tau_e direct daughter particles

               genTauEleVisPt.push_back(sumTauEleVis.Pt());
               genTauEleVisMass.push_back(sumTauEleVis.M());
           } // end for loop on gen tau_e
       } // end if pGenTauEle->size() > 0

       if (pGenTauHad->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iTau=pGenTauHad->begin(); iTau!=pGenTauHad->end(); iTau++)
           {
               genTauHadPt.push_back(iTau->pt());
               genTauHadEta.push_back(iTau->eta());
               genTauHadPhi.push_back(iTau->phi());
               genTauHadMass.push_back(iTau->mass());
               genTauHadPDGId.push_back(iTau->pdgId());
               genTauHadMotherPDGId.push_back(iTau->mother()->pdgId());

               TLorentzVector sumTauHadVis;
	       std::vector <const reco::Candidate*> daughters;
               daughters.clear();

               int nPiZeros = 0;
               int nChargedHadrons = 0;

               for (unsigned int iDau = 0; iDau < iTau->numberOfDaughters(); iDau++)
               {
                   const reco::Candidate* directDaughter = iTau->daughter(iDau);
                   daughters = findTauHadVisDaughters(directDaughter); // collect all the current daughter (if status == 1) or together its daughters (if status != 1) 

                   nPiZeros += findTauPiZeros(directDaughter);
                   nChargedHadrons += findTauChargedHadrons(directDaughter);

                   for (unsigned int jDau = 0; jDau < daughters.size(); jDau++)
                   {
                       TLorentzVector p4Daughter;
                       p4Daughter.SetPtEtaPhiM(daughters[jDau]->pt(), daughters[jDau]->eta(), daughters[jDau]->phi(), daughters[jDau]->mass());
                       sumTauHadVis = sumTauHadVis + p4Daughter;
		       //		       multipleTauHadVis = multipleTauHadVis +sumTauHadVis;
		   } // end for loop on all generations of visible daughter particles of tau_h
               } // end for loop on tau_h direct daughter particles
	       // multipleTauHadVis = multipleTauHadVis + sumTauHadVis;
               genTauHadVisPt.push_back(sumTauHadVis.Pt());
               genTauHadVisMass.push_back(sumTauHadVis.M());
               genTauHadNPionZero.push_back(nPiZeros);
               genTauHadNChargedHadrons.push_back(nChargedHadrons);
           } // end for loop on gen tau_h
       } // end if pGenTauHad->size() > 0

       if (gen_ev_info.isValid())
       {
           genEventWeight = gen_ev_info->weight();
       } // end if gen_ev_info.isValid()

       if (pileup_info.isValid())
       {
           for(edm::View<PileupSummaryInfo>::const_iterator iPileup=pileup_info->begin(); iPileup!=pileup_info->end(); iPileup++)
           {
               if (iPileup->getBunchCrossing() == 0)
               {
                   trueNInteraction = iPileup->getTrueNumInteractions();
                   recoNPU = iPileup->getPU_NumInteractions();
               } // end if iPileup->getBunchCrossing() == 0
           } // end for loop on pileup_info
       } // end if pileup_info.isValid()
   } // end if isMC == true

   // --- prepare for offline primary vertices ---
   recoNPrimaryVertex = 0; 
   if (pVertex.isValid())
   {
       for(edm::View<reco::Vertex>::const_iterator iPV=pVertex->begin(); iPV!=pVertex->end(); iPV++)
       {
           recoNPrimaryVertex++;
       } // end for loop on pVertex
   } // end if pVertex.isValid()

   eventID = iEvent.eventAuxiliary().event();

   // --- prepare muon vector ---
   if (pMu->size()>0)
   {
       int muonCounter = 0;
       for(edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
       {
           recoMuonPt.push_back(iMuon->pt());
           recoMuonEta.push_back(iMuon->eta());
           recoMuonPhi.push_back(iMuon->phi());
           recoMuonEnergy.push_back(iMuon->energy());
           recoMuonPDGId.push_back(iMuon->pdgId());
           reco::MuonPFIsolation iso = iMuon->pfIsolationR04();
           double reliso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / iMuon->pt();
           recoMuonIsolation.push_back(reliso);
           recoMuonDXY.push_back(iMuon->muonBestTrack()->dxy());
           recoMuonDZ.push_back(iMuon->muonBestTrack()->dz());
           recoMuonNTrackerLayers.push_back(iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
           recoMuonRefToElectron.push_back(0);
           recoMuonRefToTau.push_back(0);

           if (muonCounter == 0)
           {
               recoMuonTriggerFlag.push_back(1);
               muonCounter++;
           } // end if muonCounter == 0

           else{
               recoMuonTriggerFlag.push_back(0);
               muonCounter++;
           } // end else if muonCounter != 0
       } // end for loop on muons
   } // end if pMu->size()>0

   // --- prepare electron vector ---
   if (pElectron->size()>0)
   {
       for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
       {
           int isLoose = 0;
           int isMedium = 0;
           int isTight = 0;

           // ---  full5x5_sigmaIetaIeta ---
           float sigmaIetaIeta = iElectron->full5x5_sigmaIetaIeta();

           // --- fabs(dEtaSeed) ---
           float dEtaSeed = fabs(iElectron->superCluster().isNonnull() && iElectron->superCluster()->seed().isNonnull() ? iElectron->deltaEtaSuperClusterTrackAtVtx() - iElectron->superCluster()->eta() + iElectron->superCluster()->seed()->eta() : std::numeric_limits<float>::max()); 
           
           // --- fabs(dPhiIn) ---
           float dPhiIn = fabs(iElectron->deltaPhiSuperClusterTrackAtVtx());
           
           // --- variables for H/E cuts ---
           float HoE = iElectron->hadronicOverEm();
           float rho = pRho.isValid() ? (*pRho) : 0; 
           float energy = iElectron->superCluster()->energy();

           // --- variables for relIsoWithEffectiveArea ---
           float chad = iElectron->pfIsolationVariables().sumChargedHadronPt;
           float nhad = iElectron->pfIsolationVariables().sumNeutralHadronEt;
           float pho = iElectron->pfIsolationVariables().sumPhotonEt;
           float elePt = iElectron->pt();
           float eleEta = iElectron->superCluster()->eta();
           float eArea = effectiveAreas.getEffectiveArea(fabs(eleEta));
           float relIsoWithEffectiveArea = (chad + std::max(0.0f, nhad + pho - rho*eArea)) / elePt;

           // --- variables for fabs(1/E-1/p) ---
           float eInverseMinusPInverse = fabs(1.0 - iElectron->eSuperClusterOverP())*(1.0/iElectron->ecalEnergy());

           // --- expected missing inner hits ---
           int mHits = iElectron->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);

           // --- pass conversion veto ---
           bool isPassConVeto = iElectron->passConversionVeto();

           // ========= select electrons in different cut-based ID accordingly ==========
           if (fabs(eleEta) <= 1.479)
           {
               isLoose = (sigmaIetaIeta < 0.0112) &&
                         (dEtaSeed < 0.00377) &&
                         (dPhiIn < 0.0884) &&
                         (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.112 + 0.506/elePt) &&
                         (eInverseMinusPInverse < 0.193) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isMedium = (sigmaIetaIeta < 0.0106) &&
                          (dEtaSeed < 0.0032) &&
                          (dPhiIn < 0.0547) &&
                          (HoE < 0.046 + 1.16/energy + 0.0324*rho/energy) &&
                          (relIsoWithEffectiveArea < 0.0478 + 0.506/elePt) &&
                          (eInverseMinusPInverse < 0.184) &&
                          (mHits <= 1) &&
                          (isPassConVeto == true);

               isTight = (sigmaIetaIeta < 0.0104) &&
                         (dEtaSeed < 0.00255) &&
                         (dPhiIn < 0.022) &&
                         (HoE < 0.026 + 1.15/energy + 0.0324*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.0287 + 0.506/elePt) &&
                         (eInverseMinusPInverse < 0.159) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);
           }// endif (fabs(eleEta) <= 1.479)

           else{
               isLoose = (sigmaIetaIeta < 0.0425) &&
                         (dEtaSeed < 0.00674) &&
                         (dPhiIn < 0.169) &&
                         (HoE < 0.0441 + 2.54/energy + 0.183*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.108 + 0.963/elePt) &&
                         (eInverseMinusPInverse < 0.111) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isMedium = (sigmaIetaIeta < 0.0387) &&
                          (dEtaSeed < 0.00632) &&
                          (dPhiIn < 0.0394) &&
                          (HoE < 0.0275 + 2.52/energy + 0.183*rho/energy) &&
                          (relIsoWithEffectiveArea < 0.0658 + 0.963/elePt) &&
                          (eInverseMinusPInverse < 0.0721) &&
                          (mHits <= 1) &&
                          (isPassConVeto == true);

               isTight = (sigmaIetaIeta < 0.0353) &&
                         (dEtaSeed < 0.00501) &&
                         (dPhiIn < 0.0236) &&
                         (HoE < 0.0188 + 2.06/energy + 0.183*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.0445 + 0.963/elePt) &&
                         (eInverseMinusPInverse < 0.0197) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);
           } // end else (fabs(eleEta) > 1.479)

           recoElectronPt.push_back(iElectron->pt());
           recoElectronEta.push_back(iElectron->eta());
           recoElectronPhi.push_back(iElectron->phi());
           recoElectronEnergy.push_back(iElectron->energy());
           recoElectronPDGId.push_back(iElectron->pdgId());
           recoElectronIsolation.push_back(relIsoWithEffectiveArea);
           recoElectronIdLoose.push_back(isLoose);
           recoElectronIdMedium.push_back(isMedium);
           recoElectronIdTight.push_back(isTight);
           recoElectronEcalTrkEnergyPostCorr.push_back(iElectron->userFloat("ecalTrkEnergyPostCorr"));
           recoElectronEcalTrkEnergyErrPostCorr.push_back(iElectron->userFloat("ecalTrkEnergyErrPostCorr"));

           recoElectronRefToMuon.push_back(0);
           recoElectronRefToTau.push_back(0);
       } // end for loop on electrons
   } // end if pElectron->size()>0
   
   // --- prepare tau vector ---
   if (pTau->size()>0)
   {
       for(edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
       {
           recoTauPt.push_back(iTau->pt());
           recoTauEta.push_back(iTau->eta());
           recoTauPhi.push_back(iTau->phi());
           recoTauEnergy.push_back(iTau->energy());
           recoTauPDGId.push_back(iTau->pdgId());
           recoTauDecayMode.push_back(iTau->decayMode());
           recoTauRefToMuon.push_back(0);
           recoTauRefToElectron.push_back(0);

           if (iTau->isTauIDAvailable("byIsolationMVArun2017v2DBoldDMwLTraw2017"))
           {
               recoTauIsoMVArawValue.push_back(iTau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
               recoTauIsoMVAVVLoose.push_back(iTau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAVLoose.push_back(iTau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVALoose.push_back(iTau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAMedium.push_back(iTau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVATight.push_back(iTau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAVTight.push_back(iTau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAVVTight.push_back(iTau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
           } // end if TauMVA ID available

           if (iTau->isTauIDAvailable("byDeepTau2017v2p1VSjetraw"))
           {
               recoTauDeepVSeraw.push_back(iTau->tauID("byDeepTau2017v2p1VSeraw"));
	       recoTauDeepVSjetraw.push_back(iTau->tauID("byDeepTau2017v2p1VSjetraw"));
               recoTauDeepVSmuraw.push_back(iTau->tauID("byDeepTau2017v2p1VSmuraw"));

               recoTauDeepVSeLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSe"));
               recoTauDeepVSjetLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSmu"));

               recoTauDeepVSeMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSe"));
	       recoTauDeepVSjetMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSmu"));

               recoTauDeepVSeTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSe"));
	       recoTauDeepVSjetTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSmu"));

               recoTauDeepVSeVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSe"));
	       recoTauDeepVSjetVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSmu"));

               recoTauDeepVSeVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSjet"));

               recoTauDeepVSeVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSe"));
	       recoTauDeepVSjetVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSjet"));

               recoTauDeepVSeVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSe"));
	       recoTauDeepVSjetVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSjet"));

               recoTauDeepVSeVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSjet"));
           } // end if DeepTau ID available

           recoTauAntiMuMVALoose.push_back(iTau->tauID("againstMuonLoose3"));
           recoTauAntiMuMVATight.push_back(iTau->tauID("againstMuonTight3"));
       
           recoTauAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw"));
           recoTauAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA6"));
           recoTauAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA6"));
           recoTauAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA6"));
           recoTauAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA6"));
           recoTauAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA6"));
       } // end for loop on taus
   } // end if pTau->size()>0

   // --- prepare for the common source particle reference records for muon/electron candidates
   if (pElectron->size()>0 && pMu->size()>0)
   {
       int refMuonValue = 1;
       int countMuon = 0;
       for (edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
       {
           bool findMatchedMu = false;

           int refElectronValue = refMuonValue;
           int countElectron = 0;
           for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
           {
               bool findMatchedEle = false;
               for (unsigned int indexMuon = 0; indexMuon < iMuon->numberOfSourceCandidatePtrs(); indexMuon++)
               {
                   for (unsigned int indexEle = 0; indexEle < iElectron->numberOfSourceCandidatePtrs(); indexEle++)
                   {
                       if (iElectron->sourceCandidatePtr(indexEle).key() == iMuon->sourceCandidatePtr(indexMuon).key())
                       {
                           findMatchedMu = true;
                           findMatchedEle = true;
                           break;
                       } // end if find the same source of electron and muon
                   } // end for loop on electron source particles
                   if (findMatchedEle) break;
               } // end for loop on muon source particles
               
               if (findMatchedEle)
               {
                   recoElectronRefToMuon.at(countElectron) = refElectronValue;
               } // end if findMatchedEle

               countElectron += 1;
           } // end for loop on electron candidates

           if (findMatchedMu)
           {
               recoMuonRefToElectron.at(countMuon) = refMuonValue;
               refMuonValue += 1;
           } // end if findMatchedMu

           countMuon += 1;
       } // end for loop on muon candidates
   } // end if pElectron->size()>0 && pMu->size()>0

   // --- prepare muon VS tau reference for overlapped mu/tau candidates
   if (pMu->size()>0 && pTau->size()>0)
   {
       int refTauValue = 1;
       int countTau = 0;
       for (edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
       {
           bool findMatchedTau = false;

           int refMuonValue = refTauValue;
           int countMuon = 0;
           for (edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
           {
               bool findMatchedMu = false;
               for (unsigned int indexMu = 0; indexMu < iMuon->numberOfSourceCandidatePtrs(); indexMu++)
               {
                   for (unsigned int indexTau = 0; indexTau < iTau->numberOfSourceCandidatePtrs(); indexTau++)
                   {
                       if (iMuon->sourceCandidatePtr(indexMu).key() == iTau->sourceCandidatePtr(indexTau).key())
                       {
                           findMatchedMu = true;
                           findMatchedTau = true;
                           break;
                       } // end if the muon source and tau source have the same key
                   } // end for loop on tau source particles
                   if (findMatchedMu) break;
               } // end for loop on mu source particles

               if (findMatchedMu)
               {
                   recoMuonRefToTau.at(countMuon) = refMuonValue;
               } // end if findMatchedMu

               countMuon += 1;
           } // end for loop on muon candidates

           if (findMatchedTau)
           {
               recoTauRefToMuon.at(countTau) = refTauValue;
               refTauValue += 1;
           } // end if findMatchedTau

           countTau += 1;
       } // end for loop on tau candidates
   } // end if pMu->size()>0 && pTau->size()>0

   // --- prepare electron VS tau reference for overlapped ele/tau candidates
   if (pElectron->size()>0 && pTau->size()>0)
   {
       int refTauValue = 1;
       int countTau = 0;
       for (edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
       {
           bool findMatchedTau = false;

           int refElectronValue = refTauValue;
           int countElectron = 0;
           for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
           {
               bool findMatchedEle = false;
               for (unsigned int indexTau = 0; indexTau < iTau->numberOfSourceCandidatePtrs(); indexTau++)
               {
                   for (unsigned int indexEle = 0; indexEle < iElectron->numberOfSourceCandidatePtrs(); indexEle++)
                   {
                       if (iElectron->sourceCandidatePtr(indexEle).key() == iTau->sourceCandidatePtr(indexTau).key())
                       {
                           findMatchedEle = true;
                           findMatchedTau = true;
                           break;
                       } // end if find the same source of electron and tau
                   } // end for loop on electron source particles
                   if (findMatchedEle) break;
               } // end for loop on tau source particles

               if (findMatchedEle)
               {
                   recoElectronRefToTau.at(countElectron) = refElectronValue;
               } // end if findMatchedEle

               countElectron += 1;
           } // end for loop on electron candidates

           if (findMatchedTau)
           {
               recoTauRefToElectron.at(countTau) = refTauValue;
               refTauValue += 1;
           } // end if findMatchedTau

           countTau += 1;
       } // end for loop on tau candidates
   } // end if pElectron->size()>0 && pTau->size()>0

   // --- prepare jet vector ---
   /*
   if (pJet->size()>0)
   {
       for(edm::View<pat::Jet>::const_iterator iJet=pJet->begin(); iJet!=pJet->end(); iJet++)
       {
           recoJetPt.push_back(iJet->pt());
           recoJetEta.push_back(iJet->eta());
           recoJetPhi.push_back(iJet->phi());
           recoJetEnergy.push_back(iJet->energy());
           // --- btag for jet ---
           // reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Jets
           recoJetCSV.push_back(iJet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
       } // end for loop on jets
   } // end if pJet->size()>0
   */
  
   if (pslimJet->size()>0)
     {
       for (edm::View<pat::Jet>::const_iterator islimJet=pslimJet->begin(); islimJet!=pslimJet->end(); islimJet++)
	 {
	   if(fabs(islimJet->eta()) <= 2.5){
	     if (islimJet->pt() > 20){

	       oldDeepDiTauValue.push_back(islimJet->userFloat("ditau2017v1"));
	       oldDeepDiTauValueMD.push_back(islimJet->userFloat("ditau2017MDv1"));
	       DeepDiTau_boosted_nolepton_charmValue.push_back(islimJet->userFloat("DeepDiTau_boosted_nolepton_charm"));
	       DeepDiTau_boosted_nolepton_charm_massdecoValue.push_back(islimJet->userFloat("DeepDiTau_boosted_nolepton_charm_massdeco"));
	       DeepDiTau_boosted_nolepton_massdecoValue.push_back(islimJet->userFloat("DeepDiTau_boosted_nolepton_massdeco"));
	       DeepDiTau_boosted_noleptonValue.push_back(islimJet->userFloat("DeepDiTau_boosted_nolepton"));
	       DeepDiTau_nolepton_charm_massdecoValue.push_back(islimJet->userFloat("DeepDiTau_nolepton_charm_massdeco"));
	       DeepDiTau_nolepton_charmValue.push_back(islimJet->userFloat("DeepDiTau_nolepton_charm"));
	       DeepDiTau_nolepton_massdecoValue.push_back(islimJet->userFloat("DeepDiTau_nolepton_massdeco"));
	       DeepDiTau_noleptonValue.push_back(islimJet->userFloat("DeepDiTau_nolepton"));
	       //	       std::cout << "jet index: " << jetindex << std::endl;
	       //	       jet_index.push_back(islimJet->userInt(jetindex));
	       jet_index.push_back(jetindex);
	       jet_pt.push_back(islimJet->pt());
	       jet_eta.push_back(islimJet->eta());
	       jet_phi.push_back(islimJet->phi());
	       jet_energy.push_back(islimJet->energy());
	       jet_mass.push_back(islimJet->mass());
	       int nDaughters = islimJet->numberOfDaughters();
	       //	       std:: cout << "number of jet daughters: " << nDaughters << std::endl;
	       for (int iDaughter=0; iDaughter<nDaughters; iDaughter++)
		 {
		   int jetdauId = islimJet->daughter(iDaughter)->pdgId();
		   jet_daughters_index.push_back(jetindex);
		   jet_daughters_pdgid.push_back(islimJet->daughter(iDaughter)->pdgId());
		 }
	       jetindex++;
	       if (isMC){
		 int nHadronic = 0;
		 int nMuon = 0;
		 int nElectron = 0;
		 bool isTau = false;
		 bool isMuon = false;
		 bool isElectron = false;
		 bool isHadronic = false;


		 if(pGenParticles->size() > 0){
		   for(edm::View<reco::GenParticle>::const_iterator iParticle=pGenParticles->begin(); iParticle!=pGenParticles->end(); ++iParticle){		 
		     isTau = false;
		     isMuon = false;
		     isElectron = false;
		     isHadronic = false;
		     if (abs(iParticle->pdgId())==15) {
		       if (reco::deltaR(*islimJet,*iParticle)>0.4) continue;
		      		     // get decay mode
		       for (size_t d=0; d<iParticle->numberOfDaughters(); d++) {
			 if (abs(iParticle->daughter(d)->pdgId())==15) {
			   isTau = true;
			   break;
			 } else if (abs(iParticle->daughter(d)->pdgId())==13) {
			   isMuon = true;
			   break;
			 } else if (abs(iParticle->daughter(d)->pdgId())==11) {
			   isElectron = true;
			   break;
			 } else if (abs(iParticle->daughter(d)->pdgId())!=16 && abs(iParticle->daughter(d)->pdgId())!=14 && abs(iParticle->daughter(d)->pdgId())!=12) {
			 // TODO: can instead start counting neutral and charged to get DM
			   if (iParticle->daughter(d)->pdgId()!=22) {
			     isHadronic = true;
			     break;
			   } // if pdg not 22
			 } // if not neutrino
		       } // over daughter particles

		       if (isTau) {
			 continue; // not last tau before decay
		       } else if (isMuon) {
			 nMuon++;
		       } else if (isElectron) {
			 nElectron++;
		       } else if (isHadronic) {
			 nHadronic++;
		       }
		     } // if pdg is 15
		   } // loop over gen particles
		 } // if gen particles > 0

		 int nTau = nHadronic+nMuon+nElectron;
		 if (nTau==2 || nTau==1) {
		   if(nHadronic==2 && nMuon==0 && nElectron==0){
		     jet_isTauHTauH.push_back(true);
		     DeepDiTauValue_genmatched.push_back(islimJet->userFloat("ditau2017v1"));
		     DeepDiTauValueMD_genmatched.push_back(islimJet->userFloat("ditau2017MDv1"));
		   } else if(nHadronic==1 && nMuon==1 && nElectron==0){
		     jet_isTauHTauM.push_back(true);
		   } else if (nHadronic==1 && nMuon==0 && nElectron==1){
		     jet_isTauHTauE.push_back(true);
		   }else if (nHadronic==0 && nMuon==2 && nElectron==0){
		     jet_isTauMTauM.push_back(true);
		   } else if (nHadronic==0 && nMuon==1 && nElectron==1){
		     jet_isTauMTauE.push_back(true);
		   } else if (nHadronic==0 && nMuon==0 && nElectron==2){
		     jet_isTauETauE.push_back(true);
		   } else if (nHadronic==1 && nMuon==0 && nElectron==0){
		     jet_isTauH.push_back(true);
		   } else if (nHadronic==0 && nMuon==1 && nElectron==0){
		     jet_isTauM.push_back(true);
		   } else if (nHadronic==0 && nMuon==0 && nElectron==1){
		     jet_isTauE.push_back(true);
		   }
		 } // if nTau is 1 or 2
		 else if (islimJet->hadronFlavour()==5) {
		   jet_isB.push_back(true);
		 } else if (islimJet->hadronFlavour()==4) {
		   jet_isC.push_back(true);
		 } else if (std::abs(islimJet->partonFlavour())==3) {
		   jet_isS.push_back(true);
		 } else if (std::abs(islimJet->partonFlavour())==2 || std::abs(islimJet->partonFlavour())==1) {
		   jet_isUD.push_back(true);
		 } else if (islimJet->partonFlavour()==21) {
		   jet_isG.push_back(true);
		 }
	       
		  
	       } // close if MC	     
	      

	       recoJetCSV.push_back(islimJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	     } // pt cut
	   } // eta cut
	 } // loop over slimjets
     } // end if pslimJet->size() 

  
   // --- prepare MET vector ---
   if (pMet->size()>0)
   {
       for(edm::View<pat::MET>::const_iterator iMet=pMet->begin(); iMet!=pMet->end(); iMet++)
       {
           recoMET.push_back(iMet->pt());
           recoMETPhi.push_back(iMet->phi());
           recoMETPx.push_back(iMet->px());
           recoMETPy.push_back(iMet->py());
       } // end for loop on METs
   } // end if pMet->size()>0
  
   // --- fill the object tree ---
   objectTree->Fill();

   // ---- clear all the vectors for next event ----
   // --- reconstructed muons ---
   recoMuonPt.clear();
   recoMuonEta.clear();
   recoMuonPhi.clear();
   recoMuonEnergy.clear();
   recoMuonPDGId.clear();
   recoMuonIsolation.clear();
   recoMuonDXY.clear();
   recoMuonDZ.clear();
   recoMuonNTrackerLayers.clear();
   recoMuonTriggerFlag.clear();
   recoMuonRefToElectron.clear();
   recoMuonRefToTau.clear();

   // --- reconstructed electrons ---
   recoElectronPt.clear();
   recoElectronEta.clear();
   recoElectronPhi.clear();
   recoElectronEnergy.clear();
   recoElectronPDGId.clear();
   recoElectronIsolation.clear();
   recoElectronIdLoose.clear();
   recoElectronIdMedium.clear();
   recoElectronIdTight.clear();
   recoElectronEcalTrkEnergyPostCorr.clear();
   recoElectronEcalTrkEnergyErrPostCorr.clear();
   recoElectronRefToMuon.clear();
   recoElectronRefToTau.clear();

   // --- reconstructed taus ---
   recoTauPt.clear();
   recoTauEta.clear();
   recoTauPhi.clear();
   recoTauEnergy.clear();
   recoTauPDGId.clear();
   recoTauDecayMode.clear();
   recoTauRefToMuon.clear();
   recoTauRefToElectron.clear();

   recoTauIsoMVArawValue.clear();
   recoTauIsoMVAVVLoose.clear();
   recoTauIsoMVAVLoose.clear();
   recoTauIsoMVALoose.clear();
   recoTauIsoMVAMedium.clear();
   recoTauIsoMVATight.clear();
   recoTauIsoMVAVTight.clear();
   recoTauIsoMVAVVTight.clear();
   
   recoTauAntiMuMVALoose.clear();
   recoTauAntiMuMVATight.clear();

   recoTauAntiEleMVArawValue.clear();
   recoTauAntiEleMVAVLoose.clear();
   recoTauAntiEleMVALoose.clear();
   recoTauAntiEleMVAMedium.clear();
   recoTauAntiEleMVATight.clear();
   recoTauAntiEleMVAVTight.clear();

   recoTauDeepVSeraw.clear();
   recoTauDeepVSjetraw.clear();
   recoTauDeepVSmuraw.clear();
 
   recoTauDeepVSeLoose.clear();
   recoTauDeepVSjetLoose.clear();
   recoTauDeepVSmuLoose.clear();

   recoTauDeepVSeMedium.clear();
   recoTauDeepVSjetMedium.clear();
   recoTauDeepVSmuMedium.clear();

   recoTauDeepVSeTight.clear();
   recoTauDeepVSjetTight.clear();
   recoTauDeepVSmuTight.clear();

   recoTauDeepVSeVLoose.clear();
   recoTauDeepVSjetVLoose.clear();
   recoTauDeepVSmuVLoose.clear();

   recoTauDeepVSeVTight.clear();
   recoTauDeepVSjetVTight.clear();

   recoTauDeepVSeVVLoose.clear();
   recoTauDeepVSjetVVLoose.clear();

   recoTauDeepVSeVVTight.clear();
   recoTauDeepVSjetVVTight.clear();

   recoTauDeepVSeVVVLoose.clear();
   recoTauDeepVSjetVVVLoose.clear();
   
   // --- reconstructed jets ---
   recoJetPt.clear();
   recoJetEta.clear();
   recoJetPhi.clear();
   recoJetEnergy.clear();
   recoJetCSV.clear();
   
   /// ---- DeepDiTau jets ---
   oldDeepDiTauValue.clear();
   oldDeepDiTauValueMD.clear();

   DeepDiTauValue_genmatched.clear();
   DeepDiTauValueMD_genmatched.clear();

   DeepDiTau_boosted_nolepton_charmValue.clear();
   DeepDiTau_boosted_nolepton_charm_massdecoValue.clear();
   DeepDiTau_boosted_nolepton_massdecoValue.clear();
   DeepDiTau_boosted_noleptonValue.clear();

   DeepDiTau_nolepton_charm_massdecoValue.clear();
   DeepDiTau_nolepton_charmValue.clear();
   DeepDiTau_nolepton_massdecoValue.clear();
   DeepDiTau_noleptonValue.clear();

   jet_index.clear();
   jet_daughters_index.clear();
   jet_daughters_pdgid.clear();
   
   jet_pt.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_energy.clear();
   jet_mass.clear();
   jet_jetCharge.clear();

   jet_numberOfDaughters.clear();
   jet_chargedMultiplicity.clear();
   jet_neutralMultiplicity.clear();

   jet_chargedHadronMultiplicity.clear();
   jet_neutralHadronMultiplicity.clear();
   jet_muonMultiplicity.clear();
   jet_electronMultiplicity.clear();
   jet_photonMultiplicity.clear();

   jet_chargedEmEnergy.clear();
   jet_neutralEmEnergy.clear();
   jet_chargedHadronEnergy.clear();
   jet_neutralHadronEnergy.clear();
   jet_muonEnergy.clear();
   jet_electronEnergy.clear();
   jet_photonEnergy.clear();

   jet_chargedEmEnergyFraction.clear();
   jet_neutralEmEnergyFraction.clear();
   jet_chargedHadronEnergyFraction.clear();
   jet_neutralHadronEnergyFraction.clear();
   jet_muonEnergyFraction.clear();
   jet_electronEnergyFraction.clear();
   jet_photonEnergyFraction.clear();

   jet_n60.clear();
   jet_n90.clear();

   jet_hadronFlavour.clear();
   jet_partonFlavour.clear();

   jet_pfJetBProbabilityBJetTags.clear();
   jet_pfJetProbabilityBJetTags.clear();
   jet_pfTrackCountingHighEffBJetTags.clear();
   jet_pfSimpleSecondaryVertexHighEffBJetTags.clear();
   jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags.clear();
   jet_pfCombinedSecondaryVertexV2BJetTags.clear();
   jet_pfCombinedInclusiveSecondaryVertexV2BJetTags.clear();
   jet_softPFMuonBJetTags.clear();
   jet_softPFElectronBJetTags.clear();
   jet_pfCombinedMVAV2BJetTags.clear();
   jet_pfCombinedCvsLJetTags.clear();
   jet_pfCombinedCvsBJetTags.clear();
   jet_pfDeepCSVJetTags_probb.clear();
   jet_pfDeepCSVJetTags_probc.clear();
   jet_pfDeepCSVJetTags_probudsg.clear();
   jet_pfDeepCSVJetTags_probbb.clear();

   jet_isTauHTauH.clear();
   jet_isTauHTauM.clear();
   jet_isTauHTauE.clear();
   jet_isTauMTauM.clear();
   jet_isTauMTauE.clear();
   jet_isTauETauE.clear();
   jet_isTauH.clear();
   jet_isTauM.clear();
   jet_isTauE.clear();
   jet_isB.clear();
   jet_isC.clear();
   jet_isS.clear();
   jet_isUD.clear();
   jet_isG.clear();


   /*
   DeepDiTaujetPt.clear();
   DeepDiTaujetEta.clear();
   DeepDiTaujetPhi.clear();
   recoJetCSV.clear(); 
   */


   // --- reconstructed MET ---
   recoMET.clear();
   recoMETPhi.clear();
   recoMETPx.clear();
   recoMETPy.clear();
  
   // ---- gen muons ----
   genMuonPt.clear();
   genMuonEta.clear();
   genMuonPhi.clear();
   genMuonMass.clear();
   genMuonPDGId.clear();
   genMuonMotherPDGId.clear();

   // ---- gen electrons ----
   genElectronPt.clear();
   genElectronEta.clear();
   genElectronPhi.clear();
   genElectronMass.clear();
   genElectronPDGId.clear();
   genElectronMotherPDGId.clear();

   // ---- gen tau_mus ----
   genTauMuPt.clear();
   genTauMuEta.clear();
   genTauMuPhi.clear();
   genTauMuMass.clear();
   genTauMuPDGId.clear();
   genTauMuMotherPDGId.clear();
   genTauMuVisPt.clear();
   genTauMuVisMass.clear();

   // ---- gen tau_es ----
   genTauElePt.clear();
   genTauEleEta.clear();
   genTauElePhi.clear();
   genTauEleMass.clear();
   genTauElePDGId.clear();
   genTauEleMotherPDGId.clear();
   genTauEleVisPt.clear();
   genTauEleVisMass.clear();

   // ---- gen tau_hs ----
   genTauHadPt.clear();
   genTauHadEta.clear();
   genTauHadPhi.clear();
   genTauHadMass.clear();
   genTauHadPDGId.clear();
   genTauHadMotherPDGId.clear();
   genTauHadVisPt.clear();
   genTauHadVisMass.clear();
   genTauHadNPionZero.clear();
   genTauHadNChargedHadrons.clear();
}

// ------------ function for adding up all the visible daughter particles of tau_mu/tau_e ----------------
std::vector<const reco::Candidate*>
DiMuDiTauAnalyzer::findTauMuEleVisDaughters(const reco::Candidate* inputDaughter)
{
    std::vector<const reco::Candidate*> visParticles;
    if (inputDaughter->status() == 1)
    {
        if (fabs(inputDaughter->pdgId()) == 11 || fabs(inputDaughter->pdgId()) == 13 || fabs(inputDaughter->pdgId()) == 22)
        {
            visParticles.push_back(inputDaughter);
        } // end if final state particle is mu/ele/gamma
    } // end if input daughter is final state particle

    else{
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            if (grandDaughter->status() == 1)
            {
                if (fabs(grandDaughter->pdgId()) == 11 || fabs(grandDaughter->pdgId()) == 13 || fabs(grandDaughter->pdgId()) == 22)
                {
                    visParticles.push_back(grandDaughter);
                } // end if final state particle is mu/ele/gamma
            } // end if grand daughter is final state particle

            else{
                auto auxVisParticles = findTauMuEleVisDaughters(grandDaughter);
                visParticles.insert(visParticles.end(), auxVisParticles.begin(), auxVisParticles.end());
            } // end else grand daughter is(not) final state particle
        } // end for loop on grand daughters
    } // end else input daughter is(not) final state particle

    return visParticles;
}

// ------------ function for adding up all the visible daughter particles of tau_h ----------------
std::vector<const reco::Candidate*>
DiMuDiTauAnalyzer::findTauHadVisDaughters(const reco::Candidate* inputDaughter)
{
    std::vector<const reco::Candidate*> visParticles;
    if (inputDaughter->status() == 1)
    {
        if (fabs(inputDaughter->pdgId()) != 12 && fabs(inputDaughter->pdgId()) != 14 && fabs(inputDaughter->pdgId()) != 16)
        {
            visParticles.push_back(inputDaughter);
        } // end if final state particle is not neutrinos
    } // end if input daughter is final state particle

    else{
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            if (grandDaughter->status() == 1)
            {
                if (fabs(grandDaughter->pdgId()) != 12 && fabs(grandDaughter->pdgId()) != 14 && fabs(grandDaughter->pdgId()) != 16)
                {
                    visParticles.push_back(grandDaughter);
                } // end if final state particle is not neutrinos
            } // end if grand daughter is final state particle

            else{
                auto auxVisParticles = findTauHadVisDaughters(grandDaughter);
                visParticles.insert(visParticles.end(), auxVisParticles.begin(), auxVisParticles.end());
            } // end else grand daughter is(not) final state particle
        } // end for loop on grand daughters
    } // end else input daughter is(not) final state particle

    return visParticles;
}

// ------------ function for collecting all the pizeros from tau_h decay ----------------
int DiMuDiTauAnalyzer::findTauPiZeros(const reco::Candidate* inputDaughter)
{
    int numPiZero = 0;
    if (fabs(inputDaughter->pdgId()) == 111) numPiZero++;
    else if (fabs(inputDaughter->pdgId()) == 15)
    {
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            if (fabs(grandDaughter->pdgId()) == 111) numPiZero++;
            else if (fabs(grandDaughter->pdgId()) == 15)
            {
                numPiZero += findTauPiZeros(grandDaughter);
            } // end if the grand-daughter is still tau_h due to FSR
        } // end for loop on the grand-daughters
    } // end if the input daughter is still tau_h because of FSR

    return numPiZero;
}

// ------------ function for collecting all the charged hadrons from tau_h decay ----------------
int DiMuDiTauAnalyzer::findTauChargedHadrons(const reco::Candidate* inputDaughter)
{
    int numChargedHadrons = 0;
    bool chargedHadronsFromTau = fabs(inputDaughter->charge()) != 0 && fabs(inputDaughter->pdgId()) != 11 && fabs(inputDaughter->pdgId()) != 13 && fabs(inputDaughter->pdgId()) != 15; 

    if (chargedHadronsFromTau) numChargedHadrons++;
    else if (fabs(inputDaughter->pdgId()) == 15)
    {
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            chargedHadronsFromTau = fabs(grandDaughter->charge()) != 0 && fabs(grandDaughter->pdgId()) != 11 && fabs(grandDaughter->pdgId()) != 13 && fabs(grandDaughter->pdgId()) != 15;
            if (chargedHadronsFromTau) numChargedHadrons++;
            else if (fabs(grandDaughter->pdgId()) == 15)
            {
                numChargedHadrons += findTauChargedHadrons(grandDaughter);
            } // end if the grand-daughter is still tau_h due to FSR
        } // end for loop on the grand-daughters
    } // end if the input daughter is still tau_h because of FSR

    return numChargedHadrons;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DiMuDiTauAnalyzer::beginJob()
{
    // -- initialize the skimmed object vector tree --
    edm::Service<TFileService> fileService;
    objectTree = fileService->make<TTree>("objectTree","objectTree");

    objectTree->Branch("recoMuonPt", &recoMuonPt);
    objectTree->Branch("recoMuonEta", &recoMuonEta);
    objectTree->Branch("recoMuonPhi", &recoMuonPhi);
    objectTree->Branch("recoMuonEnergy", &recoMuonEnergy);
    objectTree->Branch("recoMuonPDGId", &recoMuonPDGId);
    objectTree->Branch("recoMuonIsolation", &recoMuonIsolation);
    objectTree->Branch("recoMuonDXY", &recoMuonDXY);
    objectTree->Branch("recoMuonDZ", &recoMuonDZ);
    objectTree->Branch("recoMuonNTrackerLayers", &recoMuonNTrackerLayers);
    objectTree->Branch("recoMuonTriggerFlag", &recoMuonTriggerFlag);
    objectTree->Branch("recoMuonRefToElectron", &recoMuonRefToElectron);
    objectTree->Branch("recoMuonRefToTau", &recoMuonRefToTau);

    objectTree->Branch("recoElectronPt", &recoElectronPt);
    objectTree->Branch("recoElectronEta", &recoElectronEta);
    objectTree->Branch("recoElectronPhi", &recoElectronPhi);
    objectTree->Branch("recoElectronEnergy", &recoElectronEnergy);
    objectTree->Branch("recoElectronPDGId", &recoElectronPDGId);
    objectTree->Branch("recoElectronIsolation", &recoElectronIsolation);
    objectTree->Branch("recoElectronIdLoose", &recoElectronIdLoose);
    objectTree->Branch("recoElectronIdMedium", &recoElectronIdMedium);
    objectTree->Branch("recoElectronIdTight", &recoElectronIdTight);
    objectTree->Branch("recoElectronEcalTrkEnergyPostCorr", &recoElectronEcalTrkEnergyPostCorr);
    objectTree->Branch("recoElectronEcalTrkEnergyErrPostCorr", &recoElectronEcalTrkEnergyErrPostCorr);
    objectTree->Branch("recoElectronRefToMuon", &recoElectronRefToMuon);
    objectTree->Branch("recoElectronRefToTau", &recoElectronRefToTau);

    objectTree->Branch("recoTauPt", &recoTauPt);
    objectTree->Branch("recoTauEta", &recoTauEta);
    objectTree->Branch("recoTauPhi", &recoTauPhi);
    objectTree->Branch("recoTauEnergy", &recoTauEnergy);
    objectTree->Branch("recoTauPDGId", &recoTauPDGId);
    objectTree->Branch("recoTauDecayMode", &recoTauDecayMode);
    objectTree->Branch("recoTauRefToMuon", &recoTauRefToMuon);
    objectTree->Branch("recoTauRefToElectron", &recoTauRefToElectron);

    objectTree->Branch("recoTauDeepVSeraw", &recoTauDeepVSeraw);
    objectTree->Branch("recoTauDeepVSjetraw", &recoTauDeepVSjetraw);
    objectTree->Branch("recoTauDeepVSmuraw", &recoTauDeepVSmuraw);

    objectTree->Branch("recoTauDeepVSeLoose", &recoTauDeepVSeLoose);
    objectTree->Branch("recoTauDeepVSjetLoose", &recoTauDeepVSjetLoose);
    objectTree->Branch("recoTauDeepVSmuLoose", &recoTauDeepVSmuLoose);

    objectTree->Branch("recoTauDeepVSeMedium", &recoTauDeepVSeMedium);
    objectTree->Branch("recoTauDeepVSjetMedium", &recoTauDeepVSjetMedium);
    objectTree->Branch("recoTauDeepVSmuMedium", &recoTauDeepVSmuMedium);

    objectTree->Branch("recoTauDeepVSeTight", &recoTauDeepVSeTight);
    objectTree->Branch("recoTauDeepVSjetTight", &recoTauDeepVSjetTight);
    objectTree->Branch("recoTauDeepVSmuTight", &recoTauDeepVSmuTight);

    objectTree->Branch("recoTauDeepVSeVLoose", &recoTauDeepVSeVLoose);
    objectTree->Branch("recoTauDeepVSjetVLoose", &recoTauDeepVSjetVLoose);
    objectTree->Branch("recoTauDeepVSmuVLoose", &recoTauDeepVSmuVLoose);

    objectTree->Branch("recoTauDeepVSeVTight", &recoTauDeepVSeVTight);
    objectTree->Branch("recoTauDeepVSjetVTight", &recoTauDeepVSjetVTight);

    objectTree->Branch("recoTauDeepVSeVVLoose", &recoTauDeepVSeVVLoose);
    objectTree->Branch("recoTauDeepVSjetVVLoose", &recoTauDeepVSjetVVLoose);

    objectTree->Branch("recoTauDeepVSeVVTight", &recoTauDeepVSeVVTight);
    objectTree->Branch("recoTauDeepVSjetVVTight", &recoTauDeepVSjetVVTight);

    objectTree->Branch("recoTauDeepVSeVVVLoose", &recoTauDeepVSeVVVLoose);
    objectTree->Branch("recoTauDeepVSjetVVVLoose", &recoTauDeepVSjetVVVLoose);

    objectTree->Branch("recoTauIsoMVArawValue", &recoTauIsoMVArawValue);
    objectTree->Branch("recoTauIsoMVAVVLoose", &recoTauIsoMVAVVLoose);
    objectTree->Branch("recoTauIsoMVAVLoose", &recoTauIsoMVAVLoose);
    objectTree->Branch("recoTauIsoMVALoose", &recoTauIsoMVALoose);
    objectTree->Branch("recoTauIsoMVAMedium", &recoTauIsoMVAMedium);
    objectTree->Branch("recoTauIsoMVATight", &recoTauIsoMVATight);
    objectTree->Branch("recoTauIsoMVAVTight", &recoTauIsoMVAVTight);
    objectTree->Branch("recoTauIsoMVAVVTight", &recoTauIsoMVAVVTight);

    objectTree->Branch("recoTauAntiMuMVALoose", &recoTauAntiMuMVALoose);
    objectTree->Branch("recoTauAntiMuMVATight", &recoTauAntiMuMVATight);
    
    objectTree->Branch("recoTauAntiEleMVArawValue", &recoTauAntiEleMVArawValue);
    objectTree->Branch("recoTauAntiEleMVAVLoose", &recoTauAntiEleMVAVLoose);
    objectTree->Branch("recoTauAntiEleMVALoose", &recoTauAntiEleMVALoose);
    objectTree->Branch("recoTauAntiEleMVAMedium", &recoTauAntiEleMVAMedium);
    objectTree->Branch("recoTauAntiEleMVATight", &recoTauAntiEleMVATight);
    objectTree->Branch("recoTauAntiEleMVAVTight", &recoTauAntiEleMVAVTight);
    
    objectTree->Branch("recoJetPt", &recoJetPt);
    objectTree->Branch("recoJetEta", &recoJetEta);
    objectTree->Branch("recoJetPhi", &recoJetPhi);
    objectTree->Branch("recoJetEnergy", &recoJetEnergy);
    objectTree->Branch("recoJetCSV", &recoJetCSV);
    
    objectTree->Branch("oldDeepDiTauValue", &oldDeepDiTauValue);
    objectTree->Branch("oldDeepDiTauValueMD", &oldDeepDiTauValueMD);
    /*

    objectTree->Branch("jet_index", &jet_index);
    objectTree->Branch("jet_daughters_pdgid", &jet_daughters_pdgid);
    */
    objectTree->Branch("jet_pt", &jet_pt);
    objectTree->Branch("jet_eta", &jet_eta);
    objectTree->Branch("jet_phi", &jet_phi);
    objectTree->Branch("jet_energy", &jet_energy);
    objectTree->Branch("jet_jetCharge", &jet_jetCharge);

    objectTree->Branch("jet_chargedMultiplicity", &jet_chargedMultiplicity);
    objectTree->Branch("jet_neutralMultiplicity", &jet_neutralMultiplicity);

    objectTree->Branch("jet_chargedHadronMultiplicity", &jet_chargedHadronMultiplicity);
    objectTree->Branch("jet_neutralHadronMultiplicity", &jet_neutralHadronMultiplicity);
    objectTree->Branch("jet_muonMultiplicity", &jet_muonMultiplicity);
    objectTree->Branch("jet_electronMultiplicity", &jet_electronMultiplicity);
    objectTree->Branch("jet_photonMultiplicity", &jet_photonMultiplicity);

    objectTree->Branch("jet_chargedEmEnergy", &jet_chargedEmEnergy);
    objectTree->Branch("jet_neutralEmEnergy", &jet_neutralEmEnergy);
    objectTree->Branch("jet_chargedHadronEnergy", &jet_chargedHadronEnergy);
    objectTree->Branch("jet_neutralHadronEnergy", &jet_neutralHadronEnergy);
    objectTree->Branch("jet_muonEnergy", &jet_muonEnergy);
    objectTree->Branch("jet_electronEnergy", &jet_electronEnergy);
    objectTree->Branch("jet_photonEnergy", &jet_photonEnergy);

    objectTree->Branch("jet_chargedEmEnergyFraction", &jet_chargedEmEnergyFraction);
    objectTree->Branch("jet_neutralEmEnergyFraction", &jet_neutralEmEnergyFraction);
    objectTree->Branch("jet_chargedHadronEnergyFraction", &jet_chargedHadronEnergyFraction);
    objectTree->Branch("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction);
    objectTree->Branch("jet_muonEnergyFraction", &jet_muonEnergyFraction);
    objectTree->Branch("jet_electronEnergyFraction", &jet_electronEnergyFraction);
    objectTree->Branch("jet_photonEnergyFraction", &jet_photonEnergyFraction);

    objectTree->Branch("jet_n60", &jet_n60);
    objectTree->Branch("jet_n90", &jet_n90);

    objectTree->Branch("jet_pfJetBProbabilityBJetTags", &jet_pfJetBProbabilityBJetTags);
    objectTree->Branch("jet_pfJetProbabilityBJetTags", &jet_pfJetProbabilityBJetTags);
    objectTree->Branch("jet_pfTrackCountingHighEffBJetTags", &jet_pfTrackCountingHighEffBJetTags);
    objectTree->Branch("jet_pfSimpleSecondaryVertexHighEffBJetTags", &jet_pfSimpleSecondaryVertexHighEffBJetTags);
    objectTree->Branch("jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags", &jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags);
    objectTree->Branch("jet_pfCombinedSecondaryVertexV2BJetTags", &jet_pfCombinedSecondaryVertexV2BJetTags);
    objectTree->Branch("jet_pfCombinedInclusiveSecondaryVertexV2BJetTags", &jet_pfCombinedInclusiveSecondaryVertexV2BJetTags);
    objectTree->Branch("jet_softPFMuonBJetTags", &jet_softPFMuonBJetTags);
    objectTree->Branch("jet_softPFElectronBJetTags", &jet_softPFElectronBJetTags);
    objectTree->Branch("jet_pfCombinedMVAV2BJetTags", &jet_pfCombinedMVAV2BJetTags);
    objectTree->Branch("jet_pfCombinedCvsLJetTags", &jet_pfCombinedCvsLJetTags);
    objectTree->Branch("jet_pfCombinedCvsBJetTags", &jet_pfCombinedCvsBJetTags);
    objectTree->Branch("jet_pfDeepCSVJetTags_probb", &jet_pfDeepCSVJetTags_probb);
    objectTree->Branch("jet_pfDeepCSVJetTags_probc", &jet_pfDeepCSVJetTags_probc);
    objectTree->Branch("jet_pfDeepCSVJetTags_probudsg", &jet_pfDeepCSVJetTags_probudsg);
    objectTree->Branch("jet_pfDeepCSVJetTags_probbb", &jet_pfDeepCSVJetTags_probbb);

    objectTree->Branch("jet_daughter_pt", &jet_daughter_pt);
    objectTree->Branch("jet_daughter_eta", &jet_daughter_eta);
    objectTree->Branch("jet_daughter_phi", &jet_daughter_phi);
    objectTree->Branch("jet_daughter_energy", &jet_daughter_energy);
    objectTree->Branch("jet_daughter_mass", &jet_daughter_mass);
    objectTree->Branch("jet_daughter_charge", &jet_daughter_charge);

    objectTree->Branch("jet_daughter_etaAtVtx", &jet_daughter_etaAtVtx);
    objectTree->Branch("jet_daughter_phiAtVtx", &jet_daughter_phiAtVtx);
    objectTree->Branch("jet_daughter_vx", &jet_daughter_vx);
    objectTree->Branch("jet_daughter_vy", &jet_daughter_vy);
    objectTree->Branch("jet_daughter_vz", &jet_daughter_vz);
    objectTree->Branch("jet_daughter_dxy", &jet_daughter_dxy);
    objectTree->Branch("jet_daughter_dxyError", &jet_daughter_dxyError);
    objectTree->Branch("jet_daughter_dz", &jet_daughter_dz);
    objectTree->Branch("jet_daughter_dzError", &jet_daughter_dzError);

    objectTree->Branch("jet_daughter_pdgId", &jet_daughter_pdgId);

    objectTree->Branch("jet_daughter_puppiWeight", &jet_daughter_puppiWeight);
    objectTree->Branch("jet_daughter_puppiWeightNoLep", &jet_daughter_puppiWeightNoLep);

    objectTree->Branch("jet_daughter_isIsolatedChargedHadron", &jet_daughter_isIsolatedChargedHadron);
    objectTree->Branch("jet_daughter_isStandAloneMuon", &jet_daughter_isStandAloneMuon);
    objectTree->Branch("jet_daughter_isTrackerMuon", &jet_daughter_isTrackerMuon);
    objectTree->Branch("jet_daughter_isGlobalMuon", &jet_daughter_isGlobalMuon);
    objectTree->Branch("jet_daughter_isGoodEgamma", &jet_daughter_isGoodEgamma);



    /*
    objectTree->Branch("DeepDiTaujetPt", &DeepDiTaujetPt);
    objectTree->Branch("DeepDiTaujetEta", &DeepDiTaujetEta);
    objectTree->Branch("DeepDiTaujetPhi", &DeepDiTaujetPhi);
    objectTree->Branch("DeepDiTaujetEnergy", &DeepDiTaujetEnergy);
    */    
objectTree->Branch("recoJetCSV", &recoJetCSV); 
    objectTree->Branch("recoMET", &recoMET);
    objectTree->Branch("recoMETPhi", &recoMETPhi);
    objectTree->Branch("recoMETPx", &recoMETPx);
    objectTree->Branch("recoMETPy", &recoMETPy);
    objectTree->Branch("recoNPrimaryVertex", &recoNPrimaryVertex, "recoNPrimaryVertex/I");
    objectTree->Branch("eventID", &eventID, "eventID/I");

    if (isMC)
    {
        objectTree->Branch("genMuonPt", &genMuonPt);
        objectTree->Branch("genMuonEta", &genMuonEta);
        objectTree->Branch("genMuonPhi", &genMuonPhi);
        objectTree->Branch("genMuonMass", &genMuonMass);
        objectTree->Branch("genMuonPDGId", &genMuonPDGId);
        objectTree->Branch("genMuonMotherPDGId", &genMuonMotherPDGId);

        objectTree->Branch("genElectronPt", &genElectronPt);
        objectTree->Branch("genElectronEta", &genElectronEta);
        objectTree->Branch("genElectronPhi", &genElectronPhi);
        objectTree->Branch("genElectronMass", &genElectronMass);
        objectTree->Branch("genElectronPDGId", &genElectronPDGId);
        objectTree->Branch("genElectronMotherPDGId", &genElectronMotherPDGId);

        objectTree->Branch("genTauMuPt", &genTauMuPt);
        objectTree->Branch("genTauMuEta", &genTauMuEta);
        objectTree->Branch("genTauMuPhi", &genTauMuPhi);
        objectTree->Branch("genTauMuMass", &genTauMuMass);
        objectTree->Branch("genTauMuPDGId", &genTauMuPDGId);
        objectTree->Branch("genTauMuMotherPDGId", &genTauMuMotherPDGId);
        objectTree->Branch("genTauMuVisPt", &genTauMuVisPt);
        objectTree->Branch("genTauMuVisMass", &genTauMuVisMass);

        objectTree->Branch("genTauElePt", &genTauElePt);
        objectTree->Branch("genTauEleEta", &genTauEleEta);
        objectTree->Branch("genTauElePhi", &genTauElePhi);
        objectTree->Branch("genTauEleMass", &genTauEleMass);
        objectTree->Branch("genTauElePDGId", &genTauElePDGId);
        objectTree->Branch("genTauEleMotherPDGId", &genTauEleMotherPDGId);
        objectTree->Branch("genTauEleVisPt", &genTauEleVisPt);
        objectTree->Branch("genTauEleVisMass", &genTauEleVisMass);

        objectTree->Branch("genTauHadPt", &genTauHadPt);
        objectTree->Branch("genTauHadEta", &genTauHadEta);
        objectTree->Branch("genTauHadPhi", &genTauHadPhi);
        objectTree->Branch("genTauHadMass", &genTauHadMass);
        objectTree->Branch("genTauHadPDGId", &genTauHadPDGId);
        objectTree->Branch("genTauHadMotherPDGId", &genTauHadMotherPDGId);
        objectTree->Branch("genTauHadVisPt", &genTauHadVisPt);
        objectTree->Branch("genTauHadVisMass", &genTauHadVisMass);
        objectTree->Branch("genTauHadNPionZero", &genTauHadNPionZero);
        objectTree->Branch("genTauHadNChargedHadrons", &genTauHadNChargedHadrons);

	objectTree->Branch("DeepDiTau_boosted_nolepton_charmValue", &DeepDiTau_boosted_nolepton_charmValue);
	objectTree->Branch("DeepDiTau_boosted_nolepton_charm_massdecoValue", &DeepDiTau_boosted_nolepton_charm_massdecoValue);
	objectTree->Branch("DeepDiTau_boosted_nolepton_massdecoValue", &DeepDiTau_boosted_nolepton_massdecoValue);
	objectTree->Branch("DeepDiTau_boosted_noleptonValue", &DeepDiTau_boosted_noleptonValue);

	objectTree->Branch("DeepDiTau_nolepton_charm_massdecoValue", &DeepDiTau_nolepton_charm_massdecoValue);
	objectTree->Branch("DeepDiTau_nolepton_charmValue", &DeepDiTau_nolepton_charmValue);
	objectTree->Branch("DeepDiTau_nolepton_massdecoValue", &DeepDiTau_nolepton_massdecoValue);
	objectTree->Branch("DeepDiTau_noleptonValue", &DeepDiTau_noleptonValue);

	objectTree->Branch("DeepDiTauValue_genmatched", &DeepDiTauValue_genmatched);
	objectTree->Branch("DeepDiTauValueMD_genmatched", &DeepDiTauValueMD_genmatched);

	objectTree->Branch("jet_index", &jet_index);
	objectTree->Branch("jet_daughters_index", &jet_daughters_index);
	objectTree->Branch("jet_daughters_pdgid", &jet_daughters_pdgid);

	objectTree->Branch("jet_isTauHTauH", &jet_isTauHTauH);
	objectTree->Branch("jet_isTauHTauM", &jet_isTauHTauM);
	objectTree->Branch("jet_isTauHTauE", &jet_isTauHTauE);
	objectTree->Branch("jet_isTauMTauM", &jet_isTauMTauM);
	objectTree->Branch("jet_isTauMTauE", &jet_isTauMTauE);
	objectTree->Branch("jet_isTauETauE", &jet_isTauETauE);
	objectTree->Branch("jet_isTauH", &jet_isTauH);
	objectTree->Branch("jet_isTauM", &jet_isTauM);
	objectTree->Branch("jet_isTauE", &jet_isTauE);
	objectTree->Branch("jet_isB", &jet_isB);
	objectTree->Branch("jet_isC", &jet_isC);
	objectTree->Branch("jet_isS", &jet_isS);
	objectTree->Branch("jet_isUD", &jet_isUD);
	objectTree->Branch("jet_isG", &jet_isG);

        objectTree->Branch("recoNPU", &recoNPU, "recoNPU/I");
        objectTree->Branch("trueNInteraction", &trueNInteraction, "trueNInteraction/I");
        objectTree->Branch("genEventWeight", &genEventWeight, "genEventWeight/F");
    } // end if isMC == true
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiMuDiTauAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuDiTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuDiTauAnalyzer);
