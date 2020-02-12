// -*- C++ -*-
//
// Package:    MuMuChannel/ZMuMuInclusiveAnalyzer
// Class:      ZMuMuInclusiveAnalyzer
// 
/**\class ZMuMuInclusiveAnalyzer ZMuMuInclusiveAnalyzer.cc MuMuChannel/ZMuMuInclusiveAnalyzer/plugins/ZMuMuInclusiveAnalyzer.cc

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

class ZMuMuInclusiveAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZMuMuInclusiveAnalyzer(const edm::ParameterSet&);
      ~ZMuMuInclusiveAnalyzer();

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

      // --- reconstructed electrons ---
      vector<float> recoElectronPt;
      vector<float> recoElectronEta;
      vector<float> recoElectronPhi;
      vector<float> recoElectronEnergy;
      vector<int> recoElectronPDGId;
      vector<float> recoElectronIsolation;
      vector<float> recoElectronEcalTrkEnergyPostCorr;
      vector<float> recoElectronEcalTrkEnergyErrPostCorr;

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
ZMuMuInclusiveAnalyzer::ZMuMuInclusiveAnalyzer(const edm::ParameterSet& iConfig):
    MuTag(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuTag"))),
    EleTag(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("EleTag"))),
    TauTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauTag"))),
    JetTag(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
    MetTag(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("MetTag"))),
    VertexTag(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexTag"))),
    rhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"))),
    effectiveAreas((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
    PileupTag(consumes<edm::View<PileupSummaryInfo>>(iConfig.existsAs<edm::InputTag>("PileupTag") ? iConfig.getParameter<edm::InputTag>("PileupTag") : edm::InputTag())),
    generator(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag())),
    GenMuTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenMuTag") ? iConfig.getParameter<edm::InputTag>("GenMuTag") : edm::InputTag())),
    GenEleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenEleTag") ? iConfig.getParameter<edm::InputTag>("GenEleTag") : edm::InputTag())),
    GenTauMuTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauMuTag") ? iConfig.getParameter<edm::InputTag>("GenTauMuTag") : edm::InputTag())),
    GenTauEleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauEleTag") ? iConfig.getParameter<edm::InputTag>("GenTauEleTag") : edm::InputTag())),
    GenTauHadTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauHadTag") ? iConfig.getParameter<edm::InputTag>("GenTauHadTag") : edm::InputTag()))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC = iConfig.getParameter<bool>("isMC");
}


ZMuMuInclusiveAnalyzer::~ZMuMuInclusiveAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZMuMuInclusiveAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu;
   iEvent.getByToken(MuTag, pMu);

   edm::Handle<edm::View<pat::Electron>> pElectron;
   iEvent.getByToken(EleTag, pElectron);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauTag, pTau);

   edm::Handle<edm::View<pat::Jet>> pJet;
   iEvent.getByToken(JetTag, pJet);

   edm::Handle<edm::View<pat::MET>> pMet;
   iEvent.getByToken(MetTag, pMet);

   edm::Handle<edm::View<reco::Vertex>> pVertex;
   iEvent.getByToken(VertexTag, pVertex);

   edm::Handle<double> pRho;
   iEvent.getByToken(rhoTag, pRho);

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
                   } // end for loop on all generations of visible daughter particles of tau_h
               } // end for loop on tau_h direct daughter particles

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
           // --- variables for relIsoWithEffectiveArea ---
           float chad = iElectron->pfIsolationVariables().sumChargedHadronPt;
           float nhad = iElectron->pfIsolationVariables().sumNeutralHadronEt;
           float pho = iElectron->pfIsolationVariables().sumPhotonEt;
           float elePt = iElectron->pt();
           float eleEta = iElectron->superCluster()->eta();
           float rho = pRho.isValid() ? (*pRho) : 0;
           float eArea = effectiveAreas.getEffectiveArea(fabs(eleEta));
           float relIsoWithEffectiveArea = (chad + std::max(0.0f, nhad + pho - rho*eArea)) / elePt;

           recoElectronPt.push_back(iElectron->pt());
           recoElectronEta.push_back(iElectron->superCluster()->eta());
           recoElectronPhi.push_back(iElectron->superCluster()->phi());
           recoElectronEnergy.push_back(iElectron->superCluster()->energy());
           recoElectronPDGId.push_back(iElectron->pdgId());
           recoElectronIsolation.push_back(relIsoWithEffectiveArea);
           recoElectronEcalTrkEnergyPostCorr.push_back(iElectron->userFloat("ecalTrkEnergyPostCorr"));
           recoElectronEcalTrkEnergyErrPostCorr.push_back(iElectron->userFloat("ecalTrkEnergyErrPostCorr"));
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
           recoTauIsoMVArawValue.push_back(iTau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
           recoTauIsoMVAVVLoose.push_back(iTau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
           recoTauIsoMVAVLoose.push_back(iTau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
           recoTauIsoMVALoose.push_back(iTau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
           recoTauIsoMVAMedium.push_back(iTau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
           recoTauIsoMVATight.push_back(iTau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
           recoTauIsoMVAVTight.push_back(iTau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
           recoTauIsoMVAVVTight.push_back(iTau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
       
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

   // --- prepare jet vector ---
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

   // --- reconstructed electrons ---
   recoElectronPt.clear();
   recoElectronEta.clear();
   recoElectronPhi.clear();
   recoElectronEnergy.clear();
   recoElectronPDGId.clear();
   recoElectronIsolation.clear();
   recoElectronEcalTrkEnergyPostCorr.clear();
   recoElectronEcalTrkEnergyErrPostCorr.clear();

   // --- reconstructed taus ---
   recoTauPt.clear();
   recoTauEta.clear();
   recoTauPhi.clear();
   recoTauEnergy.clear();
   recoTauPDGId.clear();
   recoTauDecayMode.clear();
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

   // --- reconstructed jets ---
   recoJetPt.clear();
   recoJetEta.clear();
   recoJetPhi.clear();
   recoJetEnergy.clear();
   recoJetCSV.clear();
   
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
ZMuMuInclusiveAnalyzer::findTauMuEleVisDaughters(const reco::Candidate* inputDaughter)
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
ZMuMuInclusiveAnalyzer::findTauHadVisDaughters(const reco::Candidate* inputDaughter)
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
int ZMuMuInclusiveAnalyzer::findTauPiZeros(const reco::Candidate* inputDaughter)
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
int ZMuMuInclusiveAnalyzer::findTauChargedHadrons(const reco::Candidate* inputDaughter)
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
ZMuMuInclusiveAnalyzer::beginJob()
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

    objectTree->Branch("recoElectronPt", &recoElectronPt);
    objectTree->Branch("recoElectronEta", &recoElectronEta);
    objectTree->Branch("recoElectronPhi", &recoElectronPhi);
    objectTree->Branch("recoElectronEnergy", &recoElectronEnergy);
    objectTree->Branch("recoElectronPDGId", &recoElectronPDGId);
    objectTree->Branch("recoElectronIsolation", &recoElectronIsolation);
    objectTree->Branch("recoElectronEcalTrkEnergyPostCorr", &recoElectronEcalTrkEnergyPostCorr);
    objectTree->Branch("recoElectronEcalTrkEnergyErrPostCorr", &recoElectronEcalTrkEnergyErrPostCorr);

    objectTree->Branch("recoTauPt", &recoTauPt);
    objectTree->Branch("recoTauEta", &recoTauEta);
    objectTree->Branch("recoTauPhi", &recoTauPhi);
    objectTree->Branch("recoTauEnergy", &recoTauEnergy);
    objectTree->Branch("recoTauPDGId", &recoTauPDGId);
    objectTree->Branch("recoTauDecayMode", &recoTauDecayMode);
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

        objectTree->Branch("recoNPU", &recoNPU, "recoNPU/I");
        objectTree->Branch("trueNInteraction", &trueNInteraction, "trueNInteraction/I");
        objectTree->Branch("genEventWeight", &genEventWeight, "genEventWeight/F");
    } // end if isMC == true
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZMuMuInclusiveAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZMuMuInclusiveAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZMuMuInclusiveAnalyzer);
