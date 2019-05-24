// -*- C++ -*-
//
// Package:    testED/MuMuTauMuTauHadAnalyzer
// Class:      MuMuTauMuTauHadAnalyzer
// 
/**\class MuMuTauMuTauHadAnalyzer MuMuTauMuTauHadAnalyzer.cc testED/MuMuTauMuTauHadAnalyzer/plugins/MuMuTauMuTauHadAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Fri, 24 May 2019 17:15:41 GMT
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
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "TTree.h"
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

class MuMuTauMuTauHadAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTauMuTauHadAnalyzer(const edm::ParameterSet&);
      ~MuMuTauMuTauHadAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu1Mu2Tag;
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu3Tag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauTag;
      edm::EDGetTokenT<edm::View<pat::Jet>> JetTag;
      edm::EDGetTokenT<edm::View<pat::Photon>> PhotonTag;
      edm::EDGetTokenT<edm::View<reco::Vertex>> VertexTag;
      bool isMC;
      float EventWeight;
      float NPVertex;
      edm::EDGetTokenT<GenEventInfoProduct> generator_;

      TTree *objectTree;
      vector<pat::Muon> recoMuons;
      vector<pat::Tau> recoTaus;
      vector<pat::Jet> recoJets;
      vector<pat::Photon> recoPhotons;

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
MuMuTauMuTauHadAnalyzer::MuMuTauMuTauHadAnalyzer(const edm::ParameterSet& iConfig):
    Mu1Mu2Tag(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2Tag"))),
    Mu3Tag(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu3Tag"))),
    TauTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauTag"))),
    JetTag(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
    PhotonTag(consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("PhotonTag"))),
    VertexTag(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexTag"))),
    generator_(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag()))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC = iConfig.getParameter<bool>("isMC");
}


MuMuTauMuTauHadAnalyzer::~MuMuTauMuTauHadAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTauMuTauHadAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2Tag, pMu1Mu2);

   edm::Handle<edm::View<pat::Muon>> pMu3;
   iEvent.getByToken(Mu3Tag, pMu3);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauTag, pTau);

   edm::Handle<edm::View<pat::Photon>> pPhoton;
   iEvent.getByToken(PhotonTag, pPhoton);

   edm::Handle<edm::View<pat::Jet>> pJet;
   iEvent.getByToken(JetTag, pJet);

   edm::Handle<edm::View<reco::Vertex>> pVertex;
   iEvent.getByToken(VertexTag, pVertex);

   if (isMC)
   {
       EventWeight = 1.0;
       edm::Handle<GenEventInfoProduct> gen_ev_info;
       iEvent.getByToken(generator_, gen_ev_info);
       if(gen_ev_info.isValid()) EventWeight = gen_ev_info->weight();
   }

   else{
       EventWeight = 1.0;
   }

   NPVertex = 0;
   if (pVertex.isValid())
   {
       for(edm::View<reco::Vertex>::const_iterator iPV=pVertex->begin(); iPV!=pVertex->end(); iPV++)
       {
           NPVertex++;
       }
   }
 
   // --- prepare muon vector ---
   recoMuons.push_back(pMu1Mu2->at(0));
   recoMuons.push_back(pMu1Mu2->at(1));
   
   for(edm::View<pat::Muon>::const_iterator iMuon=pMu3->begin(); iMuon!=pMu3->end(); iMuon++)
   {
       recoMuons.push_back(*iMuon);
   }

   // --- prepare tau vector ---
   for(edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
   {
       recoTaus.push_back(*iTau);
   }

   // --- prepare jet vector ---
   if(pJet->size()>0)
   {
       for(edm::View<pat::Jet>::const_iterator iJet=pJet->begin(); iJet!=pJet->end(); iJet++)
       {
           recoJets.push_back(*iJet);
       }
   }

   // --- prepare photon vector ---
   if(pPhoton->size()>0)
   {
       for(edm::View<pat::Photon>::const_iterator iPhoton=pPhoton->begin(); iPhoton!=pPhoton->end(); iPhoton++)
       {
           recoPhotons.push_back(*iPhoton);
       }
   }

   objectTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTauMuTauHadAnalyzer::beginJob()
{
    // -- initialize the skimmed object vector tree --
    edm::Service<TFileService> fileService;
    objectTree = fileService->make<TTree>("objectTree","objectTree");

    objectTree->Branch("recoMuons", &recoMuons);
    objectTree->Branch("recoTaus", &recoTaus);
    objectTree->Branch("recoJets", &recoJets);
    objectTree->Branch("recoPhotons", &recoPhotons);

    objectTree->Branch("EventWeight", &EventWeight, "EventWeight/F");
    objectTree->Branch("NPVertex", &NPVertex, "NPVertex/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTauMuTauHadAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTauMuTauHadAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTauMuTauHadAnalyzer);
