// -*- C++ -*-
//
// Package:    MuMuChannel/MuMuTauTauAnalyzer
// Class:      MuMuTauTauAnalyzer
// 
/**\class MuMuTauTauAnalyzer MuMuTauTauAnalyzer.cc MuMuChannel/MuMuTauTauAnalyzer/plugins/MuMuTauTauAnalyzer.cc

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
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "TTree.h"
#include <math.h>
#include <string>
#include <iostream>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuMuTauTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTauTauAnalyzer(const edm::ParameterSet&);
      ~MuMuTauTauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu1_;
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu2Mu3_;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauHad_;
      std::vector<std::string> tauDiscriminatorTags_;
      edm::EDGetTokenT<edm::View<reco::Vertex>> Vertex_;
      bool isMC_;
      float EventWeight;
      float NPVertex;
      edm::EDGetTokenT<GenEventInfoProduct> generator_;

      TTree *MuMuTauTauTree;
      double Mu1Energy;
      double Mu2Energy;
      double Mu3Energy;
      double tauHadEnergy;
      double Mu1Pt;
      double Mu2Pt;
      double Mu3Pt;
      double tauHadPt;
      double Mu1Eta;
      double Mu2Eta;
      double Mu3Eta;
      double tauHadEta;
      double Mu1Phi;
      double Mu2Phi;
      double Mu3Phi;
      double tauHadPhi;
      double Mu1Iso;
      double Mu2Iso;
      double Mu3Iso;
      // ---- below is variables for Tau discriminators ---
      double tauIsolationMVArawValue;
      double tauIsolationMVAVVLoose;
      double tauIsolationMVAVLoose;
      double tauIsolationMVALoose;
      double tauIsolationMVAMedium;
      double tauIsolationMVATight;
      double tauIsolationMVAVTight;
      double tauIsolationMVAVVTight;
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
MuMuTauTauAnalyzer::MuMuTauTauAnalyzer(const edm::ParameterSet& iConfig):
    Mu1_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu1Tag"))),
    Mu2Mu3_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu2Mu3Tag"))),
    TauHad_(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauHadTag"))),
    Vertex_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("Vertex"))),
    generator_(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag()))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC_ = iConfig.getParameter<bool>("isMC");
   tauDiscriminatorTags_ = iConfig.getParameter<std::vector<std::string>>("tauDiscriminatorTags");
   tauIsolationMVArawValue = -999.0;
   tauIsolationMVAVVLoose = -999.0;
   tauIsolationMVAVLoose = -999.0;
   tauIsolationMVALoose = -999.0;
   tauIsolationMVAMedium = -999.0;
   tauIsolationMVATight = -999.0;
   tauIsolationMVAVTight = -999.0;
   tauIsolationMVAVVTight = -999.0;
}


MuMuTauTauAnalyzer::~MuMuTauTauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTauTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu1;
   iEvent.getByToken(Mu1_, pMu1);

   edm::Handle<edm::View<pat::Muon>> pMu2Mu3;
   iEvent.getByToken(Mu2Mu3_, pMu2Mu3);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauHad_, pTau);

   edm::Handle<edm::View<reco::Vertex>> pVertex;
   iEvent.getByToken(Vertex_, pVertex);

   if (isMC_)
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

   pat::Muon Mu1 = pMu1->at(0);
   pat::Muon Mu2 = pMu2Mu3->at(0);
   pat::Muon Mu3 = pMu2Mu3->at(1);

   Mu1Energy = Mu1.energy();
   Mu1Pt = Mu1.pt();
   Mu1Eta = Mu1.eta();
   Mu1Phi = Mu1.phi();

   Mu2Energy = Mu2.energy();
   Mu2Pt = Mu2.pt();
   Mu2Eta = Mu2.eta();
   Mu2Phi = Mu2.phi();

   Mu3Energy = Mu3.energy();
   Mu3Pt = Mu3.pt();
   Mu3Eta = Mu3.eta();
   Mu3Phi = Mu3.phi();

   reco::MuonPFIsolation iso = Mu1.pfIsolationR04();
   Mu1Iso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / Mu1.pt();

   iso = Mu2.pfIsolationR04();
   Mu2Iso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / Mu2.pt();

   iso = Mu3.pfIsolationR04();
   Mu3Iso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / Mu3.pt();

   // select the tau with the smallest dR from mu3
   double smallestDR = 99.0;
   pat::Tau tauHad;
   for(edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); ++iTau)
   {
       if(deltaR(Mu3, *iTau) < smallestDR)
       {
           smallestDR = deltaR(Mu3, *iTau);
           tauHad = *iTau;
       }
   }

   tauHadEnergy = tauHad.energy();
   tauHadPt = tauHad.pt();
   tauHadEta = tauHad.eta();
   tauHadPhi = tauHad.phi();

   for (unsigned int i=0; i<tauDiscriminatorTags_.size(); i++)
   {
       if(tauDiscriminatorTags_[i].find("raw")!=std::string::npos)
       {
           tauIsolationMVArawValue = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("VVLoose")!=std::string::npos)
       {
           tauIsolationMVAVVLoose = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("VLoose")!=std::string::npos && tauDiscriminatorTags_[i].find("VVLoose")==std::string::npos)
       {
           tauIsolationMVAVLoose = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("Loose")!=std::string::npos && tauDiscriminatorTags_[i].find("VLoose")==std::string::npos && tauDiscriminatorTags_[i].find("VVLoose")==std::string::npos)
       {
           tauIsolationMVALoose = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("Medium")!=std::string::npos)
       {
           tauIsolationMVAMedium = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("Tight")!=std::string::npos && tauDiscriminatorTags_[i].find("VTight")==std::string::npos && tauDiscriminatorTags_[i].find("VVTight")==std::string::npos)
       {
           tauIsolationMVATight = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("VTight")!=std::string::npos && tauDiscriminatorTags_[i].find("VVTight")==std::string::npos)
       {
           tauIsolationMVAVTight = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else if(tauDiscriminatorTags_[i].find("VVTight")!=std::string::npos)
       {
           tauIsolationMVAVVTight = tauHad.tauID(tauDiscriminatorTags_[i]); 
       }

       else{
           std::cout << "Unknown category of tau discriminator: " << tauDiscriminatorTags_[i] << std::endl;
           std::cout << "Please make sure you are using the right name." << std::endl;
           continue;
       }
   }

   MuMuTauTauTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTauTauAnalyzer::beginJob()
{
    // -- initialize the skimmed flat tree --
    edm::Service<TFileService> fileService;
    MuMuTauTauTree = fileService->make<TTree>("MuMuTauTauTree","MuMuTauTauTree");
    
    MuMuTauTauTree->Branch("Mu1Energy", &Mu1Energy, "Mu1Energy/D");
    MuMuTauTauTree->Branch("Mu2Energy", &Mu2Energy, "Mu2Energy/D");
    MuMuTauTauTree->Branch("Mu3Energy", &Mu3Energy, "Mu3Energy/D");
    MuMuTauTauTree->Branch("Mu1Pt", &Mu1Pt, "Mu1Pt/D");
    MuMuTauTauTree->Branch("Mu2Pt", &Mu2Pt, "Mu2Pt/D");
    MuMuTauTauTree->Branch("Mu3Pt", &Mu3Pt, "Mu3Pt/D");
    MuMuTauTauTree->Branch("Mu1Eta", &Mu1Eta, "Mu1Eta/D");
    MuMuTauTauTree->Branch("Mu2Eta", &Mu2Eta, "Mu2Eta/D");
    MuMuTauTauTree->Branch("Mu3Eta", &Mu3Eta, "Mu3Eta/D");
    MuMuTauTauTree->Branch("Mu1Phi", &Mu1Phi, "Mu1Phi/D");
    MuMuTauTauTree->Branch("Mu2Phi", &Mu2Phi, "Mu2Phi/D");
    MuMuTauTauTree->Branch("Mu3Phi", &Mu3Phi, "Mu3Phi/D");
    
    MuMuTauTauTree->Branch("Mu1Iso", &Mu1Iso, "Mu1Iso/D");
    MuMuTauTauTree->Branch("Mu2Iso", &Mu2Iso, "Mu2Iso/D");
    MuMuTauTauTree->Branch("Mu3Iso", &Mu3Iso, "Mu3Iso/D");

    MuMuTauTauTree->Branch("tauHadEnergy", &tauHadEnergy, "tauHadEnergy/D");
    MuMuTauTauTree->Branch("tauHadPt", &tauHadPt, "tauHadPt/D");
    MuMuTauTauTree->Branch("tauHadEta", &tauHadEta, "tauHadEta/D");
    MuMuTauTauTree->Branch("tauHadPhi", &tauHadPhi, "tauHadPhi/D");
    MuMuTauTauTree->Branch("tauIsolationMVArawValue", &tauIsolationMVArawValue, "tauIsolationMVArawValue/D");
    MuMuTauTauTree->Branch("tauIsolationMVAVVLoose", &tauIsolationMVAVVLoose, "tauIsolationMVAVVLoose/D");
    MuMuTauTauTree->Branch("tauIsolationMVAVLoose", &tauIsolationMVAVLoose, "tauIsolationMVAVLoose/D");
    MuMuTauTauTree->Branch("tauIsolationMVALoose", &tauIsolationMVALoose, "tauIsolationMVALoose/D");
    MuMuTauTauTree->Branch("tauIsolationMVAMedium", &tauIsolationMVAMedium, "tauIsolationMVAMedium/D");
    MuMuTauTauTree->Branch("tauIsolationMVATight", &tauIsolationMVATight, "tauIsolationMVATight/D");
    MuMuTauTauTree->Branch("tauIsolationMVAVTight", &tauIsolationMVAVTight, "tauIsolationMVAVTight/D");
    MuMuTauTauTree->Branch("tauIsolationMVAVVTight", &tauIsolationMVAVVTight, "tauIsolationMVAVVTight/D");

    MuMuTauTauTree->Branch("EventWeight", &EventWeight, "EventWeight/F");
    MuMuTauTauTree->Branch("NPVertex", &NPVertex, "NPVertex/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTauTauAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTauTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTauTauAnalyzer);
