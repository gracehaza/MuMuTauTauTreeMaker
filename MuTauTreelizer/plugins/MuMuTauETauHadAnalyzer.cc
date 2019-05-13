// -*- C++ -*-
//
// Package:    testED/MuMuTauETauHadAnalyzer
// Class:      MuMuTauETauHadAnalyzer
// 
/**\class MuMuTauETauHadAnalyzer MuMuTauETauHadAnalyzer.cc testED/MuMuTauETauHadAnalyzer/plugins/MuMuTauETauHadAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Mon, 13 May 2019 07:55:33 GMT
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

class MuMuTauETauHadAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTauETauHadAnalyzer(const edm::ParameterSet&);
      ~MuMuTauETauHadAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu1Mu2_;
      edm::EDGetTokenT<edm::View<pat::Electron>> Ele_;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauHad_;
      std::vector<std::string> tauDiscriminatorTags_;
      edm::EDGetTokenT<edm::View<reco::Vertex>> Vertex_;
      bool isMC_;
      float EventWeight;
      float NPVertex;
      edm::EDGetTokenT<GenEventInfoProduct> generator_;

      TTree *MuMuTauETauHadTree;
      double Mu1Energy;
      double Mu2Energy;
      double EleEnergy;
      double tauHadEnergy;
      double Mu1Pt;
      double Mu2Pt;
      double ElePt;
      double tauHadPt;
      double Mu1Eta;
      double Mu2Eta;
      double EleEta;
      double tauHadEta;
      double Mu1Phi;
      double Mu2Phi;
      double ElePhi;
      double tauHadPhi;
      double Mu1Iso;
      double Mu2Iso;
      // ---- below is variables for electron energy scale and smearing application ---
      float ecalTrkEnergyPreCorr;
      float ecalTrkEnergyErrPreCorr;
      float ecalTrkEnergyPostCorr;
      float ecalTrkEnergyErrPostCorr;
      // ---- below is variables for Tau discriminators ---
      float tauIsolationMVArawValue;
      float tauIsolationMVAVVLoose;
      float tauIsolationMVAVLoose;
      float tauIsolationMVALoose;
      float tauIsolationMVAMedium;
      float tauIsolationMVATight;
      float tauIsolationMVAVTight;
      float tauIsolationMVAVVTight;
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
MuMuTauETauHadAnalyzer::MuMuTauETauHadAnalyzer(const edm::ParameterSet& iConfig):
    Mu1Mu2_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2Tag"))),
    Ele_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("EleTag"))),
    TauHad_(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauHadTag"))),
    Vertex_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("Vertex"))),
    generator_(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag()))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC_ = iConfig.getParameter<bool>("isMC");
   ecalTrkEnergyPreCorr = -999.0;
   ecalTrkEnergyErrPreCorr = -999.0;
   ecalTrkEnergyPostCorr = -999.0;
   ecalTrkEnergyErrPostCorr = -999.0;
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


MuMuTauETauHadAnalyzer::~MuMuTauETauHadAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTauETauHadAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);

   edm::Handle<edm::View<pat::Electron>> pElectron;
   iEvent.getByToken(Ele_, pElectron);

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

   pat::Muon Mu1 = pMu1Mu2->at(0);
   pat::Muon Mu2 = pMu1Mu2->at(1);
   pat::Electron Ele = pElectron->at(0);

   Mu1Energy = Mu1.energy();
   Mu1Pt = Mu1.pt();
   Mu1Eta = Mu1.eta();
   Mu1Phi = Mu1.phi();

   Mu2Energy = Mu2.energy();
   Mu2Pt = Mu2.pt();
   Mu2Eta = Mu2.eta();
   Mu2Phi = Mu2.phi();

   EleEnergy = Ele.energy();
   ElePt = Ele.pt();
   EleEta = Ele.eta();
   ElePhi = Ele.phi();

   reco::MuonPFIsolation iso = Mu1.pfIsolationR04();
   Mu1Iso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / Mu1.pt();

   iso = Mu2.pfIsolationR04();
   Mu2Iso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / Mu2.pt();

   // ---- electron energy scale correction & resolution smear ----
   // refer to: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
   ecalTrkEnergyPreCorr = Ele.userFloat("ecalTrkEnergyPreCorr");
   ecalTrkEnergyErrPreCorr = Ele.userFloat("ecalTrkEnergyErrPreCorr");
   ecalTrkEnergyPostCorr = Ele.userFloat("ecalTrkEnergyPostCorr");
   ecalTrkEnergyErrPostCorr = Ele.userFloat("ecalTrkEnergyErrPostCorr");

   pat::Tau tauHad = pTau->at(0);
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

   MuMuTauETauHadTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTauETauHadAnalyzer::beginJob()
{
    // -- initialize the skimmed flat tree --
    edm::Service<TFileService> fileService;
    MuMuTauETauHadTree = fileService->make<TTree>("MuMuTauETauHadTree","MuMuTauETauHadTree");
    
    MuMuTauETauHadTree->Branch("Mu1Energy", &Mu1Energy, "Mu1Energy/D");
    MuMuTauETauHadTree->Branch("Mu2Energy", &Mu2Energy, "Mu2Energy/D");
    MuMuTauETauHadTree->Branch("EleEnergy", &EleEnergy, "EleEnergy/D");
    MuMuTauETauHadTree->Branch("Mu1Pt", &Mu1Pt, "Mu1Pt/D");
    MuMuTauETauHadTree->Branch("Mu2Pt", &Mu2Pt, "Mu2Pt/D");
    MuMuTauETauHadTree->Branch("ElePt", &ElePt, "ElePt/D");
    MuMuTauETauHadTree->Branch("Mu1Eta", &Mu1Eta, "Mu1Eta/D");
    MuMuTauETauHadTree->Branch("Mu2Eta", &Mu2Eta, "Mu2Eta/D");
    MuMuTauETauHadTree->Branch("EleEta", &EleEta, "EleEta/D");
    MuMuTauETauHadTree->Branch("Mu1Phi", &Mu1Phi, "Mu1Phi/D");
    MuMuTauETauHadTree->Branch("Mu2Phi", &Mu2Phi, "Mu2Phi/D");
    MuMuTauETauHadTree->Branch("ElePhi", &ElePhi, "ElePhi/D");
    
    MuMuTauETauHadTree->Branch("Mu1Iso", &Mu1Iso, "Mu1Iso/D");
    MuMuTauETauHadTree->Branch("Mu2Iso", &Mu2Iso, "Mu2Iso/D");

    MuMuTauETauHadTree->Branch("ecalTrkEnergyPreCorr", &ecalTrkEnergyPreCorr, "ecalTrkEnergyPreCorr/F");
    MuMuTauETauHadTree->Branch("ecalTrkEnergyErrPreCorr", &ecalTrkEnergyErrPreCorr, "ecalTrkEnergyErrPreCorr/F");
    MuMuTauETauHadTree->Branch("ecalTrkEnergyPostCorr", &ecalTrkEnergyPostCorr, "ecalTrkEnergyPostCorr/F");
    MuMuTauETauHadTree->Branch("ecalTrkEnergyErrPostCorr", &ecalTrkEnergyErrPostCorr, "ecalTrkEnergyErrPostCorr/F");

    MuMuTauETauHadTree->Branch("tauHadEnergy", &tauHadEnergy, "tauHadEnergy/D");
    MuMuTauETauHadTree->Branch("tauHadPt", &tauHadPt, "tauHadPt/D");
    MuMuTauETauHadTree->Branch("tauHadEta", &tauHadEta, "tauHadEta/D");
    MuMuTauETauHadTree->Branch("tauHadPhi", &tauHadPhi, "tauHadPhi/D");

    MuMuTauETauHadTree->Branch("tauIsolationMVArawValue", &tauIsolationMVArawValue, "tauIsolationMVArawValue/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVAVVLoose", &tauIsolationMVAVVLoose, "tauIsolationMVAVVLoose/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVAVLoose", &tauIsolationMVAVLoose, "tauIsolationMVAVLoose/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVALoose", &tauIsolationMVALoose, "tauIsolationMVALoose/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVAMedium", &tauIsolationMVAMedium, "tauIsolationMVAMedium/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVATight", &tauIsolationMVATight, "tauIsolationMVATight/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVAVTight", &tauIsolationMVAVTight, "tauIsolationMVAVTight/F");
    MuMuTauETauHadTree->Branch("tauIsolationMVAVVTight", &tauIsolationMVAVVTight, "tauIsolationMVAVVTight/F");

    MuMuTauETauHadTree->Branch("EventWeight", &EventWeight, "EventWeight/F");
    MuMuTauETauHadTree->Branch("NPVertex", &NPVertex, "NPVertex/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTauETauHadAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTauETauHadAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTauETauHadAnalyzer);
