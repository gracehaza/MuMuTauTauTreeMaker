// -*- C++ -*-
//
// Package:    MuMuChannel/DiMuonAnalyzer
// Class:      DiMuonAnalyzer
// 
/**\class DiMuonAnalyzer DiMuonAnalyzer.cc MuMuChannel/DiMuonAnalyzer/plugins/DiMuonAnalyzer.cc

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TTree.h"
#include <math.h> 
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DiMuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DiMuonAnalyzer(const edm::ParameterSet&);
      ~DiMuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu1Mu2_;
      bool isMC_;
      // FIXME: we need pileup information for MC reweighting
      float EventWeight;
      edm::EDGetTokenT<GenEventInfoProduct> generator_;

      TTree *Mu1Mu2Tree;
      double Mu1Energy;
      double Mu2Energy;
      double Mu1Pt;
      double Mu2Pt;
      double Mu1Eta;
      double Mu2Eta;
      double Mu1Phi;
      double Mu2Phi;
      double invMassMu1Mu2;
      double PtMu1Mu2;
      double EtaMu1Mu2;
      double PhiMu1Mu2;
      double dRMu1Mu2;
      double dEtaMu1Mu2;
      double dPhiMu1Mu2;
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
DiMuonAnalyzer::DiMuonAnalyzer(const edm::ParameterSet& iConfig):
    Mu1Mu2_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
    generator_(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag()))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC_ = iConfig.getParameter<bool>("isMC");
}


DiMuonAnalyzer::~DiMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);

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

   pat::Muon Mu1 = pMu1Mu2->at(0);
   pat::Muon Mu2 = pMu1Mu2->at(1);

   Mu1Energy = Mu1.energy();
   Mu1Pt = Mu1.pt();
   Mu1Eta = Mu1.eta();
   Mu1Phi = Mu1.phi();

   Mu2Energy = Mu2.energy();
   Mu2Pt = Mu2.pt();
   Mu2Eta = Mu2.eta();
   Mu2Phi = Mu2.phi();

   reco::Candidate::LorentzVector p4DiMu = Mu1.p4() + Mu2.p4();
   invMassMu1Mu2 = p4DiMu.mass();
   PtMu1Mu2 = p4DiMu.pt();
   EtaMu1Mu2 = p4DiMu.eta();
   PhiMu1Mu2 = p4DiMu.phi();

   dRMu1Mu2 = deltaR(Mu1,Mu2);
   dEtaMu1Mu2 = fabs(Mu1.eta() - Mu2.eta());
   dPhiMu1Mu2 = deltaPhi(Mu1.phi(),Mu2.phi());

   Mu1Mu2Tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
DiMuonAnalyzer::beginJob()
{
    // -- initialize the skimmed flat tree --
    edm::Service<TFileService> fileService;
    Mu1Mu2Tree = fileService->make<TTree>("Mu1Mu2Tree","Mu1Mu2Tree");
    
    Mu1Mu2Tree->Branch("Mu1Energy", &Mu1Energy, "Mu1Energy/D");
    Mu1Mu2Tree->Branch("Mu2Energy", &Mu2Energy, "Mu2Energy/D");
    Mu1Mu2Tree->Branch("Mu1Pt", &Mu1Pt, "Mu1Pt/D");
    Mu1Mu2Tree->Branch("Mu2Pt", &Mu2Pt, "Mu2Pt/D");
    Mu1Mu2Tree->Branch("Mu1Eta", &Mu1Eta, "Mu1Eta/D");
    Mu1Mu2Tree->Branch("Mu2Eta", &Mu2Eta, "Mu2Eta/D");
    Mu1Mu2Tree->Branch("Mu1Phi", &Mu1Phi, "Mu1Phi/D");
    Mu1Mu2Tree->Branch("Mu2Phi", &Mu2Phi, "Mu2Phi/D");
    
    Mu1Mu2Tree->Branch("invMassMu1Mu2", &invMassMu1Mu2, "invMassMu1Mu2/D");
    Mu1Mu2Tree->Branch("PtMu1Mu2", &PtMu1Mu2, "PtMu1Mu2/D");
    Mu1Mu2Tree->Branch("EtaMu1Mu2", &EtaMu1Mu2, "EtaMu1Mu2/D");
    Mu1Mu2Tree->Branch("PhiMu1Mu2", &PhiMu1Mu2, "PhiMu1Mu2/D");
    Mu1Mu2Tree->Branch("dRMu1Mu2", &dRMu1Mu2, "dRMu1Mu2/D");
    Mu1Mu2Tree->Branch("dEtaMu1Mu2", &dEtaMu1Mu2, "dEtaMu1Mu2/D");
    Mu1Mu2Tree->Branch("dPhiMu1Mu2", &dPhiMu1Mu2, "dPhiMu1Mu2/D");

    Mu1Mu2Tree->Branch("EventWeight", &EventWeight, "EventWeight/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiMuonAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonAnalyzer);
