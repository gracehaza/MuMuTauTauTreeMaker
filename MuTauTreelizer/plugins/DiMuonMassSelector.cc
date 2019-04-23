// -*- C++ -*-
//
// Package:    testED/DiMuonMassSelector
// Class:      DiMuonMassSelector
// 
/**\class DiMuonMassSelector DiMuonMassSelector.cc testED/DiMuonMassSelector/plugins/DiMuonMassSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 16 Apr 2019 15:53:33 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
//
// class declaration
//

class DiMuonMassSelector : public edm::stream::EDFilter<> {
   public:
      explicit DiMuonMassSelector(const edm::ParameterSet&);
      ~DiMuonMassSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> mu1Tag_;
      edm::EDGetTokenT<edm::View<pat::Muon>> mu2Tag_;
      double minMass_;
      double maxMass_;
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
DiMuonMassSelector::DiMuonMassSelector(const edm::ParameterSet& iConfig) :
    mu1Tag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("mu1Tag"))),
    mu2Tag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("mu2Tag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   minMass_ = iConfig.getParameter<double>("minMass");
   maxMass_ = iConfig.getParameter<double>("maxMass");

}


DiMuonMassSelector::~DiMuonMassSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiMuonMassSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   edm::Handle<edm::View<pat::Muon>> pMu1;
   edm::Handle<edm::View<pat::Muon>> pMu2;

   iEvent.getByToken(mu1Tag_, pMu1);
   iEvent.getByToken(mu2Tag_, pMu2);

   pat::Muon mu1 = pMu1->at(0);
   pat::Muon mu2 = pMu2->at(0);

   double invMass = (mu1.p4() + mu2.p4()).M();
   if (invMass >= minMass_ && invMass <= maxMass_)
   {
       muonColl->push_back(mu1);
       muonColl->push_back(mu2);
       iEvent.put(std::move(muonColl));
       return true;
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DiMuonMassSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DiMuonMassSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DiMuonMassSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DiMuonMassSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DiMuonMassSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DiMuonMassSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuonMassSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonMassSelector);
