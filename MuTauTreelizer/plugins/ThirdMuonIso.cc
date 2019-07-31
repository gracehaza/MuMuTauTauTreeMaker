// -*- C++ -*-
//
// Package:    testED/ThirdMuonIso
// Class:      ThirdMuonIso
// 
/**\class ThirdMuonIso ThirdMuonIso.cc testED/ThirdMuonIso/plugins/ThirdMuonIso.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Wed, 31 Jul 2019 11:19:48 GMT
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
#include <math.h>
//
// class declaration
//

class ThirdMuonIso : public edm::stream::EDFilter<> {
   public:
      explicit ThirdMuonIso(const edm::ParameterSet&);
      ~ThirdMuonIso();

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
      edm::EDGetTokenT<edm::View<pat::Muon>> muonTag_;
      edm::EDGetTokenT<edm::View<pat::Muon>> mu1mu2Tag_;
      double dRCut_;
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
ThirdMuonIso::ThirdMuonIso(const edm::ParameterSet& iConfig):
    muonTag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag"))),
    mu1mu2Tag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("mu1mu2Tag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   dRCut_ = iConfig.getParameter<double>("dRCut");
}


ThirdMuonIso::~ThirdMuonIso()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ThirdMuonIso::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   edm::Handle<edm::View<pat::Muon>> pMuons;
   iEvent.getByToken(muonTag_, pMuons);
   
   if (pMuons->size() < 1) 
   {
       iEvent.put(std::move(muonColl));
       return true; // we also select the events containing no more than two muons!
   }

   edm::Handle<edm::View<pat::Muon>> pMu1Mu2;
   iEvent.getByToken(mu1mu2Tag_, pMu1Mu2);

   pat::Muon mu1 = pMu1Mu2->at(0);
   pat::Muon mu2 = pMu1Mu2->at(1);

   int CountMuon=0;
   // below we are going to select all the other muon candidates which could be the third muon

   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end(); ++iMuon)
   {
       if ((deltaR(*iMuon, mu1) < 0.0001) && (fabs(iMuon->pt()-mu1.pt()) < 0.0001) && (dRCut_ <= 0))
           continue;
       else if((deltaR(*iMuon, mu2) < 0.0001) && (fabs(iMuon->pt()-mu2.pt()) < 0.0001) && (dRCut_ <= 0)) 
           continue;
       else if((deltaR(*iMuon, mu1) < dRCut_) && (dRCut_ > 0)) 
           continue;
       else if((deltaR(*iMuon, mu2) < dRCut_) && (dRCut_ > 0)) 
           continue;
       else{
           CountMuon++;
           muonColl->push_back(*iMuon);
       }
   }

   iEvent.put(std::move(muonColl));
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ThirdMuonIso::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ThirdMuonIso::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ThirdMuonIso::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ThirdMuonIso::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ThirdMuonIso::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ThirdMuonIso::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThirdMuonIso::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ThirdMuonIso);
