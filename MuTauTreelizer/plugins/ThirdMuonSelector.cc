// -*- C++ -*-
//
// Package:    testED/ThirdMuonSelector
// Class:      ThirdMuonSelector
// 
/**\class ThirdMuonSelector ThirdMuonSelector.cc testED/ThirdMuonSelector/plugins/ThirdMuonSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Fri, 26 Apr 2019 13:15:29 GMT
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

class ThirdMuonSelector : public edm::stream::EDFilter<> {
   public:
      explicit ThirdMuonSelector(const edm::ParameterSet&);
      ~ThirdMuonSelector();

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
ThirdMuonSelector::ThirdMuonSelector(const edm::ParameterSet& iConfig):
    muonTag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag"))),
    mu1mu2Tag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("mu1mu2Tag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   dRCut_ = iConfig.getParameter<double>("dRCut");
}


ThirdMuonSelector::~ThirdMuonSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ThirdMuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   std::unique_ptr<std::vector<pat::Muon>> thirdMuonColl = std::make_unique<std::vector<pat::Muon>>();
   edm::Handle<edm::View<pat::Muon>> pMuons;
   iEvent.getByToken(muonTag_, pMuons);
   
   if (pMuons->size() < 1) return 0;
   edm::Handle<edm::View<pat::Muon>> pMu1Mu2;
   iEvent.getByToken(mu1mu2Tag_, pMu1Mu2);

   pat::Muon mu1 = pMu1Mu2->at(0);
   pat::Muon mu2 = pMu1Mu2->at(1);

   int CountMuon=0;
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

   // below we are going to select the muon with third highest pt
   if (CountMuon >= 1)
   {
       double highestPt = -1;
       pat::Muon highestPtMuon;
       for (std::vector<pat::Muon>::iterator iMuon=muonColl->begin(); iMuon!=muonColl->end(); ++iMuon)
       {
           if (iMuon->pt() > highestPt)
           {
               highestPt = iMuon->pt();
               highestPtMuon = *iMuon;
           }
       }

       thirdMuonColl->push_back(highestPtMuon);
       iEvent.put(std::move(thirdMuonColl));
       return true;
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ThirdMuonSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ThirdMuonSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ThirdMuonSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ThirdMuonSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ThirdMuonSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ThirdMuonSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThirdMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ThirdMuonSelector);
