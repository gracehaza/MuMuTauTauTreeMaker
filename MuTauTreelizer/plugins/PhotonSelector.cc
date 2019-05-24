// -*- C++ -*-
//
// Package:    testED/PhotonSelector
// Class:      PhotonSelector
// 
/**\class PhotonSelector PhotonSelector.cc testED/PhotonSelector/plugins/PhotonSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Fri, 24 May 2019 15:39:05 GMT
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

#include "DataFormats/PatCandidates/interface/Photon.h"
#include <string>
#include <math.h>
//
// class declaration
//

class PhotonSelector : public edm::stream::EDFilter<> {
   public:
      explicit PhotonSelector(const edm::ParameterSet&);
      ~PhotonSelector();

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
      edm::EDGetTokenT<edm::View<pat::Photon>> photonTag_;
      std::string relIdName_;
      bool passRelId_;
      double etaCut_;
      double ptCut_;
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
PhotonSelector::PhotonSelector(const edm::ParameterSet& iConfig):
    photonTag_(consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photonTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Photon>>();
   relIdName_ = iConfig.getParameter<std::string>("relIdName");
   passRelId_ = iConfig.getParameter<bool>("passRelId");
   etaCut_ = iConfig.getParameter<double>("etaCut");
   ptCut_ = iConfig.getParameter<double>("ptCut");
}


PhotonSelector::~PhotonSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PhotonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Photon>> photonColl = std::make_unique<std::vector<pat::Photon>>();
   edm::Handle<edm::View<pat::Photon>> pPhotons;
   iEvent.getByToken(photonTag_, pPhotons);

   if (pPhotons->size() < 1) 
   {
       iEvent.put(std::move(photonColl));
       return true;
   }

   int CountPhoton = 0;
   for(edm::View<pat::Photon>::const_iterator iPhoton=pPhotons->begin(); iPhoton!=pPhotons->end(); ++iPhoton)
   {
       if ((iPhoton->photonID(relIdName_) == passRelId_) && (iPhoton->pt() > ptCut_) && (fabs(iPhoton->eta()) < etaCut_))
       {
           CountPhoton++;
           photonColl->push_back(*iPhoton);
       }
   }

   iEvent.put(std::move(photonColl));
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PhotonSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PhotonSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PhotonSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PhotonSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PhotonSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PhotonSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonSelector);
