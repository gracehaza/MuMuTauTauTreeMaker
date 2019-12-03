// -*- C++ -*-
//
// Package:    testED/GenMuonCandSelector
// Class:      GenMuonCandSelector
// 
/**\class GenMuonCandSelector GenMuonCandSelector.cc testED/GenMuonCandSelector/plugins/GenMuonCandSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Mon, 13 May 2019 16:23:05 GMT
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <math.h>

using namespace std;
//
// class declaration
//

class GenMuonCandSelector : public edm::stream::EDFilter<> {
   public:
      explicit GenMuonCandSelector(const edm::ParameterSet&);
      ~GenMuonCandSelector();

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
      edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleTag_;
      double ptCut_;
      double etaCut_;
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
GenMuonCandSelector::GenMuonCandSelector(const edm::ParameterSet& iConfig):
    genParticleTag_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<reco::GenParticle>>();
   ptCut_ = iConfig.getParameter<double>("ptCut");
   etaCut_ = iConfig.getParameter<double>("etaCut");
}


GenMuonCandSelector::~GenMuonCandSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenMuonCandSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<reco::GenParticle>> muonColl = std::make_unique<std::vector<reco::GenParticle>>();
   edm::Handle<edm::View<reco::GenParticle>> pGenParticles;
   iEvent.getByToken(genParticleTag_, pGenParticles);

   if(pGenParticles->size() < 1)
   {
       iEvent.put(std::move(muonColl));
       return true;
   }

   int CountGenMuon = 0;

   for(edm::View<reco::GenParticle>::const_iterator iParticle=pGenParticles->begin(); iParticle!=pGenParticles->end(); ++iParticle)
   {
       if(fabs(iParticle->pdgId()) == 13 && fabs(iParticle->mother()->pdgId()) != 13 && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
       {
           CountGenMuon++;
           muonColl->push_back(*iParticle);
       }
   }

   iEvent.put(std::move(muonColl));
   return true;

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenMuonCandSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenMuonCandSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenMuonCandSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenMuonCandSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenMuonCandSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenMuonCandSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenMuonCandSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenMuonCandSelector);
