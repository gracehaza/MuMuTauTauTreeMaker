// -*- C++ -*-
//
// Package:    testED/LeadingMuonIso
// Class:      LeadingMuonIso
// 
/**\class LeadingMuonIso LeadingMuonIso.cc testED/LeadingMuonIso/plugins/LeadingMuonIso.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Thu, 11 Apr 2019 17:17:27 GMT
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

class LeadingMuonIso : public edm::stream::EDFilter<> {
   public:
      explicit LeadingMuonIso(const edm::ParameterSet&);
      ~LeadingMuonIso();

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
      edm::EDGetTokenT<edm::View<pat::Muon> > muonTag_;
      double relIsoCutVal_;
      bool passRelIso_;
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
LeadingMuonIso::LeadingMuonIso(const edm::ParameterSet& iConfig):
    muonTag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   relIsoCutVal_ = iConfig.getParameter<double>("relIsoCutVal");
   passRelIso_ = iConfig.getParameter<bool>("passRelIso");

}


LeadingMuonIso::~LeadingMuonIso()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
LeadingMuonIso::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   edm::Handle<edm::View<pat::Muon>> pMuons;
   iEvent.getByToken(muonTag_, pMuons);

   if (pMuons->size() < 1) return 0;

   int CountMuon=0;
   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end(); ++iMuon)
   {
       reco::MuonPFIsolation iso = iMuon->pfIsolationR04();
       double reliso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / iMuon->pt();
       if ((reliso < relIsoCutVal_ && passRelIso_) || (reliso > relIsoCutVal_ && !passRelIso_) || relIsoCutVal_ == -1)
       {
           CountMuon+=1;
           muonColl->push_back(*iMuon);
       }
   }

   if (CountMuon >= 1)
   {
       iEvent.put(std::move(muonColl));
       return true;
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
LeadingMuonIso::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
LeadingMuonIso::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
LeadingMuonIso::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LeadingMuonIso::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LeadingMuonIso::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LeadingMuonIso::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeadingMuonIso::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(LeadingMuonIso);
