// -*- C++ -*-
//
// Package:    MuMuChannel/MuonID
// Class:      MuonID
// 
/**\class MuonID MuonID.cc MuMuChannel/MuonID/plugins/MuonID.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 09 Apr 2019 16:34:02 GMT
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

class MuonID : public edm::stream::EDFilter<> {
   public:
      explicit MuonID(const edm::ParameterSet&);
      ~MuonID();

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
      std::string muonID_;
      int minNumObjsToPassFilter_;
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
MuonID::MuonID(const edm::ParameterSet& iConfig):
    muonTag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   muonID_ = iConfig.getParameter<std::string>("muonID");
   minNumObjsToPassFilter_ = iConfig.getParameter<int>("minNumObjsToPassFilter");
   transform(muonID_.begin(), muonID_.end(), muonID_.begin(), ::toupper);
}


MuonID::~MuonID()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuonID::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   int CountMuon=0;

   Handle<edm::View<pat::Muon>> pMuons;
   iEvent.getByToken(muonTag_, pMuons);

   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();

   if (pMuons->size() < 1) return 0;

   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end(); ++iMuon)
   {
       bool goodGlob = iMuon->isGlobalMuon() && iMuon->globalTrack()->normalizedChi2() < 3 && iMuon->combinedQuality().chi2LocalPosition < 12 && iMuon->combinedQuality().trkKink < 20;
       bool isMedium = muon::isLooseMuon(*iMuon) && iMuon->innerTrack()->validFraction() > 0.8 && muon::segmentCompatibility(*iMuon) > (goodGlob ? 0.303 : 0.451);
       bool isLoose = iMuon->isPFMuon() && (iMuon->isGlobalMuon() || iMuon->isTrackerMuon());
       bool isTight = iMuon->isGlobalMuon() && iMuon->isPFMuon() && iMuon->globalTrack()->normalizedChi2() < 10 && iMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && iMuon->numberOfMatchedStations() > 1 && fabs(iMuon->muonBestTrack()->dxy()) < 0.2 && fabs(iMuon->muonBestTrack()->dz()) < 0.5 && iMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
       if ((isMedium && muonID_=="MEDIUM") || (isLoose && muonID_=="LOOSE") || (isTight && muonID_=="TIGHT"))
       {
           CountMuon+=1;
           muonColl->push_back(*iMuon);
       }
   }

   if (CountMuon>=1)
   {
       iEvent.put(std::move(muonColl));
       return (CountMuon >= minNumObjsToPassFilter_);
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MuonID::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuonID::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuonID::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MuonID::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonID::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonID::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonID::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonID);
