// -*- C++ -*-
//
// Package:    testED/TrigMuMatcher
// Class:      TrigMuMatcher
// 
/**\class TrigMuMatcher TrigMuMatcher.cc testED/TrigMuMatcher/plugins/TrigMuMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Mon, 15 Apr 2019 13:37:25 GMT
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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"
//
// class declaration
//

class TrigMuMatcher : public edm::stream::EDFilter<> {
   public:
      explicit TrigMuMatcher(const edm::ParameterSet&);
      ~TrigMuMatcher();

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
      edm::EDGetTokenT<edm::View<pat::Muon>> muonsTag_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
      std::vector<std::string> trigNames_;
      double dRCut_;
      double mu1PtCut_;
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
TrigMuMatcher::TrigMuMatcher(const edm::ParameterSet& iConfig):
    muonsTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("triggerObjects")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   trigNames_ = iConfig.getParameter<std::vector<std::string>>("trigNames");
   dRCut_ = iConfig.getParameter<double>("dRCut");
   mu1PtCut_ = iConfig.getParameter<double>("mu1PtCut");

}


TrigMuMatcher::~TrigMuMatcher()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrigMuMatcher::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   std::unique_ptr<std::vector<pat::Muon>> leadMuonColl = std::make_unique<std::vector<pat::Muon>>();
   edm::Handle<edm::View<pat::Muon>> pMuons;
   edm::Handle<edm::TriggerResults> pTriggerBits;
   edm::Handle<std::vector<pat::TriggerObjectStandAlone> > pTriggerObjects;
   
   iEvent.getByToken(muonsTag_, pMuons);
   iEvent.getByToken(triggerBits_, pTriggerBits);
   iEvent.getByToken(triggerObjects_, pTriggerObjects);

   const edm::TriggerNames &names = iEvent.triggerNames(*pTriggerBits);
   std::vector<unsigned int> corrTrigSpot;

   for (unsigned int i = 0, n = pTriggerBits->size(); i<n; ++i)
   {
       for ( std::string iName : trigNames_ ){
           if (names.triggerName(i).find(iName) != std::string::npos )
           {
               corrTrigSpot.push_back(i);
           }
       }
   }

   pat::TriggerObjectStandAlone TO;
   bool checkPassEvent = false;
   for ( uint iobj = 0; iobj <pTriggerObjects->size(); iobj++ )
   {
       TO = pTriggerObjects->at(iobj);
       TO.unpackPathNames(names); // get the list of trigger paths in the input file
       bool checkObjMatch = false; // This variable will be used to check if this object(muon) passes at least one customized trigger path

       for (unsigned int num : corrTrigSpot) // Loop through the customized trigger paths
       {
           const std::string& name = names.triggerName(num);
           if (TO.hasPathName(name, true) && !checkObjMatch) // if there is no yet matched objects matched with the reconstructed muon candidate
           {
               for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon) // loop through the muons
               {
                   double dRCurr = deltaR(*iMuon, TO);
                   if (iMuon->pt() > mu1PtCut_ && dRCurr < dRCut_) // select the muon closest to the reconstructed muon and having highest pt
                   {
                       checkPassEvent = true;
                       checkObjMatch = true;
                       muonColl->push_back(*iMuon);
                       break;
                   } // end if pt/dR requirement
               } // end loop of reconstructed muons 
           } // end if there is matched objects with the reco-muon candidate
       } // end loop of customized trigger muon paths
   } // end loop of pat trigger objects

   // we are going to find the leading muon 
   if (muonColl->size() > 0){
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

      leadMuonColl->push_back(highestPtMuon);
      iEvent.put(std::move(leadMuonColl));
   }

   return (checkPassEvent);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TrigMuMatcher::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TrigMuMatcher::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TrigMuMatcher::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TrigMuMatcher::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TrigMuMatcher::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TrigMuMatcher::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrigMuMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrigMuMatcher);
