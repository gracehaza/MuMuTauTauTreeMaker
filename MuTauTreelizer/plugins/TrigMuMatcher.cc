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
      double muPtCut_;
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
   muPtCut_ = iConfig.getParameter<double>("muPtCut");

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
   edm::Handle<edm::View<pat::Muon>> pMuons;
   edm::Handle<edm::TriggerResults> pTriggerBits;
   edm::Handle<std::vector<pat::TriggerObjectStandAlone> > pTriggerObjects;
   
   iEvent.getByToken(muonsTag_, pMuons);
   iEvent.getByToken(triggerBits_, pTriggerBits);
   iEvent.getByToken(triggerObjects_, pTriggerObjects);

   const edm::TriggerNames &names = iEvent.triggerNames(*pTriggerBits);
   std::vector<unsigned int> corrTrigSpot;

   // --- find the user customized trigger path names in the existing trigger path list in TriggerResults of HLT ---
   for (unsigned int i = 0, n = pTriggerBits->size(); i<n; ++i)
   {
       for ( std::string iName : trigNames_ ){
           if (names.triggerName(i).find(iName) != std::string::npos )
           {
               corrTrigSpot.push_back(i);
           }
       }
   }

   // --- prepare for the reco-muon object for trigger matching --- 
   pat::Muon recoLeadingMu;
   pat::TriggerObjectStandAlone trigObj;
   bool checkPassEvent = false;

   for (uint iobj = 0; iobj < pTriggerObjects->size(); iobj++)
   {
       trigObj = pTriggerObjects->at(iobj);
       trigObj.unpackPathNames(names); // get the list of trigger paths of each triggered object
       bool checkObjMatch = false; // This variable will be used to check if this object(muon) passes at least one customized trigger path

       for (unsigned int num : corrTrigSpot) // Loop over the trigger paths of TriggerResults cluster that are compatible with customized trigger paths
       {
           const std::string& name = names.triggerName(num);
           if (trigObj.hasPathName(name, true) && !checkObjMatch) // if there is no yet reco-muon candidates matched with the triggered object
           {
               for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end(); ++iMuon) // loop over all the reco-muons
               {
                   double dRCurr = deltaR(*iMuon, trigObj); // use dR to match the reco-muon and the triggered object
                   if (iMuon->pt() > muPtCut_ && dRCurr < dRCut_ && fabs(iMuon->pdgId()) == fabs(trigObj.pdgId())) // match reco-mu with trigger-mu
                   {
                       recoLeadingMu = *iMuon;
                       checkPassEvent = true;
                       checkObjMatch = true;
                       break;
                   } // end if reco-mu matched with trigger-mu
               } // end loop of reconstructed muons 
           } // end if there is at least one reco-muon candidate matched with the triggered object
       } // end loop of trigger paths in TriggerResults cluster that are compatible with customized trigger paths
   } // end loop of all the triggered objects

   // --- resort reco-mu (the one matched with trigger-mu will be the leading muon) --- 
   if (checkPassEvent)
   {
       muonColl->push_back(recoLeadingMu);

       for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end(); ++iMuon)
       {
           if (deltaR(*iMuon, recoLeadingMu) < 0.0001 && (iMuon->pt() - recoLeadingMu.pt()) < 0.00001) continue; 
           muonColl->push_back(*iMuon);
       } // end for loop on reco-muons

       iEvent.put(std::move(muonColl));
   } // end if checkPassEvent == true

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
