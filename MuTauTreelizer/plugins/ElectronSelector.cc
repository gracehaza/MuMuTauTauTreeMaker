// -*- C++ -*-
//
// Package:    testED/ElectronSelector
// Class:      ElectronSelector
// 
/**\class ElectronSelector ElectronSelector.cc testED/ElectronSelector/plugins/ElectronSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 07 May 2019 14:58:43 GMT
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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include <string>
//
// class declaration
//

class ElectronSelector : public edm::stream::EDFilter<> {
   public:
      explicit ElectronSelector(const edm::ParameterSet&);
      ~ElectronSelector();

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
      edm::EDGetTokenT<edm::View<pat::Electron>> electronTag_;
      std::string relIdName_;
      bool passRelId_;
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
ElectronSelector::ElectronSelector(const edm::ParameterSet& iConfig):
    electronTag_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Electron>>();
   relIdName_ = iConfig.getParameter<std::string>("relIdName");
   passRelId_ = iConfig.getParameter<bool>("passRelId");
}


ElectronSelector::~ElectronSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ElectronSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Electron>> electronColl = std::make_unique<std::vector<pat::Electron>>();
   edm::Handle<edm::View<pat::Electron>> pElectrons;
   iEvent.getByToken(electronTag_, pElectrons);

   if (pElectrons->size() < 1) return 0;

   int CountElectron = 0;
   for(edm::View<pat::Electron>::const_iterator iElectron=pElectrons->begin(); iElectron!=pElectrons->end(); ++iElectron)
   {
       if (iElectron->electronID(relIdName_) == passRelId_)
       {
           CountElectron++;
           electronColl->push_back(*iElectron);
       }
   }

   if (CountElectron >= 1)
   {
       iEvent.put(std::move(electronColl));
       return true;
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ElectronSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ElectronSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ElectronSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ElectronSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ElectronSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ElectronSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronSelector);
