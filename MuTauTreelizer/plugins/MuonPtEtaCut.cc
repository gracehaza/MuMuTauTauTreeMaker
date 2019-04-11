// -*- C++ -*-
//
// Package:    MuMuChannel/MuonPtEtaCut
// Class:      MuonPtEtaCut
// 
/**\class MuonPtEtaCut MuonPtEtaCut.cc MuMuChannel/MuonPtEtaCut/plugins/MuonPtEtaCut.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Wed, 10 Apr 2019 16:09:53 GMT
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

class MuonPtEtaCut : public edm::stream::EDFilter<> {
   public:
      explicit MuonPtEtaCut(const edm::ParameterSet&);
      ~MuonPtEtaCut();

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
      double Eta_;
      double Pt_;
      unsigned int minNumObjsToPassFilter_;
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
MuonPtEtaCut::MuonPtEtaCut(const edm::ParameterSet& iConfig):
    muonTag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   minNumObjsToPassFilter_ = iConfig.getParameter<unsigned int>("minNumObjsToPassFilter");
   Eta_ = iConfig.getParameter<double>("Eta");
   Pt_ = iConfig.getParameter<double>("Pt");

}


MuonPtEtaCut::~MuonPtEtaCut()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuonPtEtaCut::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   unsigned int nPassingMuons=0;
   edm::Handle<edm::View<pat::Muon>> recoObjs;
   iEvent.getByToken(muonTag_, recoObjs);

   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   for (edm::View<pat::Muon>::const_iterator iRecoObj = recoObjs->begin(); iRecoObj != recoObjs->end();  ++iRecoObj)
   {
       if((fabs(iRecoObj->eta())<Eta_ || Eta_==-1) && iRecoObj->pt()>Pt_)
       {
           muonColl->push_back(*iRecoObj);
           nPassingMuons++;
       }
   }

   iEvent.put(std::move(muonColl), "muonColls");

   return (nPassingMuons >= minNumObjsToPassFilter_);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MuonPtEtaCut::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuonPtEtaCut::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuonPtEtaCut::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MuonPtEtaCut::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonPtEtaCut::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonPtEtaCut::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonPtEtaCut::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonPtEtaCut);
