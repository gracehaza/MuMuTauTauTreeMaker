// -*- C++ -*-
//
// Package:    testED/TauHadSelector
// Class:      TauHadSelector
// 
/**\class TauHadSelector TauHadSelector.cc testED/TauHadSelector/plugins/TauHadSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Mon, 29 Apr 2019 12:55:27 GMT
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

#include "DataFormats/PatCandidates/interface/Tau.h"
#include <math.h>
//
// class declaration
//

class TauHadSelector : public edm::stream::EDFilter<> {
   public:
      explicit TauHadSelector(const edm::ParameterSet&);
      ~TauHadSelector();

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
      edm::EDGetTokenT<edm::View<pat::Tau>> tauTag_;
      std::vector<std::string> tauDiscriminatorTag_;
      double pTMin_;
      double etaMax_;
      bool passDiscriminator_;
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
TauHadSelector::TauHadSelector(const edm::ParameterSet& iConfig):
    tauTag_(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("tauTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Tau>>();
   tauDiscriminatorTag_ = iConfig.getParameter<std::vector<std::string>>("tauDiscriminatorTag");
   passDiscriminator_ = iConfig.getParameter<bool>("passDiscriminator");
   pTMin_ = iConfig.getParameter<double>("pTMin");
   etaMax_ = iConfig.getParameter<double>("etaMax");
}


TauHadSelector::~TauHadSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TauHadSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Tau>> tauColl = std::make_unique<std::vector<pat::Tau>>();

   edm::Handle<edm::View<pat::Tau>> pTaus;
   iEvent.getByToken(tauTag_, pTaus);

   if (pTaus->size() < 1) return 0;

   int CountTau=0;
   for (edm::View<pat::Tau>::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); ++iTau)
   {
       bool passIso = false;
       for (unsigned int i=0; i<tauDiscriminatorTag_.size(); i++)
       {
           if (iTau->tauID(tauDiscriminatorTag_[i])>0.5) passIso = true;
       } // end for loop on tau discriminator

       if (((passIso && passDiscriminator_) || (!passIso && !passDiscriminator_)) && 
          ((etaMax_ == -1.0) || (fabs((iTau)->eta()) < etaMax_)) &&
          ((pTMin_ == -1.0)  || ((iTau)->pt() > pTMin_)))
       {
           tauColl->push_back(*iTau);
           CountTau++;
       } // end if tau decay mode finding & pt/eta cut
   } // end for loop on tau candidates

   if (CountTau >= 1)
   {
       iEvent.put(std::move(tauColl));
       return true;
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TauHadSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TauHadSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TauHadSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TauHadSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TauHadSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TauHadSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauHadSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauHadSelector);
