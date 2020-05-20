// -*- C++ -*-
//
// Package:    testED/JetSelector
// Class:      JetSelector
// 
/**\class JetSelector JetSelector.cc testED/JetSelector/plugins/JetSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Fri, 24 May 2019 13:54:34 GMT
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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include <string>
#include <math.h>
//
// class declaration
//

class JetSelector : public edm::stream::EDFilter<> {
   public:
      explicit JetSelector(const edm::ParameterSet&);
      ~JetSelector();

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
      edm::EDGetTokenT<edm::View<pat::Jet>> jetTag_;
      std::string jetIdName_;
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
JetSelector::JetSelector(const edm::ParameterSet& iConfig):
    jetTag_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Jet>>();
   jetIdName_ = iConfig.getParameter<std::string>("jetIdName");
   etaCut_ = iConfig.getParameter<double>("etaCut");
   ptCut_ = iConfig.getParameter<double>("ptCut");
   transform(jetIdName_.begin(), jetIdName_.end(), jetIdName_.begin(), ::toupper);
}


JetSelector::~JetSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JetSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Jet>> jetColl = std::make_unique<std::vector<pat::Jet>>();
   edm::Handle<edm::View<pat::Jet>> pJets;
   iEvent.getByToken(jetTag_, pJets);

   if (pJets->size() < 1) 
   {
       iEvent.put(std::move(jetColl));
       return true; // we also select the events containing zero jet !
   }

   int CountJet = 0;
   for(edm::View<pat::Jet>::const_iterator iJet=pJets->begin(); iJet!=pJets->end(); ++iJet)
   {
       if ((iJet->pt() > ptCut_) && (fabs(iJet->eta()) < etaCut_) && (iJet->isPFJet()))
       {
           // -- reference: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h ---
           double jetEnergyUncorrected = iJet->chargedHadronEnergy() + iJet->neutralHadronEnergy() + iJet->photonEnergy() + iJet->electronEnergy() + iJet->muonEnergy() + iJet->HFEMEnergy();
           double muonFraction = (jetEnergyUncorrected > 0. ? (iJet->muonEnergy()/jetEnergyUncorrected) : 1.0);
           if ((jetIdName_ == "TIGHTLEPVETO") && 
               (iJet->chargedHadronEnergyFraction()>0) && 
               (iJet->neutralHadronEnergyFraction()<0.9) &&
               (iJet->chargedEmEnergyFraction()<0.8) &&
               (iJet->neutralEmEnergyFraction()<0.9) &&
               (iJet->chargedMultiplicity()>0) &&
               (iJet->numberOfDaughters()>1) &&
               (muonFraction < 0.8))
           {
               CountJet++;
               jetColl->push_back(*iJet);
           }

           else if ((jetIdName_ == "TIGHT") &&
                    (iJet->chargedHadronEnergyFraction()>0) &&
                    (iJet->neutralHadronEnergyFraction()<0.9) &&
                    (iJet->neutralEmEnergyFraction()<0.9) &&
                    (iJet->chargedMultiplicity()>0) &&
                    (iJet->numberOfDaughters()>1))
           {
               CountJet++;
               jetColl->push_back(*iJet);
           }
       } // end if jet candidate pass pt/eta cut and isPFJet
   } // end loop over all jet candidates

   iEvent.put(std::move(jetColl));
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
JetSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
JetSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
JetSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
JetSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
JetSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
JetSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(JetSelector);
