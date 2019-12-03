// -*- C++ -*-
//
// Package:    GenMuTauTreelizer/GenTauHadCandSelector
// Class:      GenTauHadCandSelector
// 
/**\class GenTauHadCandSelector GenTauHadCandSelector.cc GenMuTauTreelizer/GenTauHadCandSelector/plugins/GenTauHadCandSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 26 Nov 2019 16:31:19 GMT
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
#include "TLorentzVector.h"
#include <math.h>

using namespace std;
//
// class declaration
//

class GenTauHadCandSelector : public edm::stream::EDFilter<> {
   public:
      explicit GenTauHadCandSelector(const edm::ParameterSet&);
      ~GenTauHadCandSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void checkTauDecayMode(const reco::Candidate*, std::vector<const reco::Candidate*>&);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

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
GenTauHadCandSelector::GenTauHadCandSelector(const edm::ParameterSet& iConfig):
    genParticleTag_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<reco::GenParticle>>();
   ptCut_ = iConfig.getParameter<double>("ptCut");
   etaCut_ = iConfig.getParameter<double>("etaCut");
}


GenTauHadCandSelector::~GenTauHadCandSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool GenTauHadCandSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<reco::GenParticle>> pGenParticles;
   std::unique_ptr<std::vector<reco::GenParticle>> tauHadColl = std::make_unique<std::vector<reco::GenParticle>>();

   iEvent.getByToken(genParticleTag_, pGenParticles);

   if(pGenParticles->size() < 1)
   {
       iEvent.put(std::move(tauHadColl));
       return true;
   }

   int CountTauHad = 0;

   for(edm::View<reco::GenParticle>::const_iterator iParticle=pGenParticles->begin(); iParticle!=pGenParticles->end(); ++iParticle)
   {
       if (fabs(iParticle->pdgId()) == 15 && fabs(iParticle->mother()->pdgId()) != 15)
       {
           int nDaughters = iParticle->numberOfDaughters();
           for (int iDaughter=0; iDaughter<nDaughters; iDaughter++)
           {
               int dauId = iParticle->daughter(iDaughter)->pdgId();
               std::vector<const reco::Candidate*> tauHadCand;
               tauHadCand.clear();

               if (fabs(dauId) == 15)
               {
                   checkTauDecayMode(iParticle->daughter(iDaughter), tauHadCand);

                   if (tauHadCand.size() > 0 && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
                   {
                       CountTauHad++;
                       tauHadColl->push_back(*iParticle);
                       break;
                   } // end if tauHadCand vector is filled
               } // end if daughter particle is tau (for exclusing tau->tau+gamma->tau+gamma+gamma ... -> mu/ele + ngammas)

               else if (fabs(dauId) != 11 && fabs(dauId) != 12 && fabs(dauId) != 13 && fabs(dauId) != 14 && fabs(dauId) != 16 && fabs(dauId) != 22 && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
               {
                   CountTauHad++;
                   tauHadColl->push_back(*iParticle);
                   break;
               } // else if the daughters are not mu/ele or neutrinos (leptonic decay of tau)
           } // end for loop on the daughter particles
       } // end if tau_h candidates
   }// end for loop on genParticles

   iEvent.put(std::move(tauHadColl));
   return true;

}

// ------------ regressive function for excluding leptonic decay of taus that have several final state radiations before decay -----------------
void GenTauHadCandSelector::checkTauDecayMode(const reco::Candidate* inputDaughter, std::vector<const reco::Candidate*>& daughterCand)
{
    bool tauLepDecayVeto = false;
    int nGrandDaughters = inputDaughter->numberOfDaughters();
    for (int iGrandDaughter = 0; iGrandDaughter < nGrandDaughters; iGrandDaughter++)
    {
        int grandDauId = inputDaughter->daughter(iGrandDaughter)->pdgId(); 
        if (fabs(grandDauId) == 15)
        {
            checkTauDecayMode(inputDaughter->daughter(iGrandDaughter), daughterCand);
        } // end if granddaughter is still tau (FSR)

        else{
            if (fabs(grandDauId) != 11 && fabs(grandDauId) != 12 && fabs(grandDauId) != 13 && fabs(grandDauId) != 14)
            {
                daughterCand.push_back(inputDaughter->daughter(iGrandDaughter));
            } // end if not leptonic decay of input daughter tau

            else{
                tauLepDecayVeto = true;
                break;
            } // end else if leptonic decay of input daughter tau
        } // end else granddaughter not tau
    } // end for loop on granddaughters

    if (tauLepDecayVeto == true)
    {
        daughterCand.clear();
    } // end if tauLepDecayVeto == true
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenTauHadCandSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenTauHadCandSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenTauHadCandSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenTauHadCandSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenTauHadCandSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenTauHadCandSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenTauHadCandSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenTauHadCandSelector);
