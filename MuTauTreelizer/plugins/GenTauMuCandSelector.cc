// -*- C++ -*-
//
// Package:    GenMuTauTreelizer/GenTauMuCandSelector
// Class:      GenTauMuCandSelector
// 
/**\class GenTauMuCandSelector GenTauMuCandSelector.cc GenMuTauTreelizer/GenTauMuCandSelector/plugins/GenTauMuCandSelector.cc

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

class GenTauMuCandSelector : public edm::stream::EDFilter<> {
   public:
      explicit GenTauMuCandSelector(const edm::ParameterSet&);
      ~GenTauMuCandSelector();

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
GenTauMuCandSelector::GenTauMuCandSelector(const edm::ParameterSet& iConfig):
    genParticleTag_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<reco::GenParticle>>();
   ptCut_ = iConfig.getParameter<double>("ptCut");
   etaCut_ = iConfig.getParameter<double>("etaCut");
}


GenTauMuCandSelector::~GenTauMuCandSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool GenTauMuCandSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<reco::GenParticle>> pGenParticles;
   std::unique_ptr<std::vector<reco::GenParticle>> tauMuColl = std::make_unique<std::vector<reco::GenParticle>>();

   iEvent.getByToken(genParticleTag_, pGenParticles);

   if(pGenParticles->size() < 1)
   {
       iEvent.put(std::move(tauMuColl));
       return true;
   }

   int CountTauMu = 0;

   for(edm::View<reco::GenParticle>::const_iterator iParticle=pGenParticles->begin(); iParticle!=pGenParticles->end(); ++iParticle)
   {
       if (fabs(iParticle->pdgId()) == 15 && fabs(iParticle->mother()->pdgId()) != 15)
       {
           int nDaughters = iParticle->numberOfDaughters();
           for (int iDaughter=0; iDaughter<nDaughters; iDaughter++)
           {
               int dauId = iParticle->daughter(iDaughter)->pdgId();
               std::vector<const reco::Candidate*> tauMuCand;
               tauMuCand.clear();

               if (fabs(dauId) == 15)
               {
                   checkTauDecayMode(iParticle->daughter(iDaughter), tauMuCand);

                   if (tauMuCand.size() > 0 && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
                   {
                       CountTauMu++;
                       tauMuColl->push_back(*iParticle);
                       break;
                   } // end if tauMuCand vector is filled
               } // end if daughter particle is tau (for exclusing tau->tau+gamma->tau+gamma+gamma ... -> hadrons + ngammas) -- FSR

               else if ((fabs(dauId) == 11 || fabs(dauId) == 12 || fabs(dauId) == 13 || fabs(dauId) == 14) && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
               {
                   CountTauMu++;
                   tauMuColl->push_back(*iParticle);
                   break;
               } // else if the daughters are mu/ele or neutrinos (leptonic decay of tau)
           } // end for loop on the daughter particles
       } // end if tau_mu candidates
   }// end for loop on genParticles

   iEvent.put(std::move(tauMuColl));
   return true;

}

// ------------ regressive function for excluding hadronic decay of taus that have several final state radiations before decay -----------------
void GenTauMuCandSelector::checkTauDecayMode(const reco::Candidate* inputDaughter, std::vector<const reco::Candidate*>& daughterCand)
{
    bool tauHadDecayVeto = false;
    int nGrandDaughters = inputDaughter->numberOfDaughters();
    for (int iGrandDaughter = 0; iGrandDaughter < nGrandDaughters; iGrandDaughter++)
    {
        int grandDauId = inputDaughter->daughter(iGrandDaughter)->pdgId(); 
        if (fabs(grandDauId) == 15)
        {
            checkTauDecayMode(inputDaughter->daughter(iGrandDaughter), daughterCand);
        } // end if granddaughter is still tau (FSR)

        else{
            if (fabs(grandDauId) == 11 || fabs(grandDauId) == 12 || fabs(grandDauId) == 13 || fabs(grandDauId) == 14)
            {
                daughterCand.push_back(inputDaughter->daughter(iGrandDaughter));
            } // end if leptonic decay of input daughter tau

            else if (fabs(grandDauId) == 16)
            {
                continue;
            } // end if a tau neutrino is found, since the tau neutrino is produced in both leptonic and hadronic decays

            else{
                tauHadDecayVeto = true;
                break;
            } // end else if leptonic decay of input daughter tau
        } // end else granddaughter not tau
    } // end for loop on granddaughters

    if (tauHadDecayVeto == true)
    {
        daughterCand.clear();
    } // end if tauHadDecayVeto == true
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenTauMuCandSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenTauMuCandSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenTauMuCandSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenTauMuCandSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenTauMuCandSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenTauMuCandSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenTauMuCandSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenTauMuCandSelector);
