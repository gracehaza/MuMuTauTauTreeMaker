// -*- C++ -*-
//
// Package:    GenMuTauTreelizer/GenTauEleCandSelector
// Class:      GenTauEleCandSelector
// 
/**\class GenTauEleCandSelector GenTauEleCandSelector.cc GenMuTauTreelizer/GenTauEleCandSelector/plugins/GenTauEleCandSelector.cc

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
#include <math.h>

using namespace std;
//
// class declaration
//

class GenTauEleCandSelector : public edm::stream::EDFilter<> {
   public:
      explicit GenTauEleCandSelector(const edm::ParameterSet&);
      ~GenTauEleCandSelector();

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
GenTauEleCandSelector::GenTauEleCandSelector(const edm::ParameterSet& iConfig):
    genParticleTag_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesTag")))
{
   //now do what ever initialization is needed
   produces<std::vector<reco::GenParticle>>();
   ptCut_ = iConfig.getParameter<double>("ptCut");
   etaCut_ = iConfig.getParameter<double>("etaCut");
}


GenTauEleCandSelector::~GenTauEleCandSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool GenTauEleCandSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<reco::GenParticle>> pGenParticles;
   std::unique_ptr<std::vector<reco::GenParticle>> tauEleColl = std::make_unique<std::vector<reco::GenParticle>>();

   iEvent.getByToken(genParticleTag_, pGenParticles);

   if(pGenParticles->size() < 1)
   {
       iEvent.put(std::move(tauEleColl));
       return true;
   }

   int CountTauEle = 0;

   for(edm::View<reco::GenParticle>::const_iterator iParticle=pGenParticles->begin(); iParticle!=pGenParticles->end(); ++iParticle)
   {
       if (fabs(iParticle->pdgId()) == 15 && fabs(iParticle->mother()->pdgId()) != 15)
       {
           int nDaughters = iParticle->numberOfDaughters();
           for (int iDaughter=0; iDaughter<nDaughters; iDaughter++)
           {
               int dauId = iParticle->daughter(iDaughter)->pdgId();
               std::vector<const reco::Candidate*> tauEleCand;
               tauEleCand.clear();

               if (fabs(dauId) == 15)
               {
                   checkTauDecayMode(iParticle->daughter(iDaughter), tauEleCand);

                   if (tauEleCand.size() > 0 && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
                   {
                       CountTauEle++;
                       tauEleColl->push_back(*iParticle);
                       break;
                   } // end if tauEleCand vector is filled
               } // end if daughter particle is tau (for exclusing tau->tau+gamma->tau+gamma+gamma ... -> hadrons + ngammas) -- FSR

               else if ((fabs(dauId) == 11 || fabs(dauId) == 12) && iParticle->pt() >= ptCut_ && fabs(iParticle->eta()) <= etaCut_)
               {
                   CountTauEle++;
                   tauEleColl->push_back(*iParticle);
                   break;
               } // else if the daughters are ele or neutrinos (leptonic decay of tau)
           } // end for loop on the daughter particles
       } // end if tau_e candidates
   }// end for loop on genParticles

   iEvent.put(std::move(tauEleColl));
   return true;

}

// ------------ regressive function for excluding hadronic decay of taus that have several final state radiations before decay -----------------
void GenTauEleCandSelector::checkTauDecayMode(const reco::Candidate* inputDaughter, std::vector<const reco::Candidate*>& daughterCand)
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
            if (fabs(grandDauId) == 11 || fabs(grandDauId) == 12)
            {
                daughterCand.push_back(inputDaughter->daughter(iGrandDaughter));
            } // end if leptonic decay of input daughter tau (e-nu)

            else if (fabs(grandDauId) == 13 || fabs(grandDauId) == 14 || fabs(grandDauId) == 16 || fabs(grandDauId) == 22)
            {
                continue;
            } // end if a tau neutrino is found, since the tau neutrino is produced in both leptonic and hadronic decays; or if mu-nu from tau decay; or photon from FSR

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
GenTauEleCandSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenTauEleCandSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenTauEleCandSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenTauEleCandSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenTauEleCandSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenTauEleCandSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenTauEleCandSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenTauEleCandSelector);
