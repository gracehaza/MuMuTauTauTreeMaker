// -*- C++ -*-
//
// Package:    testED/ElectronCandSelector
// Class:      ElectronCandSelector
// 
/**\class ElectronCandSelector ElectronCandSelector.cc testED/ElectronCandSelector/plugins/ElectronCandSelector.cc

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
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include <string>
#include <math.h>
//
// class declaration
//

class ElectronCandSelector : public edm::stream::EDFilter<> {
   public:
      explicit ElectronCandSelector(const edm::ParameterSet&);
      ~ElectronCandSelector();

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
      edm::EDGetTokenT<double> rhoTag_;
      std::string relIdName_;
      bool passRelIso_;
      double etaCut_;
      double ptCut_;
      EffectiveAreas effectiveAreas_;
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
ElectronCandSelector::ElectronCandSelector(const edm::ParameterSet& iConfig):
    electronTag_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronTag"))),
    rhoTag_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"))),
    effectiveAreas_((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath())
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Electron>>();
   relIdName_ = iConfig.getParameter<std::string>("relIdName");
   passRelIso_ = iConfig.getParameter<bool>("passRelIso");
   etaCut_ = iConfig.getParameter<double>("etaCut");
   ptCut_ = iConfig.getParameter<double>("ptCut");
}


ElectronCandSelector::~ElectronCandSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ElectronCandSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Electron>> electronColl = std::make_unique<std::vector<pat::Electron>>();
   edm::Handle<edm::View<pat::Electron>> pElectrons;
   iEvent.getByToken(electronTag_, pElectrons);

   edm::Handle<double> pRho;
   iEvent.getByToken(rhoTag_, pRho);

   if (pElectrons->size() < 1) 
   {
       iEvent.put(std::move(electronColl));
       return true; // In this case for control region & fake rate study, the events without electrons are also kept
   }

   int CountElectron = 0;
   if (relIdName_ != "Veto" && relIdName_ != "Loose" && relIdName_ != "Medium" && relIdName_ != "Tight")
   {
       std::cout << "Error: Unknown electron ID type! Please verify your option." << std::endl;
   }

   for(edm::View<pat::Electron>::const_iterator iElectron=pElectrons->begin(); iElectron!=pElectrons->end(); ++iElectron)
   {
       bool isVeto = false;
       bool isLoose = false;
       bool isMedium = false;
       bool isTight = false;

       // ====== implement impact parameter cuts =======
       // reference: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria_for_V

       // ---  full5x5_sigmaIetaIeta ---
       double sigmaIetaIeta = iElectron->full5x5_sigmaIetaIeta();

       // --- fabs(dEtaSeed) ---
       double dEtaSeed = fabs(iElectron->superCluster().isNonnull() && iElectron->superCluster()->seed().isNonnull() ? iElectron->deltaEtaSuperClusterTrackAtVtx() - iElectron->superCluster()->eta() + iElectron->superCluster()->seed()->eta() : std::numeric_limits<float>::max()); 
       
       // --- fabs(dPhiIn) ---
       double dPhiIn = fabs(iElectron->deltaPhiSuperClusterTrackAtVtx());
       
       // --- variables for H/E cuts ---
       double HoE = iElectron->hadronicOverEm();
       double rho = pRho.isValid() ? (*pRho) : 0; 
       double energy = iElectron->superCluster()->energy();

       // --- variables for relIsoWithEffectiveArea ---
       double chad = iElectron->pfIsolationVariables().sumChargedHadronPt;
       double nhad = iElectron->pfIsolationVariables().sumNeutralHadronEt;
       double pho = iElectron->pfIsolationVariables().sumPhotonEt;
       double elePt = iElectron->pt();
       double eleEta = iElectron->superCluster()->eta();
       double eArea = effectiveAreas_.getEffectiveArea(fabs(eleEta));
       double relIsoWithEffectiveArea = (chad + std::max(0.0, nhad + pho - rho*eArea)) / elePt;

       // --- variables for fabs(1/E-1/p) ---
       double eInverseMinusPInverse = fabs(1.0 - iElectron->eSuperClusterOverP())*(1.0/iElectron->ecalEnergy());

       // --- expected missing inner hits ---
       int mHits = iElectron->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);

       // --- pass conversion veto ---
       bool isPassConVeto = iElectron->passConversionVeto();

       // ========= select electrons in different cut-based ID accordingly ==========
       if (fabs(eleEta) <= 1.479)
       {
           isVeto = (sigmaIetaIeta < 0.0126) && 
                    (dEtaSeed < 0.00463) && 
                    (dPhiIn < 0.148) && 
                    (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
                    (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.198 + 0.506/elePt)) &&
                    (eInverseMinusPInverse < 0.209) &&
                    (mHits <= 2) &&
                    (isPassConVeto == true);

           isLoose = (sigmaIetaIeta < 0.0112) &&
                     (dEtaSeed < 0.00377) &&
                     (dPhiIn < 0.0884) &&
                     (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
                     (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.112 + 0.506/elePt)) &&
                     (eInverseMinusPInverse < 0.193) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);

           isMedium = (sigmaIetaIeta < 0.0106) &&
                      (dEtaSeed < 0.0032) &&
                      (dPhiIn < 0.0547) &&
                      (HoE < 0.046 + 1.16/energy + 0.0324*rho/energy) &&
                      (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.0478 + 0.506/elePt)) &&
                      (eInverseMinusPInverse < 0.184) &&
                      (mHits <= 1) &&
                      (isPassConVeto == true);

           isTight = (sigmaIetaIeta < 0.0104) &&
                     (dEtaSeed < 0.00255) &&
                     (dPhiIn < 0.022) &&
                     (HoE < 0.026 + 1.15/energy + 0.0324*rho/energy) &&
                     (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.0287 + 0.506/elePt)) &&
                     (eInverseMinusPInverse < 0.159) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);

           if(((relIdName_ == "Veto" && isVeto) || (relIdName_ == "Loose" && isLoose) || (relIdName_ == "Medium" && isMedium) || (relIdName_ == "Tight" && isTight)) && (elePt > ptCut_) && (fabs(eleEta) < etaCut_))
           {
               CountElectron++;
               electronColl->push_back(*iElectron);
           } //endif (customized ID && ptCut && etaCut)

       }// endif (fabs(eleEta) <= 1.479)

       else{

           isVeto = (sigmaIetaIeta < 0.0457) && 
                    (dEtaSeed < 0.00814) && 
                    (dPhiIn < 0.19) && 
                    (HoE < 0.05 + 2.54/energy + 0.183*rho/energy) &&
                    (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.203 + 0.963/elePt)) &&
                    (eInverseMinusPInverse < 0.132) &&
                    (mHits <= 3) &&
                    (isPassConVeto == true);

           isLoose = (sigmaIetaIeta < 0.0425) &&
                     (dEtaSeed < 0.00674) &&
                     (dPhiIn < 0.169) &&
                     (HoE < 0.0441 + 2.54/energy + 0.183*rho/energy) &&
                     (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.108 + 0.963/elePt)) &&
                     (eInverseMinusPInverse < 0.111) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);

           isMedium = (sigmaIetaIeta < 0.0387) &&
                      (dEtaSeed < 0.00632) &&
                      (dPhiIn < 0.0394) &&
                      (HoE < 0.0275 + 2.52/energy + 0.183*rho/energy) &&
                      (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.0658 + 0.963/elePt)) &&
                      (eInverseMinusPInverse < 0.0721) &&
                      (mHits <= 1) &&
                      (isPassConVeto == true);

           isTight = (sigmaIetaIeta < 0.0353) &&
                     (dEtaSeed < 0.00501) &&
                     (dPhiIn < 0.0236) &&
                     (HoE < 0.0188 + 2.06/energy + 0.183*rho/energy) &&
                     (passRelIso_ == false || (passRelIso_ == true && relIsoWithEffectiveArea < 0.0445 + 0.963/elePt)) &&
                     (eInverseMinusPInverse < 0.0197) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);

           if(((relIdName_ == "Veto" && isVeto) || (relIdName_ == "Loose" && isLoose) || (relIdName_ == "Medium" && isMedium) || (relIdName_ == "Tight" && isTight)) && (elePt > ptCut_) && (fabs(eleEta) < etaCut_))
           {
               CountElectron++;
               electronColl->push_back(*iElectron);
           } //endif (customized ID && ptCut && etaCut)

       } // end else (fabs(eleEta) > 1.479)

   } // end for loop for pat::Electron candidates

   iEvent.put(std::move(electronColl));
   return true;

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ElectronCandSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ElectronCandSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ElectronCandSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ElectronCandSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ElectronCandSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ElectronCandSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronCandSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronCandSelector);
