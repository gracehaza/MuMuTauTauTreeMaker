// -*- C++ -*-
//
// Package:    slimmedClean/ElectronClusterClean
// Class:      ElectronClusterClean
// 
/**\class ElectronClusterClean ElectronClusterClean.cc slimmedClean/ElectronClusterClean/plugins/ElectronClusterClean.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Fri, 15 May 2020 17:45:22 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include <math.h>

using namespace std;
//
// class declaration
//

class ElectronClusterClean : public edm::stream::EDProducer<> {
   public:
      explicit ElectronClusterClean(const edm::ParameterSet&);
      ~ElectronClusterClean();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<edm::View<pat::Electron>> EleTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauTag;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>> ParticleFlowCandTag;
      std::string electronID;
      edm::EDGetTokenT<double> rhoTag;
      double ptCut;
      bool passRelIso;
      EffectiveAreas effectiveAreas;
      // ----------member data ---------------------------
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
ElectronClusterClean::ElectronClusterClean(const edm::ParameterSet& iConfig):
    EleTag(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("EleTag"))),
    TauTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauTag"))),
    ParticleFlowCandTag(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("ParticleFlowCandTag"))),
    rhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"))),
    effectiveAreas((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath())
{
   produces<std::vector<pat::Electron>>("slimmedElectronsCleaned").setBranchAlias("electronColl");
   produces<std::vector<pat::PackedCandidate>>("packedPFCandidatesElectronCleaned").setBranchAlias("pfCandColl");
   electronID = iConfig.getParameter<std::string>("electronID");
   ptCut = iConfig.getParameter<double>("ptCut");
   passRelIso = iConfig.getParameter<bool>("passRelIso");
   //now do what ever other initialization is needed
   transform(electronID.begin(), electronID.end(), electronID.begin(), ::toupper);
}


ElectronClusterClean::~ElectronClusterClean()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ElectronClusterClean::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Electron>> electronColl = std::make_unique<std::vector<pat::Electron>>();
   std::unique_ptr<std::vector<pat::PackedCandidate>> pfCandColl = std::make_unique<std::vector<pat::PackedCandidate>>();

   std::vector<pat::Electron> cleanedElectrons;
   cleanedElectrons.clear();

   edm::Handle<edm::View<pat::Electron>> pEle;
   iEvent.getByToken(EleTag, pEle);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauTag, pTau);

   edm::Handle<edm::View<pat::PackedCandidate>> pPFCand;
   iEvent.getByToken(ParticleFlowCandTag, pPFCand);

   edm::Handle<double> pRho;
   iEvent.getByToken(rhoTag, pRho);

   for(edm::View<pat::Electron>::const_iterator iElectron = pEle->begin(); iElectron != pEle->end(); ++iElectron)
   {
       bool findMatchedElectron = false;
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
       double eArea = effectiveAreas.getEffectiveArea(fabs(eleEta));
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
           isLoose = (sigmaIetaIeta < 0.0112) &&
                     (dEtaSeed < 0.00377) &&
                     (dPhiIn < 0.0884) &&
                     (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
                     (passRelIso == false || (passRelIso == true && relIsoWithEffectiveArea < 0.112 + 0.506/elePt)) &&
                     (eInverseMinusPInverse < 0.193) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);

           isMedium = (sigmaIetaIeta < 0.0106) &&
                      (dEtaSeed < 0.0032) &&
                      (dPhiIn < 0.0547) &&
                      (HoE < 0.046 + 1.16/energy + 0.0324*rho/energy) &&
                      (passRelIso == false || (passRelIso == true && relIsoWithEffectiveArea < 0.0478 + 0.506/elePt)) &&
                      (eInverseMinusPInverse < 0.184) &&
                      (mHits <= 1) &&
                      (isPassConVeto == true);

           isTight = (sigmaIetaIeta < 0.0104) &&
                     (dEtaSeed < 0.00255) &&
                     (dPhiIn < 0.022) &&
                     (HoE < 0.026 + 1.15/energy + 0.0324*rho/energy) &&
                     (passRelIso == false || (passRelIso == true && relIsoWithEffectiveArea < 0.0287 + 0.506/elePt)) &&
                     (eInverseMinusPInverse < 0.159) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);
       }// endif (fabs(eleEta) <= 1.479)

       else{
           isLoose = (sigmaIetaIeta < 0.0425) &&
                     (dEtaSeed < 0.00674) &&
                     (dPhiIn < 0.169) &&
                     (HoE < 0.0441 + 2.54/energy + 0.183*rho/energy) &&
                     (passRelIso == false || (passRelIso == true && relIsoWithEffectiveArea < 0.108 + 0.963/elePt)) &&
                     (eInverseMinusPInverse < 0.111) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);

           isMedium = (sigmaIetaIeta < 0.0387) &&
                      (dEtaSeed < 0.00632) &&
                      (dPhiIn < 0.0394) &&
                      (HoE < 0.0275 + 2.52/energy + 0.183*rho/energy) &&
                      (passRelIso == false || (passRelIso == true && relIsoWithEffectiveArea < 0.0658 + 0.963/elePt)) &&
                      (eInverseMinusPInverse < 0.0721) &&
                      (mHits <= 1) &&
                      (isPassConVeto == true);

           isTight = (sigmaIetaIeta < 0.0353) &&
                     (dEtaSeed < 0.00501) &&
                     (dPhiIn < 0.0236) &&
                     (HoE < 0.0188 + 2.06/energy + 0.183*rho/energy) &&
                     (passRelIso == false || (passRelIso == true && relIsoWithEffectiveArea < 0.0445 + 0.963/elePt)) &&
                     (eInverseMinusPInverse < 0.0197) &&
                     (mHits <= 1) &&
                     (isPassConVeto == true);
       } // end else (fabs(eleEta) > 1.479)

       bool passElectronID = (isLoose && electronID == "LOOSE") || (isMedium && electronID == "MEDIUM") || (isTight && electronID == "TIGHT");

       for (edm::View<pat::Tau>::const_iterator iTau = pTau->begin(); iTau != pTau->end(); ++iTau)
       {
           if (deltaR(*iTau, *iElectron) < 0.8 && passElectronID && elePt > ptCut)
           {
               findMatchedElectron = true;
           } // end if find matched electron
       } // end loop on taus

       if (findMatchedElectron)
       {
           cleanedElectrons.push_back(*iElectron);
       } // end if findMatchedElectron == true

       else{
           electronColl->push_back(*iElectron);
       } // end if findMatchedElectron == false
   } // end loop on electrons

   for(edm::View<pat::PackedCandidate>::const_iterator iCand = pPFCand->begin(); iCand != pPFCand->end(); ++iCand)
   {
       bool findMatchedPFElectron = false;
       if (fabs(iCand->pdgId()) != 11) 
       {
           pfCandColl->push_back(*iCand);
           continue;
       } // end if non-electron candidates

       for (unsigned int iElectron = 0; iElectron < cleanedElectrons.size(); iElectron++)
       {
           if (deltaR(*iCand, cleanedElectrons.at(iElectron)) < 0.01)
           {
               findMatchedPFElectron = true;
           } // end if deltaR requirement
       } // end loop on cleaned electrons

       if (!findMatchedPFElectron)
       {
           pfCandColl->push_back(*iCand);
       } // end if findMatchedPFElectron == true
   } // end loop on PackedPFCandidate

   iEvent.put(std::move(electronColl), "slimmedElectronsCleaned");
   iEvent.put(std::move(pfCandColl), "packedPFCandidatesElectronCleaned");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ElectronClusterClean::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ElectronClusterClean::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ElectronClusterClean::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ElectronClusterClean::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ElectronClusterClean::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ElectronClusterClean::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronClusterClean::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronClusterClean);
