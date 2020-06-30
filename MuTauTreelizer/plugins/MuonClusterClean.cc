// -*- C++ -*-
//
// Package:    slimmedClean/MuonClusterClean
// Class:      MuonClusterClean
// 
/**\class MuonClusterClean MuonClusterClean.cc slimmedClean/MuonClusterClean/plugins/MuonClusterClean.cc

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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <math.h>

using namespace std;
//
// class declaration
//

class MuonClusterClean : public edm::stream::EDProducer<> {
   public:
      explicit MuonClusterClean(const edm::ParameterSet&);
      ~MuonClusterClean();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<edm::View<pat::Muon>> MuTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauTag;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>> ParticleFlowCandTag;
      std::string muonID;
      double ptCut;

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
MuonClusterClean::MuonClusterClean(const edm::ParameterSet& iConfig):
    MuTag(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuTag"))),
    TauTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauTag"))),
    ParticleFlowCandTag(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("ParticleFlowCandTag")))
{
   produces<std::vector<pat::Muon>>("slimmedMuonsCleaned").setBranchAlias("muonColl");
   produces<std::vector<pat::PackedCandidate>>("packedPFCandidatesMuonCleaned").setBranchAlias("pfCandColl");
   muonID = iConfig.getParameter<std::string>("muonID");
   ptCut = iConfig.getParameter<double>("ptCut");
   //now do what ever other initialization is needed
   transform(muonID.begin(), muonID.end(), muonID.begin(), ::toupper);
}


MuonClusterClean::~MuonClusterClean()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonClusterClean::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   std::unique_ptr<std::vector<pat::PackedCandidate>> pfCandColl = std::make_unique<std::vector<pat::PackedCandidate>>();

   std::vector<pat::Muon> cleanedMuons;
   cleanedMuons.clear();

   edm::Handle<edm::View<pat::Muon>> pMu;
   iEvent.getByToken(MuTag, pMu);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauTag, pTau);

   edm::Handle<edm::View<pat::PackedCandidate>> pPFCand;
   iEvent.getByToken(ParticleFlowCandTag, pPFCand);

   for(edm::View<pat::Muon>::const_iterator iMuon = pMu->begin(); iMuon != pMu->end(); ++iMuon)
   {
       bool findMatchedMuon = false;
       bool goodGlob = iMuon->isGlobalMuon() && iMuon->globalTrack()->normalizedChi2() < 3 && iMuon->combinedQuality().chi2LocalPosition < 12 && iMuon->combinedQuality().trkKink < 20;
       bool isMedium = muon::isLooseMuon(*iMuon) && iMuon->innerTrack()->validFraction() > 0.8 && muon::segmentCompatibility(*iMuon) > (goodGlob ? 0.303 : 0.451);
       bool isLoose = iMuon->isPFMuon() && (iMuon->isGlobalMuon() || iMuon->isTrackerMuon());
       bool isTight = iMuon->isGlobalMuon() && iMuon->isPFMuon() && iMuon->globalTrack()->normalizedChi2() < 10 && iMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && iMuon->numberOfMatchedStations() > 1 && fabs(iMuon->muonBestTrack()->dxy()) < 0.2 && fabs(iMuon->muonBestTrack()->dz()) < 0.5 && iMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
       bool passMuonID = (isLoose && muonID == "LOOSE") || (isMedium && muonID == "MEDIUM") || (isTight && muonID == "TIGHT");

       for (edm::View<pat::Tau>::const_iterator iTau = pTau->begin(); iTau != pTau->end(); ++iTau)
       {
           if (deltaR(*iTau, *iMuon) < 0.8 && passMuonID && iMuon->pt() > ptCut) 
           {
               findMatchedMuon = true;
           } // end if find matched muon
       } // end loop on taus

       if (findMatchedMuon)
       {
           cleanedMuons.push_back(*iMuon);
       } // end if findMatchedMuon == true

       else{
           muonColl->push_back(*iMuon);
       } // end if findMatchedMuon == false
   } // end loop on muons

   for(edm::View<pat::PackedCandidate>::const_iterator iCand = pPFCand->begin(); iCand != pPFCand->end(); ++iCand)
   {
       bool findMatchedPFMuon = false;
       if (fabs(iCand->pdgId()) != 13) continue;
       for (unsigned int iMuon = 0; iMuon < cleanedMuons.size(); iMuon++)
       {
           if (deltaR(*iCand, cleanedMuons.at(iMuon)) < 0.01)
           {
               findMatchedPFMuon = true;
           } // end if deltaR requirement
       } // end loop on cleaned muons

       if (!findMatchedPFMuon)
       {
           pfCandColl->push_back(*iCand);
       } // end if findMatchedPFMuon == true
   } // end loop on PackedPFCandidate

   iEvent.put(std::move(muonColl), "slimmedMuonsCleaned");
   iEvent.put(std::move(pfCandColl), "packedPFCandidatesMuonCleaned");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MuonClusterClean::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuonClusterClean::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuonClusterClean::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MuonClusterClean::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonClusterClean::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonClusterClean::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonClusterClean::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonClusterClean);
