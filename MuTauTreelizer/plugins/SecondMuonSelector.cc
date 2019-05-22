// -*- C++ -*-
//
// Package:    testED/SecondMuonSelector
// Class:      SecondMuonSelector
// 
/**\class SecondMuonSelector SecondMuonSelector.cc testED/SecondMuonSelector/plugins/SecondMuonSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 16 Apr 2019 14:46:19 GMT
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

class SecondMuonSelector : public edm::stream::EDFilter<> {
   public:
      explicit SecondMuonSelector(const edm::ParameterSet&);
      ~SecondMuonSelector();

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
      edm::EDGetTokenT<edm::View<pat::Muon>> muonTag_;
      edm::EDGetTokenT<edm::View<pat::Muon>> mu1Tag_;
      double relIsoCutVal_;
      bool passRelIso_;
      double dRCut_;
      bool oppositeSign_;
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
SecondMuonSelector::SecondMuonSelector(const edm::ParameterSet& iConfig):
    muonTag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag"))),
    mu1Tag_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("mu1Tag")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon>>();
   relIsoCutVal_ = iConfig.getParameter<double>("relIsoCutVal");
   passRelIso_ = iConfig.getParameter<bool>("passRelIso");
   dRCut_ = iConfig.getParameter<double>("dRCut");
   oppositeSign_ = iConfig.getParameter<bool>("oppositeSign");
}


SecondMuonSelector::~SecondMuonSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SecondMuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Muon>> muonColl = std::make_unique<std::vector<pat::Muon>>();
   std::unique_ptr<std::vector<pat::Muon>> secondMuonColl = std::make_unique<std::vector<pat::Muon>>();
   edm::Handle<edm::View<pat::Muon>> pMuons;
   iEvent.getByToken(muonTag_, pMuons);

   if (pMuons->size() < 1) return 0;
   edm::Handle<edm::View<pat::Muon>> pMu1;
   iEvent.getByToken(mu1Tag_, pMu1);
   pat::Muon mu1 = pMu1->at(0);

   int CountMuon=0;
   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end(); ++iMuon)
   {
       double dR = deltaR(*iMuon, mu1);
       if (dRCut_ > 0 && dR < dRCut_) continue;
       reco::MuonPFIsolation iso = iMuon->pfIsolationR04();
       double reliso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / iMuon->pt();
       if ((reliso < relIsoCutVal_ && passRelIso_) || (reliso > relIsoCutVal_ && !passRelIso_) || relIsoCutVal_ == -1)
       {
           if ((oppositeSign_ && mu1.pdgId() == (-1)* iMuon->pdgId()) || (!oppositeSign_ && mu1.pdgId() == iMuon->pdgId())){
              CountMuon+=1;
              muonColl->push_back(*iMuon);
           }
       }
   }

   // below we are going to select the second muon with smallest dR from mu1
   if (CountMuon >= 1)
   {
       double smallestDR = 99.0;
       pat::Muon smallestDRMuon;
       for (std::vector<pat::Muon>::iterator iMuon=muonColl->begin(); iMuon!=muonColl->end(); ++iMuon)
       {
           if ((deltaR(*iMuon, mu1) < 0.0001) && (fabs(iMuon->pt()-mu1.pt()) < 0.0001)) // exclude mu1 from the collection
               continue;

           else if (deltaR(*iMuon, mu1) < smallestDR)
           {
               smallestDR = deltaR(*iMuon, mu1);
               smallestDRMuon = *iMuon;
           }
       }

       secondMuonColl->push_back(smallestDRMuon);
       iEvent.put(std::move(secondMuonColl));
       return true;
   }

   return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SecondMuonSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SecondMuonSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
SecondMuonSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SecondMuonSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SecondMuonSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SecondMuonSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SecondMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SecondMuonSelector);
