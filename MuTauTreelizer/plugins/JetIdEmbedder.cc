#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

class JetIdEmbedder : public edm::stream::EDProducer<> {
public:
  JetIdEmbedder(const edm::ParameterSet& pset);
  virtual ~JetIdEmbedder(){}
  void produce(edm::Event& evt, const edm::EventSetup& es);
private:
  std::string puDisc_;
  edm::EDGetTokenT<edm::View<pat::Jet> > srcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ditau2017v1Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> ditau2017MDv1Token_;
};

JetIdEmbedder::JetIdEmbedder(const edm::ParameterSet& pset):
  puDisc_(pset.exists("discriminator") ? pset.getParameter<std::string>("discriminator") : "pileupJetId:fullDiscriminant"),
  srcToken_(consumes<edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("slimmedJetTag")))
{
  if (pset.exists("ditau2017v1")) { ditau2017v1Token_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ditau2017v1")); }
  if (pset.exists("ditau2017MDv1")) { ditau2017MDv1Token_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ditau2017MDv1")); }
  produces<pat::JetCollection>();
}

void JetIdEmbedder::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::unique_ptr<pat::JetCollection> output(new pat::JetCollection);

  edm::Handle<edm::View<pat::Jet> > input;
  evt.getByToken(srcToken_, input);

  edm::Handle<edm::ValueMap<float> > ditau2017v1;
  bool ditau2017v1Valid = evt.getByToken(ditau2017v1Token_, ditau2017v1);

  edm::Handle<edm::ValueMap<float> > ditau2017MDv1;
  bool ditau2017MDv1Valid = evt.getByToken(ditau2017MDv1Token_, ditau2017MDv1);

  output->reserve(input->size());
  for (size_t i = 0; i < input->size(); ++i) {
    pat::Jet jet = input->at(i);

    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    bool loose = true;
    bool tight = true;
    bool tightLepVeto = true;
    if (std::abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.99) {
        loose = false;
      }
      if (jet.neutralHadronEnergyFraction() >= 0.90) {
        tight = false;
        tightLepVeto = false;
      }

      if (jet.neutralEmEnergyFraction() >= 0.99) {
        loose = false;
      }
      if (jet.neutralEmEnergyFraction() >= 0.90) {
        tight = false;
        tightLepVeto = false;
      }

      if (jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1){
        loose = false;
        tight = false;
        tightLepVeto = false;
      }

      if (jet.muonEnergyFraction() >= 0.8)
        {
          tightLepVeto = false;
        }

      if (std::abs(jet.eta()) < 2.4) {
        if (jet.chargedHadronEnergyFraction() <= 0) {
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.chargedHadronMultiplicity() <= 0) {
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.chargedEmEnergyFraction() >= 0.99) {
          loose = false;
          tight = false;
        }
        if (jet.chargedEmEnergyFraction() >= 0.90) {
          tightLepVeto = false;
        }
      }
    }
    if (std::abs(jet.eta()) > 2.7 && std::abs(jet.eta()) <= 3.0) {
      if (jet.neutralEmEnergyFraction() >= 0.90) {
        loose = false;
        tight = false;
      }
      if (jet.neutralMultiplicity()<=2) {
        loose = false;
        tight = false;
      }
    }
    if (std::abs(jet.eta()) > 3.0) {
      if (jet.neutralEmEnergyFraction() >= 0.90) {
        loose = false;
        tight = false;
      }
      if (jet.neutralMultiplicity()<=10) {
        loose = false;
        tight = false;
      }
    }
    jet.addUserInt("idLoose", loose);
    jet.addUserInt("idTight", tight);
    jet.addUserInt("idTightLepVeto", tightLepVeto);

    // Pileup discriminant
    bool passPU = true;
    float jpumva = jet.userFloat(puDisc_);
    if(jet.pt() > 20)
      {
        if(fabs(jet.eta()) > 3.)
          {
            if(jpumva <= -0.45) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.75)
          {
            if(jpumva <= -0.55) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.5)
          {
            if(jpumva <= -0.6) passPU = false;
          }
        else if(jpumva <= -0.63) passPU = false;
      }
    else
      {
        if(fabs(jet.eta()) > 3.)
          {
            if(jpumva <= -0.95) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.75)
          {
            if(jpumva <= -0.94) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.5)
          {
            if(jpumva <= -0.96) passPU = false;
          }
        else if(jpumva <= -0.95) passPU = false;
      }

    jet.addUserInt("puID", passPU);
    float ditau2017v1Value = 0;
    float ditau2017MDv1Value = 0;
    edm::Ref<edm::View<pat::Jet> > jRef(input,i);
    if (ditau2017v1Valid) ditau2017v1Value = (*ditau2017v1)[jRef];
    if (ditau2017MDv1Valid) ditau2017MDv1Value = (*ditau2017MDv1)[jRef];
    jet.addUserFloat("ditau2017v1",ditau2017v1Value);
    jet.addUserFloat("ditau2017MDv1",ditau2017MDv1Value);
    output->push_back(jet);
  }

  evt.put(std::move(output));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetIdEmbedder);
