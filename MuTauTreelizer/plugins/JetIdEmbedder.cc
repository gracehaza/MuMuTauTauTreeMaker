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
  if (pset.exists("DeepDiTau_boosted")) { DeepDiTau_boostedToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_boosted")); }
  if (pset.exists("DeepDiTau_boosted_massdeco")) {DeepDiTau_boosted_massdecoToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_boosted_massdeco")); }
  if (pset.exists("DeepDiTau_boosted_nolepton_charm")) { DeepDiTau_boosted_nolepton_charmToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_boosted_nolepton_charm")); }
  if (pset.exists("DeepDiTau_boosted_nolepton_charm_massdeco")) { DeepDiTau_boosted_nolepton_charm_massdecoToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_boosted_nolepton_charm_massdeco")); }
  if (pset.exists("DeepDiTau_boosted_nolepton_massdeco")) { DeepDiTau_boosted_nolepton_massdecoToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_boosted_nolepton_massdeco")); }
  if (pset.exists("DeepDiTau_boosted_nolepton")) { DeepDiTau_boosted_noleptonToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_boosted_nolepton")); }
  if (pset.exists("DeepDiTau_massdeco")) {DeepDiTau_massdecoToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_massdeco")); }
  if (pset.exists("DeepDiTau")) {DeepDiTauToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau")); }
  if (pset.exists("DeepDiTau_nolepton_charm_massdeco")) { DeepDiTau_nolepton_charm_massdecoToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_nolepton_charm_massdeco")); }
  if (pset.exists("DeepDiTau_nolepton_charm")) { DeepDiTau_nolepton_charmToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_nolepton_charm")); }
  if (pset.exists("DeepDiTau_nolepton_massdeco")) { DeepDiTau_nolepton_massdecoToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_nolepton_massdeco")); }
  if (pset.exists("DeepDiTau_nolepton")) { DeepDiTau_noleptonToken = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("DeepDiTau_nolepton")); }
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

  edm::Handle<edm::ValueMap<float> > DeepDiTau_boosted;
  bool DeepDiTau_boostedValid = evt.getByToken(DeepDiTau_boostedToken_, DeepDiTau_boosted);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_boosted_massdeco;
  bool DeepDiTau_boosted_massdecoValid = evt.getByToken(DeepDiTau_boosted_massdecoToken_, DeepDiTau_boosted_massdeco);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_boosted_nolepton_charm;
  bool DeepDiTau_boosted_nolepton_charmValid = evt.getByToken(DeepDiTau_boosted_nolepton_charmToken_, DeepDiTau_boosted_nolepton_charm);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_boosted_nolepton_charm_massdeco;
  bool DeepDiTau_boosted_nolepton_charm_massdecoValid = evt.getByToken(DeepDiTau_boosted_nolepton_charm_massdecoToken_, DeepDiTau_boosted_nolepton_charm_massdeco);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_boosted_nolepton_massdeco;
  bool DeepDiTau_boosted_nolepton_massdecoValid = evt.getByToken(DeepDiTau_boosted_nolepton_massdecoToken_, DeepDiTau_boosted_nolepton_massdeco);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_boosted_nolepton;
  bool DeepDiTau_boosted_noleptonValid = evt.getByToken(DeepDiTau_boosted_nolepton_Token_, DeepDiTau_boosted_nolepton);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_massdeco;
  bool DeepDiTau_massdecoValid = evt.getByToken(DeepDiTau_massdecoToken_, DeepDiTau_massdeco);

  edm::Handle<edm::ValueMap<float> > DeepDiTau;
  bool DeepDiTauValid = evt.getByToken(DeepDiTauToken_, DeepDiTau);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_nolepton_charm_massdeco;
  bool DeepDiTau_nolepton_charm_massdecoValid = evt.getByToken(DeepDiTau_nolepton_charm_massdecoToken_, DeepDiTau_nolepton_charm_massdeco);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_nolepton_charm;
  bool DeepDiTau_nolepton_charmValid = evt.getByToken(DeepDiTau_nolepton_charmToken_, DeepDiTau_nolepton_charm);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_nolepton_massdeco;
  bool DeepDiTau_nolepton_massdecoValid = evt.getByToken(DeepDiTau_nolepton_massdecoToken_, DeepDiTau_nolepton_massdeco);

  edm::Handle<edm::ValueMap<float> > DeepDiTau_nolepton;
  bool DeepDiTau_noleptonValid = evt.getByToken(DeepDiTau_noleptonToken_, DeepDiTau_nolepton);			


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
    float DeepDiTau_boostedValue= 0;
    float DeepDiTau_boosted_massdecoValue = 0;
    float DeepDiTau_boosted_nolepton_charmValue = 0;
    float DeepDiTau_boosted_nolepton_charm_massdecoValue = 0;
    float DeepDiTau_boosted_nolepton_massdecoValue = 0;
    float DeepDiTau_boosted_noleptonValue = 0;
    float DeepDiTau_massdecoValue = 0;
    float DeepDiTauValue = 0;
    float DeepDiTau_nolepton_charm_massdecoValue = 0;
    float DeepDiTau_nolepton_charmValue = 0;
    float DeepDiTau_nolepton_massdecoValue = 0;
    float DeepDiTau_noleptonValue = 0;

    edm::Ref<edm::View<pat::Jet> > jRef(input,i);
    if (ditau2017v1Valid) ditau2017v1Value = (*ditau2017v1)[jRef];
    if (ditau2017MDv1Valid) ditau2017MDv1Value = (*ditau2017MDv1)[jRef];
    if (DeepDiTau_boostedValid) DeepDiTau_boostedValue = (*DeepDiTau_boosted)[jRef];
    if (DeepDiTau_boosted_massdecoValid) DeepDiTau_boosted_massdecoValue = (*DeepDiTau_boosted_massdeco)[jRef];
    if (DeepDiTau_boosted_nolepton_charmValid) DeepDiTau_boosted_nolepton_charmValue = (*DeepDiTau_boosted_nolepton_charmValue)[jRef];
    if (DeepDiTau_boosted_nolepton_charm_massdecoValid) DeepDiTau_boosted_nolepton_charm_massdecoValue = (*DeepDiTau_boosted_nolepton_charm_massdecoValue)[jRef];
    if (DeepDiTau_boosted_nolepton_massdecoValid) DeepDiTau_boosted_nolepton_massdecoValue = (*DeepDiTau_boosted_nolepton_massdecoValue)[jRef];
    if (DeepDiTau_boosted_noleptonValid) DeepDiTau_boosted_noleptonValue = (*DeepDiTau_boosted_noleptonValue)[jRef];
    if (DeepDiTau_massdecoValid) DeepDiTau_massdecoValue = (*DeepDiTau_massdecoValue)[jRef];
    if (DeepDiTauValid) DeepDiTauValue = (*DeepDiTauValue)[jRef];
    if (DeepDiTau_nolepton_charm_massdecoValid) = (*DeepDiTau_nolepton_charm_massdecoValue)[jRef];
    if (DeepDiTau_nolepton_charmValid) = (*DeepDiTau_nolepton_charmValue)[jRef];
    if (DeepDiTau_nolepton_massdecoValid) = (*DeepDiTau_nolepton_massdecoValue)[jRef];
    if (DeepDiTau_noleptonValid) = (*DeepDiTau_noleptonValue)[jRef];

    jet.addUserFloat("ditau2017v1",ditau2017v1Value);
    jet.addUserFloat("ditau2017MDv1",ditau2017MDv1Value);
    jet.addUserFloat("DeepDiTau_boosted",DeepDiTau_boostedValue);
    jet.addUserFloat("DeepDiTau_boosted_massdeco",DeepDiTau_boosted_massdecoValue);
    jet.addUserFloat("DeepDiTau_boosted_nolepton_charm",DeepDiTau_boosted_nolepton_charmValue);
    jet.addUserFloat("DeepDiTau_boosted_nolepton_charm_massdeco",DeepDiTau_boosted_nolepton_charm_massdecoValue);
    jet.addUserFloat("DeepDiTau_boosted_nolepton_massdeco",DeepDiTau_boosted_nolepton_massdecoValue);
    jet.addUserFloat("DeepDiTau_boosted_nolepton",DeepDiTau_boosted_noleptonValue);
    jet.addUserFloat("DeepDiTau_massdeco",DeepDiTau_massdecoValue);
    jet.addUserFloat("DeepDiTau",DeepDiTauValue);
    jet.addUserFloat("DeepDiTau_nolepton_charm_massdeco",DeepDiTau_nolepton_charm_massdecoValue);
    jet.addUserFloat("DeepDiTau_nolepton_charm",DeepDiTau_nolepton_charmValue);
    jet.addUserFloat("DeepDiTau_nolepton_massdeco",DeepDiTau_nolepton_massdecoValue);
    jet.addUserFloat("DeepDiTau_nolepton",DeepDiTau_noleptonValue)

    output->push_back(jet);
  }

  evt.put(std::move(output));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetIdEmbedder);
