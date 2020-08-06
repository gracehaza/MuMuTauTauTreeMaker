// -*- C++ -*-
//
// Class:      DeepDiTauProducer
//
//
// Original Author:  Devin Taylor
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "MuMuTauTauTreeMaker/MuTauTreelizer/interface/DeepDiTau.h"
#include "MuMuTauTauTreeMaker/MuTauTreelizer/interface/DeepCache.h"
#include "DataFormats/Common/interface/ValueMap.h"

class DeepDiTauProducer : public edm::stream::EDProducer<edm::GlobalCache<DeepCache> > {
public:
  explicit DeepDiTauProducer(const edm::ParameterSet&, const DeepCache*);

  ~DeepDiTauProducer() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static std::unique_ptr<DeepCache> initializeGlobalCache(const edm::ParameterSet& cfg);
  static void globalEndJob(const DeepCache* cache) {}

private:
  edm::EDGetTokenT<std::vector<pat::Jet> > jetCollectionToken_;

  std::map<std::string, DeepDiTau*> deepDiTaus_;

  std::vector<std::string> deepDiTauLabels_;
  const DeepCache* cache_;

};

DeepDiTauProducer::DeepDiTauProducer(const edm::ParameterSet& iConfig, const DeepCache* cache)
    : cache_(cache) {
  auto inputLabel = iConfig.getParameter<edm::InputTag>("slimmedJetTag");

  // Load DeepDiTau
  auto deepCfg = iConfig.getParameter<edm::ParameterSet>("DeepDiTauConfiguration");
  auto graphDefinitions = deepCfg.getParameter<std::vector<edm::ParameterSet> >("graphDefinitions");
  for (auto graphDefinition : graphDefinitions) {
    std::string graphName = graphDefinition.getParameter<std::string>("name");
    deepDiTaus_[graphName] = new DeepDiTau(cache_);
    deepDiTaus_[graphName]->configure(graphDefinition);
    deepDiTauLabels_.push_back(graphName);
  }

  jetCollectionToken_ = consumes<std::vector<pat::Jet> >(inputLabel);

  // output will be association maps for deep id discriminants
  for (auto name : deepDiTauLabels_) {
    produces<edm::ValueMap<float> >(name);
  }
}

DeepDiTauProducer::~DeepDiTauProducer() {}

void DeepDiTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get jet collection
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByToken(jetCollectionToken_, jets);

  // create output map to store deep id discriminants
  std::map<std::string, std::vector<float> > deepDiTauScores;
  for (auto name : deepDiTauLabels_) {
    deepDiTauScores[name] = std::vector<float>();
  }

  // loop over jets and calculate the deep id discriminants
  for (const auto& jet : *jets) {
    for (auto pair : deepDiTaus_) {
      deepDiTauScores.at(pair.first).push_back(pair.second->evaluate(jet));
    }
  }


  // store them all in the event as value maps
  for (auto name : deepDiTauLabels_) {
    std::unique_ptr<edm::ValueMap<float> > output(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler output_filler(*output);
    output_filler.insert(jets, deepDiTauScores.at(name).begin(), deepDiTauScores.at(name).end());
    output_filler.fill();
    iEvent.put(std::move(output),name);
  }
}
std::unique_ptr<DeepCache> DeepDiTauProducer::initializeGlobalCache(const edm::ParameterSet& cfg) {
  std::map<std::string, std::string> graphNames;
  if (!cfg.exists("DeepDiTauConfiguration")) {
    return std::make_unique<DeepCache>(graphNames, false);
  }
  edm::ParameterSet deepCfg = cfg.getParameter<edm::ParameterSet>("DeepDiTauConfiguration");
  auto graphDefinitions = deepCfg.getParameter<std::vector<edm::ParameterSet> >("graphDefinitions");
  bool memmapped = deepCfg.getParameter<bool>("memmapped");
  for (auto graphDefinition : graphDefinitions) {
    std::string graphName = graphDefinition.getParameter<std::string>("name");
    edm::FileInPath graphPath = graphDefinition.getParameter<edm::FileInPath>("path");
    std::string graphFullPath = graphPath.fullPath();
    graphNames[graphName] = graphFullPath;
  }
  return std::make_unique<DeepCache>(graphNames, memmapped);
}

void DeepDiTauProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setAllowAnything();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DeepDiTauProducer);
