/*
 * \class DeepDiTau
 *
 * DiTau identification using DNN.
 *
 * \author Devin Taylor, UC Davis
 */

#include "MuMuTauTauTreeMaker/MuTauTreelizer/interface/DeepDiTau.h"
#include <fstream>
#include <sstream>

// Define the inputs to the DNN
namespace {

  namespace empty {
    constexpr int NumberOfOutputs = 3; // Defines the number of elements in the output vector
  }

  namespace ditauInputs_2017_v1 {
    constexpr int NumberOfOutputs = 3;
    constexpr int NumberOfChargedHadrons = 10;
    constexpr int NumberOfNeutralHadrons = 10;
    constexpr int NumberOfMuons = 4;
    constexpr int NumberOfElectrons = 4;
    constexpr int NumberOfPhotons = 4;
    namespace JetBlockInputs {
      enum vars {
        jet_pt = 0,
        jet_eta,
        jet_phi,
        jet_mass,
        jet_jetCharge,
        jet_chargedMultiplicity,
        jet_neutralMultiplicity,
        jet_chargedHadronMultiplicity,
        jet_neutralHadronMultiplicity,
        jet_muonMultiplicity,
        jet_electronMultiplicity,
        jet_photonMultiplicity,
        jet_chargedEmEnergy,
        jet_neutralEmEnergy,
        jet_chargedHadronEnergy,
        jet_neutralHadronEnergy,
        jet_muonEnergy,
        jet_electronEnergy,
        jet_photonEnergy,
        jet_chargedEmEnergyFraction,
        jet_neutralEmEnergyFraction,
        jet_chargedHadronEnergyFraction,
        jet_neutralHadronEnergyFraction,
        jet_muonEnergyFraction,
        jet_electronEnergyFraction,
        jet_photonEnergyFraction,
        jet_pfJetBProbabilityBJetTags,
        jet_pfJetProbabilityBJetTags,
        jet_pfTrackCountingHighEffBJetTags,
        jet_pfSimpleSecondaryVertexHighEffBJetTags,
        jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags,
        jet_pfCombinedSecondaryVertexV2BJetTags,
        jet_pfCombinedInclusiveSecondaryVertexV2BJetTags,
        jet_pfCombinedMVAV2BJetTags,
        jet_pfCombinedCvsLJetTags,
        jet_pfCombinedCvsBJetTags,
        jet_pfDeepCSVJetTags_probb,
        jet_pfDeepCSVJetTags_probc,
        jet_pfDeepCSVJetTags_probudsg,
        jet_pfDeepCSVJetTags_probbb,
        NumberOfInputs
      };
    }

    namespace ChargedHadronBlockInputs {
      enum vars {
        charged_hadron_pt = 0,
        charged_hadron_eta,
        charged_hadron_phi,
        charged_hadron_charge,
        charged_hadron_etaAtVtx,
        charged_hadron_phiAtVtx,
        charged_hadron_vx,
        charged_hadron_vy,
        charged_hadron_vz,
        charged_hadron_dxy,
        charged_hadron_dz,
        charged_hadron_pixelLayersWithMeasurement,
        charged_hadron_stripLayersWithMeasurement,
        charged_hadron_trackerLayersWithMeasurement,
        charged_hadron_trackHighPurity,
        charged_hadron_puppiWeight,
        charged_hadron_puppiWeightNoLep,
        charged_hadron_isIsolatedChargedHadron,
        NumberOfInputs
      };
    }

    namespace NeutralHadronBlockInputs {
      enum vars {
        neutral_hadron_pt = 0,
        neutral_hadron_eta,
        neutral_hadron_phi,
        neutral_hadron_puppiWeight,
        neutral_hadron_puppiWeightNoLep,
        NumberOfInputs
      };
    }

    namespace MuonBlockInputs {
      enum vars {
        muon_pt = 0,
        muon_eta,
        muon_phi,
        muon_charge,
        muon_etaAtVtx,
        muon_phiAtVtx,
        muon_vx,
        muon_vy,
        muon_vz,
        muon_dxy,
        muon_dz,
        muon_pixelLayersWithMeasurement,
        muon_stripLayersWithMeasurement,
        muon_trackerLayersWithMeasurement,
        muon_trackHighPurity,
        muon_puppiWeight,
        muon_isStandAloneMuon,
        muon_isGlobalMuon,
        NumberOfInputs
      };
    }

    namespace ElectronBlockInputs {
      enum vars {
        electron_pt = 0,
        electron_eta,
        electron_phi,
        electron_charge,
        electron_etaAtVtx,
        electron_phiAtVtx,
        electron_vx,
        electron_vy,
        electron_vz,
        electron_dxy,
        electron_dz,
        electron_pixelLayersWithMeasurement,
        electron_stripLayersWithMeasurement,
        electron_trackerLayersWithMeasurement,
        electron_trackHighPurity,
        electron_puppiWeight,
        NumberOfInputs
      };
    }

    namespace PhotonBlockInputs {
      enum vars {
        photon_pt = 0,
        photon_eta,
        photon_phi,
        photon_puppiWeight,
        photon_puppiWeightNoLep,
        photon_isGoodEgamma,
        NumberOfInputs
      };
    }

  }  // namespace ditauInputs_2017_v1

}  // anonymous namespace

  
void DeepDiTau::configure(const edm::ParameterSet& iConfig)
{

    if (isConfigured_) return;

    name_ = iConfig.getParameter<std::string>("name");
    path_ = iConfig.getParameter<edm::FileInPath>("path");
    meanPath_ = iConfig.getParameter<edm::FileInPath>("means");

    // load the means and sigmas
    std::string fullMeansPath = meanPath_.fullPath();
    std::ifstream inputMeansFile(fullMeansPath);
    std::string line, word, name;
    float mean, sigma;
    while (std::getline(inputMeansFile, line)) {
      std::istringstream ss(line);
      std::getline(ss,word,',');
      name = word;
      std::getline(ss,word,',');
      mean = std::stof(word);
      std::getline(ss,word,',');
      sigma = std::stof(word);
      names_.push_back(name);
      means_.push_back(mean);
      sigmas_.push_back(sigma);
    }
    inputMeansFile.close();

    // define graphs
    if (name_=="ditau2017v1") {
      inputNames_.push_back("input_1");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::JetBlockInputs::NumberOfInputs});
      kJet_ = 0;
      inputNames_.push_back("input_2");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::ChargedHadronBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfChargedHadrons}); 
      kChargedHadron_ = 1;
      inputNames_.push_back("input_3");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::NeutralHadronBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfNeutralHadrons}); 
      kNeutralHadron_ = 2;
      inputNames_.push_back("input_4");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::MuonBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfMuons}); 
      kMuon_ = 3;
      inputNames_.push_back("input_5");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::ElectronBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfElectrons}); 
      kElectron_ = 4;
      inputNames_.push_back("input_6");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::PhotonBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfPhotons}); 
      kPhoton_ = 5;
      outputName_ = "ID_pred/Softmax";
    }
    else if (name_=="ditau2017MDv1") {
      inputNames_.push_back("input_1");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::JetBlockInputs::NumberOfInputs});
      kJet_ = 0;
      inputNames_.push_back("input_2");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::ChargedHadronBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfChargedHadrons}); 
      kChargedHadron_ = 1;
      inputNames_.push_back("input_3");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::NeutralHadronBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfNeutralHadrons}); 
      kNeutralHadron_ = 2;
      inputNames_.push_back("input_4");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::MuonBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfMuons}); 
      kMuon_ = 3;
      inputNames_.push_back("input_5");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::ElectronBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfElectrons}); 
      kElectron_ = 4;
      inputNames_.push_back("input_6");
      inputShapes_.push_back(tensorflow::TensorShape{1, ditauInputs_2017_v1::PhotonBlockInputs::NumberOfInputs, ditauInputs_2017_v1::NumberOfPhotons}); 
      kPhoton_ = 5;
      //inputNames_.push_back("input_7");
      //inputShapes_.push_back(tensorflow::TensorShape{1, 1});
      //kMass_ = 6;
      //outputName_ = "concatenate_2/concat"; // TODO lookup, and probably rerun and give a more consistent name
      outputName_ = "ID_pred/Softmax";
    }

    inputTensors_.resize(inputShapes_.size());
    for (size_t i=0; i<inputShapes_.size(); i++) {
      inputTensors_[i] = tensorflow::NamedTensor(inputNames_[i], tensorflow::Tensor(tensorflow::DT_FLOAT, inputShapes_.at(i)));
    }

    // now validate the graph
    auto graph = cache_->getGraph(name_);
    for (size_t i=0; i<inputShapes_.size(); i++){
      const auto& name = graph.node(i).name();
      // not necessary to be in same order in the input graph
      auto it = std::find(inputNames_.begin(), inputNames_.end(), name);
      if (it==inputNames_.end()) {
        throw cms::Exception("DeepDiTau")
          << "Processing graph " << name_ << ".\n"
          << "Unknown input name " << name;
      }
      const auto& shape = graph.node(i).attr().at("shape").shape();
      int j = std::distance(inputNames_.begin(),it);
      for (int d=1; d<inputShapes_.at(j).dims(); d++) { // skip first dim since it should be -1 and not 1 like we define here for evaluation
        if (shape.dim(d).size() != inputShapes_.at(j).dim_size(d)) {
          throw cms::Exception("DeepDiTau")
            << "Number of inputs in graph does not match those expected for " << name_ << ".\n"
            << "Expected input " << j << " dim " << d << " = " << inputShapes_.at(j).dim_size(d) << "."
            << " Found " << shape.dim(d).size() << ".";
        }
      }
    }
    const auto& outName = graph.node(graph.node_size() - 1).name();
    if (outName!=outputName_) {
      throw cms::Exception("DeepDiTau")
        << "Processing graph " << name_ << ".\n"
        << "Unexpected output name. Expected " << outputName_ << " found " << name << ".";
    }

    isConfigured_ = true;
}

float DeepDiTau::evaluate(const pat::Jet& jet) {
    const tensorflow::Tensor pred = getPrediction(jet);
    float ditau_score = pred.matrix<float>()(0,2); // hard coded for now
    return ditau_score;
}


// get the prediction for the selected DNN
tensorflow::Tensor DeepDiTau::getPrediction(const pat::Jet& jet) {

  std::vector<tensorflow::Tensor> pred_vector;
  tensorflow::Tensor prediction;

  if (name_=="ditau2017v1") {
    getPrediction_2017_v1(jet, pred_vector);
    prediction = tensorflow::Tensor(tensorflow::DT_FLOAT, {1, ditauInputs_2017_v1::NumberOfOutputs});
    for (int k = 0; k < ditauInputs_2017_v1::NumberOfOutputs; ++k) {
      const float pred = pred_vector[0].flat<float>()(k); // just one prediction vector for now
      if (!(pred >= 0 && pred <= 1)) {
        throw cms::Exception("DeepDiTau")
            << "invalid prediction = " << pred << " for pred_index = " << k;
      } 
      prediction.matrix<float>()(0, k) = pred;
    }
  } // ditau2017v1
  else if (name_=="ditau2017MDv1") {
    getPrediction_2017_md_v1(jet, pred_vector);
    prediction = tensorflow::Tensor(tensorflow::DT_FLOAT, {1, ditauInputs_2017_v1::NumberOfOutputs+1});
    for (int k = 0; k < ditauInputs_2017_v1::NumberOfOutputs; ++k) {
      const float pred = pred_vector[0].flat<float>()(k); // just one prediction vector for now
      if (!(pred >= 0 && pred <= 1)) {
        throw cms::Exception("DeepDiTau")
            << "invalid prediction = " << pred << " for pred_index = " << k;
      } 
      prediction.matrix<float>()(0, k) = pred;
    }
  } // ditau2017MDv1
  else {
    prediction = tensorflow::Tensor(tensorflow::DT_FLOAT, {1, empty::NumberOfOutputs});
    prediction.matrix<float>().setZero();
    prediction.matrix<float>()(0, 2) = -1.0; // hard coded for now
  }
  
  return prediction;
}

// ditau2017v1
void DeepDiTau::getPrediction_2017_v1(const pat::Jet& jet, std::vector<tensorflow::Tensor>& pred_vector) {
  createJetBlockInputs(jet);
  createChargedHadronBlockInputs(jet);
  createNeutralHadronBlockInputs(jet);
  createMuonBlockInputs(jet);
  createElectronBlockInputs(jet);
  createPhotonBlockInputs(jet);

  //std::cout << "input jet " << inputTensors_.at(kJet_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kJet_).second.matrix<float>() << std::endl;
  //std::cout << "input charged hadron " << inputTensors_.at(kChargedHadron_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kChargedHadron_).second.tensor<float,3>() << std::endl;
  //std::cout << "input neutral hadron " << inputTensors_.at(kNeutralHadron_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kNeutralHadron_).second.tensor<float,3>() << std::endl;
  //std::cout << "input muon " << inputTensors_.at(kMuon_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kMuon_).second.tensor<float,3>() << std::endl;
  //std::cout << "input electron " << inputTensors_.at(kElectron_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kElectron_).second.tensor<float,3>() << std::endl;
  //std::cout << "input photon " << inputTensors_.at(kPhoton_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kPhoton_).second.tensor<float,3>() << std::endl;

  tensorflow::run(&(cache_->getSession(name_)),
                  inputTensors_,
                  {outputName_},
                  &pred_vector);

  //std::cout << "prediction " << pred_vector[0].DebugString() << std::endl;
  //std::cout << pred_vector[0].matrix<float>() << std::endl;

}

// ditau2017MDv1
void DeepDiTau::getPrediction_2017_md_v1(const pat::Jet& jet, std::vector<tensorflow::Tensor>& pred_vector) {
  createJetBlockInputs(jet);
  createChargedHadronBlockInputs(jet);
  createNeutralHadronBlockInputs(jet);
  createMuonBlockInputs(jet);
  createElectronBlockInputs(jet);
  createPhotonBlockInputs(jet);
  //createMassInput(jet);

  //std::cout << "input jet " << inputTensors_.at(kJet_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kJet_).second.matrix<float>() << std::endl;
  //std::cout << "input charged hadron " << inputTensors_.at(kChargedHadron_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kChargedHadron_).second.tensor<float,3>() << std::endl;
  //std::cout << "input neutral hadron " << inputTensors_.at(kNeutralHadron_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kNeutralHadron_).second.tensor<float,3>() << std::endl;
  //std::cout << "input muon " << inputTensors_.at(kMuon_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kMuon_).second.tensor<float,3>() << std::endl;
  //std::cout << "input electron " << inputTensors_.at(kElectron_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kElectron_).second.tensor<float,3>() << std::endl;
  //std::cout << "input photon " << inputTensors_.at(kPhoton_).second.DebugString() << std::endl;
  //std::cout << inputTensors_.at(kPhoton_).second.tensor<float,3>() << std::endl;

  tensorflow::run(&(cache_->getSession(name_)),
                  inputTensors_,
                  {outputName_},
                  &pred_vector);

  //std::cout << "prediction " << pred_vector[0].DebugString() << std::endl;
  //std::cout << pred_vector[0].matrix<float>() << std::endl;

}

void DeepDiTau::createMassInput(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kMass_).second;
    inputs.flat<float>().setZero();

    namespace dnn = ditauInputs_2017_v1::JetBlockInputs;
    int v;
    v = dnn::jet_mass;
    inputs.matrix<float>()(0, 0) = getValueNorm(jet.mass(), means_[v], sigmas_[v]);
}

void DeepDiTau::createJetBlockInputs(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kJet_).second;
    inputs.flat<float>().setZero();

    namespace dnn = ditauInputs_2017_v1::JetBlockInputs;
    int v;
    v = dnn::jet_pt;
    inputs.matrix<float>()(0, v) = getValueLogLinear(jet.pt(), means_[v], sigmas_[v]);
    v = dnn::jet_eta;
    inputs.matrix<float>()(0, v) = getValueLinear(jet.eta(), means_[v], sigmas_[v]);
    v = dnn::jet_phi;
    inputs.matrix<float>()(0, v) = getValueLinear(jet.phi(), means_[v], sigmas_[v]);
    v = dnn::jet_mass;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.mass(), means_[v], sigmas_[v]);
    v = dnn::jet_jetCharge;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.jetCharge(), means_[v], sigmas_[v]);
    v = dnn::jet_chargedMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.chargedMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_neutralMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.neutralMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_chargedHadronMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.chargedHadronMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_neutralHadronMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.neutralHadronMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_muonMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.muonMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_electronMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.electronMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_photonMultiplicity;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.photonMultiplicity(), means_[v], sigmas_[v]);
    v = dnn::jet_chargedEmEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.chargedEmEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_neutralEmEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.neutralEmEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_chargedHadronEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.chargedHadronEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_neutralHadronEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.neutralHadronEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_muonEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.muonEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_electronEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.electronEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_photonEnergy;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.photonEnergy(), means_[v], sigmas_[v]);
    v = dnn::jet_chargedEmEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.chargedEmEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_neutralEmEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.neutralEmEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_chargedHadronEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.chargedHadronEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_neutralHadronEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.neutralHadronEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_muonEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.muonEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_electronEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.electronEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_photonEnergyFraction;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.photonEnergyFraction(), means_[v], sigmas_[v]);
    v = dnn::jet_pfJetBProbabilityBJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfJetBProbabilityBJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfJetProbabilityBJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfJetProbabilityBJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfTrackCountingHighEffBJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfTrackCountingHighEffBJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfSimpleSecondaryVertexHighEffBJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfCombinedSecondaryVertexV2BJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfCombinedInclusiveSecondaryVertexV2BJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfCombinedMVAV2BJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfCombinedMVAV2BJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfCombinedCvsLJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfCombinedCvsLJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfCombinedCvsBJetTags;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfCombinedCvsBJetTags"), means_[v], sigmas_[v]);
    v = dnn::jet_pfDeepCSVJetTags_probb;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfDeepCSVJetTags:probb"), means_[v], sigmas_[v]);
    v = dnn::jet_pfDeepCSVJetTags_probc;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfDeepCSVJetTags:probc"), means_[v], sigmas_[v]);
    v = dnn::jet_pfDeepCSVJetTags_probudsg;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfDeepCSVJetTags:probudsg"), means_[v], sigmas_[v]);
    v = dnn::jet_pfDeepCSVJetTags_probbb;
    inputs.matrix<float>()(0, v) = getValueNorm(jet.bDiscriminator("pfDeepCSVJetTags:probbb"), means_[v], sigmas_[v]);
}

void DeepDiTau::createChargedHadronBlockInputs(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kChargedHadron_).second;
    inputs.flat<float>().setZero();

    namespace dnnj = ditauInputs_2017_v1::JetBlockInputs;
    namespace dnn = ditauInputs_2017_v1::ChargedHadronBlockInputs;
    // in the means list, the charged hadron are after the jet
    int n = dnnj::NumberOfInputs;
    int v;
    int idh = 0;
    for (size_t d=0; d<jet.numberOfDaughters(); d++) {
      // limit the number of daughters
      if (idh>=ditauInputs_2017_v1::NumberOfChargedHadrons) continue;
      const pat::PackedCandidate * dau = (pat::PackedCandidate*)jet.daughter(d);
      // require charged hadron
      if ((std::abs(dau->pdgId())==13) || (std::abs(dau->pdgId())==11) || (std::abs(dau->pdgId())==22) || (std::abs(dau->charge())==0)) continue;
      v = dnn::charged_hadron_pt;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pt()/jet.pt(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_eta;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->eta()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_phi;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phi()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_charge;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->charge(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_etaAtVtx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->etaAtVtx()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_phiAtVtx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phiAtVtx()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_vx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vx(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_vy;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vy(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_vz;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vz(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_dxy;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->dxy(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_dz;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->dz(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_pixelLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pixelLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_stripLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->stripLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_trackerLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->trackerLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_trackHighPurity;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->trackHighPurity(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_puppiWeight;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeight(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_puppiWeightNoLep;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeightNoLep(), means_[v+n], sigmas_[v+n]);
      v = dnn::charged_hadron_isIsolatedChargedHadron;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->isIsolatedChargedHadron(), means_[v+n], sigmas_[v+n]);
      idh++;
    }

}

void DeepDiTau::createNeutralHadronBlockInputs(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kNeutralHadron_).second;
    inputs.flat<float>().setZero();

    namespace dnnj = ditauInputs_2017_v1::JetBlockInputs;
    namespace dnnch = ditauInputs_2017_v1::ChargedHadronBlockInputs;
    namespace dnn = ditauInputs_2017_v1::NeutralHadronBlockInputs;
    // in the means list, the neutral hadron are after the jet and charged hadron
    int n = dnnj::NumberOfInputs + dnnch::NumberOfInputs;
    int v;
    int idh = 0;
    for (size_t d=0; d<jet.numberOfDaughters(); d++) {
      // limit the number of daughters
      if (idh>=ditauInputs_2017_v1::NumberOfNeutralHadrons) continue;
      const pat::PackedCandidate * dau = (pat::PackedCandidate*)jet.daughter(d);
      // require neutral hadron
      if ((std::abs(dau->pdgId())==13) || (std::abs(dau->pdgId())==11) || (std::abs(dau->pdgId())==22) || (std::abs(dau->charge())>0)) continue;
      v = dnn::neutral_hadron_pt;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pt()/jet.pt(), means_[v+n], sigmas_[v+n]);
      v = dnn::neutral_hadron_eta;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->eta()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::neutral_hadron_phi;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phi()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::neutral_hadron_puppiWeight;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeight(), means_[v+n], sigmas_[v+n]);
      v = dnn::neutral_hadron_puppiWeightNoLep;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeightNoLep(), means_[v+n], sigmas_[v+n]);
      idh++;
    }

}


void DeepDiTau::createMuonBlockInputs(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kMuon_).second;
    inputs.flat<float>().setZero();

    namespace dnnj = ditauInputs_2017_v1::JetBlockInputs;
    namespace dnnch = ditauInputs_2017_v1::ChargedHadronBlockInputs;
    namespace dnnnh = ditauInputs_2017_v1::NeutralHadronBlockInputs;
    namespace dnn = ditauInputs_2017_v1::MuonBlockInputs;
    // in the means list, the muon are after the jet, charged hadron, neutral hadron
    int n = dnnj::NumberOfInputs + dnnch::NumberOfInputs + dnnnh::NumberOfInputs;
    int v;
    int idh = 0;
    for (size_t d=0; d<jet.numberOfDaughters(); d++) {
      // limit the number of daughters
      if (idh>=ditauInputs_2017_v1::NumberOfMuons) continue;
      const pat::PackedCandidate * dau = (pat::PackedCandidate*)jet.daughter(d);
      // require muon
      if (std::abs(dau->pdgId())!=13) continue;
      v = dnn::muon_pt;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pt()/jet.pt(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_eta;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->eta()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_phi;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phi()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_charge;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->charge(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_etaAtVtx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->etaAtVtx()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_phiAtVtx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phiAtVtx()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_vx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vx(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_vy;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vy(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_vz;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vz(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_dxy;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->dxy(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_dz;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->dz(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_pixelLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pixelLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_stripLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->stripLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_trackerLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->trackerLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_trackHighPurity;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->trackHighPurity(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_puppiWeight;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeight(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_isStandAloneMuon;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->isStandAloneMuon(), means_[v+n], sigmas_[v+n]);
      v = dnn::muon_isGlobalMuon;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->isGlobalMuon(), means_[v+n], sigmas_[v+n]);
      idh++;
    }

}

void DeepDiTau::createElectronBlockInputs(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kElectron_).second;
    inputs.flat<float>().setZero();

    namespace dnnj = ditauInputs_2017_v1::JetBlockInputs;
    namespace dnnch = ditauInputs_2017_v1::ChargedHadronBlockInputs;
    namespace dnnnh = ditauInputs_2017_v1::NeutralHadronBlockInputs;
    namespace dnnm = ditauInputs_2017_v1::MuonBlockInputs;
    namespace dnn = ditauInputs_2017_v1::ElectronBlockInputs;
    // in the means list, the electron are after the jet, charged hadron, neutral hadron, muon
    int n = dnnj::NumberOfInputs + dnnch::NumberOfInputs + dnnnh::NumberOfInputs + dnnm::NumberOfInputs;
    int v;
    int idh = 0;
    for (size_t d=0; d<jet.numberOfDaughters(); d++) {
      // limit the number of daughters
      if (idh>=ditauInputs_2017_v1::NumberOfElectrons) continue;
      const pat::PackedCandidate * dau = (pat::PackedCandidate*)jet.daughter(d);
      // require electron
      if (std::abs(dau->pdgId())!=11) continue;
      v = dnn::electron_pt;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pt()/jet.pt(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_eta;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->eta()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_phi;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phi()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_charge;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->charge(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_etaAtVtx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->etaAtVtx()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_phiAtVtx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phiAtVtx()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_vx;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vx(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_vy;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vy(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_vz;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->vz(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_dxy;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->dxy(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_dz;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->dz(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_pixelLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pixelLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_stripLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->stripLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_trackerLayersWithMeasurement;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->trackerLayersWithMeasurement(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_trackHighPurity;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->trackHighPurity(), means_[v+n], sigmas_[v+n]);
      v = dnn::electron_puppiWeight;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeight(), means_[v+n], sigmas_[v+n]);
      idh++;
    }

}

void DeepDiTau::createPhotonBlockInputs(const pat::Jet& jet) {
    tensorflow::Tensor& inputs = inputTensors_.at(kPhoton_).second;
    inputs.flat<float>().setZero();

    namespace dnnj = ditauInputs_2017_v1::JetBlockInputs;
    namespace dnnch = ditauInputs_2017_v1::ChargedHadronBlockInputs;
    namespace dnnnh = ditauInputs_2017_v1::NeutralHadronBlockInputs;
    namespace dnnm = ditauInputs_2017_v1::MuonBlockInputs;
    namespace dnne = ditauInputs_2017_v1::ElectronBlockInputs;
    namespace dnn = ditauInputs_2017_v1::PhotonBlockInputs;
    // in the means list, the photon are after the jet, charged hadron, neutral hadron, muon, photon
    int n = dnnj::NumberOfInputs + dnnch::NumberOfInputs + dnnnh::NumberOfInputs + dnnm::NumberOfInputs + dnne::NumberOfInputs;
    int v;
    int idh = 0;
    for (size_t d=0; d<jet.numberOfDaughters(); d++) {
      // limit the number of daughters
      if (idh>=ditauInputs_2017_v1::NumberOfPhotons) continue;
      const pat::PackedCandidate * dau = (pat::PackedCandidate*)jet.daughter(d);
      // require photon
      if (std::abs(dau->pdgId())!=22) continue;
      v = dnn::photon_pt;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->pt()/jet.pt(), means_[v+n], sigmas_[v+n]);
      v = dnn::photon_eta;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->eta()-jet.eta(), means_[v+n], sigmas_[v+n]);
      v = dnn::photon_phi;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->phi()-jet.phi(), means_[v+n], sigmas_[v+n]);
      v = dnn::photon_puppiWeight;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeight(), means_[v+n], sigmas_[v+n]);
      v = dnn::photon_puppiWeightNoLep;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->puppiWeightNoLep(), means_[v+n], sigmas_[v+n]);
      v = dnn::photon_isGoodEgamma;
      inputs.tensor<float,3>()(0, v, idh) = getValueNorm(dau->isGoodEgamma(), means_[v+n], sigmas_[v+n]);
      idh++;
    }

}

