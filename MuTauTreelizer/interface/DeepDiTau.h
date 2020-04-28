//#ifndef DevTools_Ntuplizer_DeepDiTau_h
//#define DevTools_Ntuplizer_DeepDiTau_h
#ifndef MuMuTauTauTreeMaker_MuTauTreelizer_DeepDiTau_h
#define MuMuTauTauTreeMaker_MuTauTreelizer_DeepDiTau_h

/*
 * \class DeepDiTau
 *
 * Base class for ditau identification using TensorFlow DNN.
 *
 * \author Devin Taylor, UC Davis
 */

#include <Math/VectorUtil.h>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "MuMuTauTauTreeMaker/MuTauTreelizer/interface/DeepCache.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

class DeepDiTau {
  public:
    DeepDiTau(const DeepCache* cache)
        : isConfigured_(false),
        cache_(cache)
        {}
    void configure(const edm::ParameterSet&);
    float evaluate( const pat::Jet& );
    ~DeepDiTau() {}

  private:
    static constexpr float default_value = -999.;
    // Utility to convert value to a normalized float input for the DNN
    template <typename T>
    static float getValue(T value) {
      return std::isnormal(value) ? static_cast<float>(value) : 0.f;
    }

    template <typename T>
    static float getValueLinear(T value, float min_value, float max_value) {
      const float fixed_value = getValue(value);
      const float clamped_value = std::clamp(fixed_value, min_value, max_value);
      float transformed_value = (clamped_value - min_value) / (max_value - min_value);
      //std::cout << "linearized " << value << " " << min_value << " " << max_value << " " << transformed_value << std::endl;
      return transformed_value;
    }

    template <typename T>
    static float getValueLogLinear(T value, float min_value, float max_value) {
      const float fixed_value = getValue(value);
      const float clamped_value = std::clamp(fixed_value, min_value, max_value);
      float transformed_value = (log(clamped_value) - log(min_value)) / (log(max_value) - log(min_value));
      //std::cout << "loglinearized " << value << " " << min_value << " " << max_value << " " << transformed_value << std::endl;
      return transformed_value;
    }

    template <typename T>
    static float getValueNorm(T value, float mean, float sigma, float n_sigmas_max = -1) {
      const float fixed_value = getValue(value);
      const float norm_value = (fixed_value - mean) / sigma;
      float result;
      if (n_sigmas_max>0) {
          result = std::clamp(norm_value, -n_sigmas_max, n_sigmas_max);
      }
      else {
          result = norm_value;
      }
      //std::cout << "normalized " << value << " " << mean << " " << sigma << " " << result << std::endl;
      return result;
    }

    void createJetBlockInputs(const pat::Jet&);
    void createChargedHadronBlockInputs(const pat::Jet&);
    void createNeutralHadronBlockInputs(const pat::Jet&);
    void createMuonBlockInputs(const pat::Jet&);
    void createElectronBlockInputs(const pat::Jet&);
    void createPhotonBlockInputs(const pat::Jet&);
    void createMassInput(const pat::Jet&);
    tensorflow::Tensor getPrediction(const pat::Jet&);
    void getPrediction_2017_v1(const pat::Jet&, std::vector<tensorflow::Tensor>&);
    void getPrediction_2017_md_v1(const pat::Jet&, std::vector<tensorflow::Tensor>&);

  protected:
    bool isConfigured_;
    const DeepCache* cache_;
  private:
    std::string name_;
    edm::FileInPath path_;
    edm::FileInPath meanPath_;
    std::vector<std::string> names_;
    std::vector<float> means_;
    std::vector<float> sigmas_;
    std::vector<tensorflow::TensorShape> inputShapes_;
    std::vector<std::string> inputNames_;
    tensorflow::NamedTensorList inputTensors_;
    size_t kJet_;
    size_t kChargedHadron_;
    size_t kNeutralHadron_;
    size_t kMuon_;
    size_t kElectron_;
    size_t kPhoton_;
    size_t kMass_;
    std::string outputName_;
};

#endif
