/*#ifndef DevTools_Ntuplizer_DeepCache_h
#define DevTools_Ntuplizer_DeepCache_h
*/
#ifndef MuMuTauTauTreeMaker_MuTauTreelizer_DeepCache_h
#define MuMuTauTauTreeMaker_MuTauTreelizer_DeepCache_h
/*
 * \class DeepCache
 *
 * Cache for tensorflow networks
 *
 * \author Devin Taylor, UC Davis
 */

#include <Math/VectorUtil.h>
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "tensorflow/core/util/memmapped_file_system.h"

class DeepCache {
  public:
    DeepCache(const std::map<std::string, std::string>& graph_names, bool mem_mapped);
    ~DeepCache();

    tensorflow::Session& getSession(const std::string& name = "") const { return *sessions_.at(name); }
    const tensorflow::GraphDef& getGraph(const std::string& name = "") const { return *graphs_.at(name); }

  private:
    std::map<std::string, std::shared_ptr<tensorflow::GraphDef> > graphs_;
    std::map<std::string, tensorflow::Session*> sessions_;
    std::map<std::string, std::unique_ptr<tensorflow::MemmappedEnv> > memmappedEnv_;
};
#endif
