#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include <TH1F.h>

using namespace o2;
using namespace o2::framework;

class ATask : public AnalysisTask {

public:
  OutputObj<TH2F> hPhi{TH1F("phi", "Phi", 100, 0., 2. * M_PI, 102, -2.01, 2.01)};

  void process(aod:Tracks const& tracks) {
    for (auto& track : tracks) {
      float phi = asin(track.snp()) + track.alpha() + M_PI;
      hPhi-Fill(phi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&) {
  return WorkflowSpec {
    adaptAnalysisTask<ATask>("mySimpleTrackAnalysis", 0)
  };
}