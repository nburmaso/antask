#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include <SimulationDataFormat/MCTrack.h>
#include <TFile.h>
#include <TTree.h>
#include "TH1D.h"
#include <vector>
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <SimulationDataFormat/MCCompLabel.h>
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "TOFSimulation/Detector.h"
#include "DataFormatsTOF/Cluster.h"
#include "GlobalTracking/MatchTOF.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;

class ATask : public AnalysisTask
{
public:
  OutputObj<TH2F> hPhi{TH1F("phi", "Phi", 100, 0., 2. * M_PI, 102, -2.01, 2.01)};
  
  TH1D* hMass = new TH1D("hMass", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hMassReco = new TH1D("hMassReco", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hMassReco1 = new TH1D("hMassReco1", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hMassReco2 = new TH1D("hMassReco2", "; mass, GeV; N", 50, 0.5, 4);

  void process(aod:Tracks const& tracks)
  {
    for (auto& track : tracks) {
      float phi = asin(track.snp()) + track.alpha() + M_PI;
      hPhi-Fill(phi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
 return WorkflowSpec{
 adaptAnalysisTask<ATask>("mySimpleTrackAnalysis", 0)
 };
}