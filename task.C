#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "RootTreeReader.h"
#include "RootTreeWriter.h"

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

#define PI 3.14159265

class ATask : public AnalysisTask
{
public:
  void process(aod::Tracks const& tracks)
  {
    // mc
    TFile f("o2sim.root");
    TTree* tree = (TTree*) f.Get("o2sim");
    vector<MCTrack>* mctracks = nullptr;
    tree->SetBranchAddress("MCTrack",&mctracks);

    // tof clusters
    TFile fTof("tofclusters.root");
    TTree* treeTof = (TTree*) fTof.Get("o2sim");
    std::vector<o2::tof::Cluster>* clusters = nullptr;
    treeTof->SetBranchAddress("TOFCluster",&clusters);

    // tof match
    TFile tofMatch("o2match_tof.root");
    TTree* treeTofMatch = (TTree*) tofMatch.Get("matchTOF");
    std::vector<o2::dataformats::MatchInfoTOF>* TOFMatchInfo = new std::vector<o2::dataformats::MatchInfoTOF>;
    treeTofMatch->SetBranchAddress("TOFMatchInfo", &TOFMatchInfo);

    //tpc-its
    TFile* fmatchITSTPC = new TFile("o2match_itstpc.root");
    TTree* tracksTree = (TTree*)fmatchITSTPC->Get("matchTPCITS");
    std::vector<o2::dataformats::TrackTPCITS>* tpcItsTracks = new std::vector<o2::dataformats::TrackTPCITS>;
    tracksTree->SetBranchAddress("TPCITS", &tpcItsTracks);
    std::vector<o2::MCCompLabel>* mcMatchTPC = new std::vector<o2::MCCompLabel>;
    tracksTree->SetBranchAddress("MatchTPCMCTruth", &mcMatchTPC);

    TH1D* hMass = new TH1D("hMass", "; mass, GeV; N", 50, 0.5, 4);
    TH1D* hMassReco = new TH1D("hMassReco", "; mass, GeV; N", 50, 0.5, 4);
    TH1D* hMassReco1 = new TH1D("hMassReco1", "; mass, GeV; N", 50, 0.5, 4);
    TH1D* hMassReco2 = new TH1D("hMassReco2", "; mass, GeV; N", 50, 0.5, 4);

    hMass->Sumw2();
    hMassReco->Sumw2();
    hMassReco1->Sumw2();
    hMassReco2->Sumw2();

    TDatabasePDG *pdg = new TDatabasePDG();
    std::vector<int> mcID;
    std::vector<int> recoInd;
    int countRec = 0;
    int countTof = 0;

    Double_t px = 0.0;
    Double_t py = 0.0;
    Double_t pz = 0.0;
    Double_t m = 0.0;
    Double_t invM = 0.0;
    TLorentzVector pl1, pl2, pl;
    TVector3 p;

  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    EffAnlysis<ATask>("EffAnalysis")
  };
}

// #include "Framework/runDataProcessing.h"
// #include "Framework/AnalysisTask.h"
// #include <TH1F.h>
// using namespace o2;
// using namespace o2:framework;
// class ATask : public AnalysisTask
// {
// public:
//  OutputObj<TH2F> hPhi{TH1F("phi", "Phi", 100, 0., 2. * M_PI, 102, -2.01, 2.01)};

//  void process(aod:Tracks const& tracks)
//  {
//  for (auto& track : tracks) {
//  float phi = asin(track.snp()) + track.alpha() + M_PI;
//  hPhi-Fill(phi);
//  }
//  }
// };
// WorkflowSpec defineDataProcessing(ConfigContext const&)
// {
//  return WorkflowSpec{
//  adaptAnalysisTask<ATask>("mySimpleTrackAnalysis", 0)
//  };
// }