#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <SimulationDataFormat/MCTrack.h>
#include <TFile.h>
#include <TTree.h>
#include "TH1D.h"
#include <vector>
#include "TDatabasePDG.h"
#include "TROOT.h"
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

struct ATask{
  TDatabasePDG *pdg = new TDatabasePDG();

  void run(ProcessingContext& pc) {
    // mc
    TFile f("o2sim.root");
    TTree* tree = (TTree*) f.Get("o2sim");
    std::vector<MCTrack>* mctracks;
    tree->SetBranchAddress("MCTrack", &mctracks);

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

    tracksTree->GetEvent(0);
    treeTofMatch->GetEvent(0);

    // loop over mc events
    for (Int_t ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      Int_t nTracks = mctracks->size();

      // loop over mc tracks
      for (UInt_t i=0; i<nTracks; i++){
        MCTrack &track = (*mctracks)[i];

        if (track.getMotherTrackId() != -1)
          continue;

        TParticlePDG *particle = pdg->GetParticle(track.GetPdgCode());
        if ((int)(particle->Charge()) == 0)
          continue;

        if (abs(track.GetRapidity())<1)
          mcID.push_back(i);
      } // mc tracks

      if (mcID.size() == 2) {
        MCTrack &track1 = (*mctracks)[mcID[0]];
        MCTrack &track2 = (*mctracks)[mcID[1]];
        px = track1.GetStartVertexMomentumX();
        py = track1.GetStartVertexMomentumY();
        pz = track1.GetStartVertexMomentumZ();
        m = track1.GetMass();

        p[0] = px;
        p[1] = py;
        p[2] = pz;
        pl1.SetVectM(p, m);

        px = track2.GetStartVertexMomentumX();
        py = track2.GetStartVertexMomentumY();
        pz = track2.GetStartVertexMomentumZ();
        m = track2.GetMass();

        p[0] = px;
        p[1] = py;
        p[2] = pz;
        pl2.SetVectM(p, m);

        pl = pl1 + pl2;

        invM = pl.M();
        hMass->Fill(invM);

        for (int imc = 0; imc<mcID.size(); imc++) {
          for (int j=0; j<tpcItsTracks->size(); j++){
            o2::MCCompLabel trLabel = mcMatchTPC->at(j);
            int EvID = trLabel.getEventID();
            if (EvID != ev)
              continue;

            int ID = trLabel.getTrackID();
            if (ID == mcID[imc]) {
              recoInd.push_back(j);
              countRec++;

              for (int k=0; k<TOFMatchInfo->size(); k++){
                o2::dataformats::MatchInfoTOF infoTOF = TOFMatchInfo->at(k);
                int tofTrInd = infoTOF.getTrackIndex();
                if (j == tofTrInd) {
                  countTof++;
                }
              }
            }
          }
        }
      }

      if (countRec==2 && countTof==2) {
        hMassReco2->Fill(invM);
      }

      if (countRec==2 && countTof>=1) {
        hMassReco1->Fill(invM);
      }

      if (countRec==2) {
        hMassReco->Fill(invM);
      }

      countRec = 0;
      countTof = 0;
      mcID.clear();
      recoInd.clear();
      pl.Clear();
      pl1.Clear();
      pl2.Clear();
      p.Clear();
      px = 0.0;
      py = 0.0;
      pz = 0.0;
      m = 0.0;
      invM = 0.0;
    } // mc events

    TFile fig("eff.root","recreate");
    hMass->Write();
    hMassReco->Write();
    hMassReco1->Write();
    hMassReco2->Write();
    fig.Write();

    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&) {
  return WorkflowSpec{
    adaptAnalysisTask<ATask>("efficiency-analysis")
  };
}