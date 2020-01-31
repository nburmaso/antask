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

using namespace std;
using namespace o2;

#define PI 3.14159265

void checkGen() {
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

  gROOT->cd();
  gStyle->SetOptStat(0);

  TH1D* hPt = new TH1D("hPt", "; Pt, GeV; N", 50, 0, 2);
  TH1D* hPhi = new TH1D("hPhi", "; phi, deg; N", 180, -90, 90);
  TH1D* hMass = new TH1D("hMass", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hEta = new TH1D("hEta", "; eta; N", 1000, -1, 1);

  TH1D* hPtReco = new TH1D("hPtReco", "; Pt, GeV; N", 50, 0, 2);
  TH1D* hPhiReco = new TH1D("hPhiReco", "; phi, deg; N", 180, -90, 90);
  TH1D* hMassReco = new TH1D("hMassReco", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hMassReco1 = new TH1D("hMassReco1", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hMassReco2 = new TH1D("hMassReco2", "; mass, GeV; N", 50, 0.5, 4);
  TH1D* hEtaReco = new TH1D("hEtaReco", "; eta; N", 1000, -1, 1);

  hPt->Sumw2();
  hPtReco->Sumw2();
  hMass->Sumw2();
  hMassReco->Sumw2();

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

      hPt->Fill( sqrt( pl[0]*pl[0] + pl[1]*pl[1]) );
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
      hPtReco->Fill(sqrt( pl[0]*pl[0] + pl[1]*pl[1]) );
      hMassReco2->Fill(invM);
    }

    if (countRec==2 && countTof>=1) {
      hPtReco->Fill(sqrt( pl[0]*pl[0] + pl[1]*pl[1]) );
      hMassReco1->Fill(invM);
    }

    if (countRec==2) {
      hPtReco->Fill(sqrt( pl[0]*pl[0] + pl[1]*pl[1]) );
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

// histograms
//  TCanvas* c1 = new TCanvas("c1", "Pt", 1800, 800);
//  c1->Divide(1,2);
//  c1->cd(1);
//  gPad->SetGrid(1,1);
  hPt->SetTitle("MC Pt");
  hPt->SetFillColor(kBlue);
  hPt->SetLineColor(kBlue);
  hPt->SetMarkerStyle(21);
  hPt->SetMarkerSize(0.7);

  // hPt->Draw("elp");

//  c1->cd(2);
//  gPad->SetGrid(1,1);
  hPtReco->SetTitle("Reco Pt");
  hPtReco->SetFillColor(kBlue);
  hPtReco->SetLineColor(kBlue);
  hPtReco->SetMarkerStyle(21);
  hPtReco->SetMarkerSize(0.7);

  // hPtReco->Draw("elp");
  // c1->Print("../graphs/pt_mc_vs_reco.png");

  // TH1D* hEff = (TH1D*)hPtReco->Clone("hEff");
  // TCanvas* c10 = new TCanvas("c10", "Eff", 1800, 800);
  // c10->cd();
  // gPad->SetGrid(1,1);
  // hEff->Divide(hEff, hPt, 1, 1, "b");
  // hEff->SetTitle("Eff.");
  // hEff->GetXaxis()->SetTitle("Pt, GeV");
  // hEff->GetYaxis()->SetTitle("eff.");
  // hEff->GetXaxis()->SetRange(1, 30);
  // hEff->SetFillColor(kBlue);
  // hEff->SetLineColor(kBlue);
  // hEff->SetMarkerStyle(21);
  // hEff->SetMarkerSize(0.7);

  // hEff->Draw("elp");
  // c10->Print("../graphs/pt_eff.png");

//  TCanvas* c2 = new TCanvas("c2", "Mass", 1800, 800);
//  c2->Divide(1,2);
//  c2->cd(1);
//  gPad->SetGrid(1,1);
//  gPad->SetLeftMargin(0.025);
//  gPad->SetRightMargin(0.03);
//  gPad->SetBottomMargin(0.10);
//  gPad->SetTopMargin(0.08);
  hMass->SetTitle("MC inv. mass");
  hMass->SetFillColor(kBlue);
  hMass->SetLineColor(kBlue);
  hMass->SetMarkerStyle(21);
  hMass->SetMarkerSize(0.7);

  hMass->GetYaxis()->SetTickLength(0.01);
  hMass->GetXaxis()->SetTitleSize(0.05);
  hMass->GetYaxis()->SetTitleSize(0.05);
  hMass->GetXaxis()->SetLabelSize(0.05);
  hMass->GetYaxis()->SetLabelSize(0.05);

  // hMass->Draw("elp");
/*
  c2->cd(2);
  gPad->SetGrid(1,1);
  gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.025);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetTopMargin(0.08);
*/
  hMassReco->SetTitle("Reco inv. mass");
  hMassReco->SetFillColor(kBlue);
  hMassReco->SetLineColor(kBlue);
  hMassReco->SetMarkerStyle(21);
  hMassReco->SetMarkerSize(0.7);

  hMassReco->GetYaxis()->SetTickLength(0.01);
  hMassReco->GetXaxis()->SetTitleSize(0.05);
  hMassReco->GetYaxis()->SetTitleSize(0.05);
  hMassReco->GetXaxis()->SetLabelSize(0.05);
  hMassReco->GetYaxis()->SetLabelSize(0.05);

  // hMassReco->Draw("elp");
  // c2->Print("../graphs/mass_mc_vs_reco.png");

  // TH1D* hMassEff = (TH1D*)hMassReco->Clone("hMassEff");
  // TCanvas* c20 = new TCanvas("c20", "Eff.", 1800, 800);
  // c20->cd();
  // gPad->SetGrid(1,1);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.04);
  // gPad->SetBottomMargin(0.10);
  // gPad->SetTopMargin(0.08);
  // hMassEff->Divide(hMassEff, hMass, 1, 1, "b");
  // hMassEff->SetTitle("Effectiveness of reconstruction");
  // hMassEff->GetXaxis()->SetTitle("inv. mass, GeV");
  // hMassEff->GetYaxis()->SetTitle("eff.");
  // hMassEff->GetXaxis()->SetRange(1, 50);
  // hMassEff->SetFillColor(kBlue);
  // hMassEff->SetLineColor(kBlue);
  // hMassEff->SetMarkerStyle(21);
  // hMassEff->SetMarkerSize(0.7);

  // hMassEff->GetYaxis()->SetTickLength(0.01);
  // hMassEff->GetXaxis()->SetTitleSize(0.04);
  // hMassEff->GetYaxis()->SetTitleSize(0.04);
  // hMassEff->GetXaxis()->SetLabelSize(0.025);
  // hMassEff->GetYaxis()->SetLabelSize(0.025);

  // hMassEff->Draw("elp");
  // c20->Print("../graphs/mass_eff.png");

  // TH1D* hMassEff1 = (TH1D*)hMassReco1->Clone("hMassEff1");
  // TCanvas* c21 = new TCanvas("c21", "Eff. 1 tof", 1800, 800);
  // c21->cd();
  // gPad->SetGrid(1,1);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.04);
  // gPad->SetBottomMargin(0.10);
  // gPad->SetTopMargin(0.08);
  // hMassEff1->Divide(hMassEff1, hMass, 1, 1, "b");
  // hMassEff1->SetTitle("Effectiveness of reconstruction, 1 track");
  // hMassEff1->GetXaxis()->SetTitle("inv. mass, GeV");
  // hMassEff1->GetYaxis()->SetTitle("eff.");
  // hMassEff1->GetXaxis()->SetRange(1, 50);
  // hMassEff1->SetFillColor(kBlue);
  // hMassEff1->SetLineColor(kBlue);
  // hMassEff1->SetMarkerStyle(21);
  // hMassEff1->SetMarkerSize(0.7);

  // hMassEff1->GetYaxis()->SetTickLength(0.01);
  // hMassEff1->GetXaxis()->SetTitleSize(0.04);
  // hMassEff1->GetYaxis()->SetTitleSize(0.04);
  // hMassEff1->GetXaxis()->SetLabelSize(0.025);
  // hMassEff1->GetYaxis()->SetLabelSize(0.025);

  // hMassEff1->Draw("elp");
  // // c21->Print("../graphs/mass_eff_1tof.png");

  // TH1D* hMassEff2 = (TH1D*)hMassReco2->Clone("hMassEff2");
  // TCanvas* c22 = new TCanvas("c22", "Eff. 2 tof", 1800, 800);
  // c22->cd();
  // gPad->SetGrid(1,1);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.04);
  // gPad->SetBottomMargin(0.10);
  // gPad->SetTopMargin(0.08);
  // hMassEff2->Divide(hMassEff2, hMass, 1, 1, "b");
  // hMassEff2->SetTitle("Effectiveness of reconstruction, 2 tracks");
  // hMassEff2->GetXaxis()->SetTitle("inv. mass, GeV");
  // hMassEff2->GetYaxis()->SetTitle("eff.");
  // hMassEff2->GetXaxis()->SetRange(1, 50);
  // hMassEff2->SetFillColor(kBlue);
  // hMassEff2->SetLineColor(kBlue);
  // hMassEff2->SetMarkerStyle(21);
  // hMassEff2->SetMarkerSize(0.7);

  // hMassEff2->GetYaxis()->SetTickLength(0.01);
  // hMassEff2->GetXaxis()->SetTitleSize(0.04);
  // hMassEff2->GetYaxis()->SetTitleSize(0.04);
  // hMassEff2->GetXaxis()->SetLabelSize(0.025);
  // hMassEff2->GetYaxis()->SetLabelSize(0.025);

  // hMassEff2->Draw("elp");
  // c22->Print("../graphs/mass_eff_2tof.png");
//

// saving to file
  TFile fig("eff.root","new");
  hPt->Write();
  hPtReco->Write();
  // hEff->Write();
  hMass->Write();
  hMassReco->Write();
  hMassReco1->Write();
  hMassReco2->Write();
  // hMassEff->Write();
  // hMassEff1->Write();
  // hMassEff2->Write();
  fig.Write();
//

}
