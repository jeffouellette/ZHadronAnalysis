#ifndef __MinbiasAnalysis_h__
#define __MinbiasAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace atlashi;

class MinbiasAnalysis : public FullAnalysis {

  private:
  TFile* zMixFile = nullptr;
  TFile* mixedEventsFile = nullptr;
  TTree* mixedEventsTree = nullptr;

  TTree* LoadEventMixingTree (const char* _inFile, const char* _treeName);

  public:

  MinbiasAnalysis (const char* _name = "bkg") : FullAnalysis () {
    name = _name;
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    hasBkg = false;
    backgroundSubtracted = true;
    histsUnfolded = true;
    iaaCalculated = true;
    //icpCalculated = true;
  }

  void LoadEventWeights () override;
  void GenerateWeights (const char* inFileName = "outFile.root", const char* outFileName = "eventWeightsFile.root");

  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  //void CreateHists () override;
  //void CopyAnalysis (PhysicsAnalysis* a, const bool copyBkgs = false) override;
  //void CombineHists ();
  //void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true);
  //void SaveHists (const char* histFileName = "savedHists.root");
  //void ScaleHists ();
  void Execute (const char* inFileName, const char* outFileName) override {
    cout << "Error: In MinbiasAnalysis :: Execute: Called invalid function! A third argument is required to specify the mixing file name. Exiting." << endl;
  }
  void Execute (const char* inFileName, const char* mbInFileName, const char* outFileName);
};




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load histograms into memory, then combine channels.
//////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  FullAnalysis :: LoadHists (histFileName, _finishHists);

  PhysicsAnalysis :: CombineHists ();
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load the event weights into memory
//////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: LoadEventWeights () {
  if (eventWeightsLoaded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory;

  eventWeightsFile = new TFile (Form ("%s/MinbiasAnalysis/Nominal/eventWeightsFile.root", rootPath.Data ()), "read");

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iFinerCent = 0; iFinerCent < numFinerCentBins; iFinerCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_PbPbQ2_weights[iSpc][iFinerCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_z_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFinerCent, iPtZ, name.c_str ()));
        h_PbPbPsi2_weights[iSpc][iFinerCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_z_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFinerCent, iPtZ, name.c_str ()));
      }
    }
  }

  eventWeightsLoaded = true;

  _gDirectory-> cd ();
  return;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Create new histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: CreateHists () {
//
//  FullAnalysis :: CreateHists ();
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_pt_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_z_trk_xzh_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", nXHZBins[iPtZ], xHZBins[iPtZ]);
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
//        }
//      }
//    }
//  }
//
//  histsLoaded = true;
//  return;
//}
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
//// Fill combined species histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: CombineHists () {
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 2; iSpc++) {
//      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          if (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_pt_corr[2][iPtZ][iPhi][iCent]->Add (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]);
//          if (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_xzh_corr[2][iPtZ][iPhi][iCent]->Add (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]);
//
//          if (iPhi != 0) {
//            if (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zpt_corr[iSpc][iPtZ][iCent]->Add (h_z_trk_zpt_corr[iSpc][iPtZ][iPhi][iCent]);
//            if (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zpt_corr[2][iPtZ][iCent]->Add (h_z_trk_zpt_corr[iSpc][iPtZ][iPhi][iCent]);
//            if (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zxzh_corr[iSpc][iPtZ][iCent]->Add (h_z_trk_zxzh_corr[iSpc][iPtZ][iPhi][iCent]);
//            if (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]) h_z_trk_zxzh_corr[2][iPtZ][iCent]->Add (h_z_trk_zxzh_corr[iSpc][iPtZ][iPhi][iCent]);
//          }
//        } // end loop over phi
//      } // end loop over pT^Z
//    } // end loop over species
//  } // end loop over centralities
//  return;
//}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Scale histograms for plotting, calculating signals, etc.
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: ScaleHists () {
//  if (histsScaled || !histsLoaded)
//    return;
//
//  FullAnalysis :: ScaleHists ();
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
//        const float counts = countsHist->GetBinContent (1);
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          const double countsdPhi = counts * (phiHighBins[iPhi]-phiLowBins[iPhi]);
//
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]->Scale (1. / countsdPhi, "width");
//        }
//      }
//    }
//  }
//}
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load pre-filled histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
//  if (histsLoaded)
//    return;
//
//  FullAnalysis :: LoadHists (histFileName, false);
//
//  if (!histFile) {
//    SetupDirectories ("", "ZTrackAnalysis/");
//    histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");
//  }
//  TDirectory* _gDirectory = gDirectory;
//  if (!histFile->IsOpen ()) {
//    cout << "Error in MinbiasAnalysis :: LoadHists: histFile not open after calling parent function, exiting." << endl;
//    return;
//  }
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//          h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_corr_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        }
//      }
//    }
//  }
//
//  _gDirectory->cd ();
//
//  histsLoaded = true;
//
//  if (_finishHists) {
//    MinbiasAnalysis :: CombineHists ();
//    MinbiasAnalysis :: ScaleHists ();
//  }
//
//  return;
//}
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
//// Save histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void MinbiasAnalysis :: SaveHists (const char* histFileName) {
//  FullAnalysis :: SaveHists (histFileName);
//
//  if (!histFile) {
//    SetupDirectories ("", "ZTrackAnalysis/");
//    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "update");
//    histFile->cd ();
//  }
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          SafeWrite (h_z_trk_pt_corr[iSpc][iPtZ][iPhi][iCent]);
//          SafeWrite (h_z_trk_xzh_corr[iSpc][iPtZ][iPhi][iCent]);
//        }
//      }
//    }
//  }
//
//  histFile->Close ();
//  histFile = nullptr;
//  histsLoaded = false;
//  return;
//}




////////////////////////////////////////////////////////////////////////////////////////////////
// Loads the file for mixing Z's into minimum bias events
////////////////////////////////////////////////////////////////////////////////////////////////
TTree* MinbiasAnalysis :: LoadEventMixingTree (const char* _inFile, const char* _treeName) {
  if (zMixFile && zMixFile->IsOpen ())
    zMixFile->Close ();

  SetupDirectories ("", "ZTrackAnalysis/");

  TString inFile = TString (_inFile);
  if (!inFile.Contains (".root"))
    inFile = inFile + ".root";
  inFile.Replace (0, 7, "Data");

  zMixFile = new TFile (Form ("%s/%s", rootPath.Data (), inFile.Data ()), "read");

  if (!zMixFile || !zMixFile->IsOpen ()) {
    cout << "Cannot find zMixFile! Will return null tree!" << endl;
    return nullptr;
  }

  mixedEventsFile = new TFile (Form ("%s/%s_mixed.root", rootPath.Data (), inFile.Data ()), "recreate");

  TTree* _tree = (TTree*) zMixFile->Get (_treeName);
  return _tree;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: Execute (const char* inFileName, const char* mbInFileName, const char* outFileName) {

  LoadEventWeights ();

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* mbInFile = new TFile (Form ("%s/%s", rootPath.Data (), mbInFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), mbInFileName) << endl;

  TTree* mbPbPbTree = (TTree*) mbInFile->Get ("PbPbZTrackTree");
  TTree* mbppTree = (TTree*) mbInFile->Get ("ppZTrackTree");

  CreateHists ();

  unsigned int event_number = 0, z_event_number = 0, lumi_block = 0;
  int run_number = 0, z_run_number = 0, ntrk = 0, z_ntrk = 0;
  bool isEE = false;//, passes_toroid = false;
  float z_event_weight = 1;// q2_weight = 1, psi2_weight = 1; // vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, zdcEnergy = 0, z_fcal_et = 0, z_q2 = 0, z_psi2 = 0, z_vz = 0, z_zdcEnergy = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  vector<float>* z_trk_pt = nullptr, *z_trk_eta = nullptr, *z_trk_phi = nullptr, *z_trk_charge = nullptr;
  vector<float> out_z_trk_pt (0), out_z_trk_eta (0), out_z_trk_phi (0), out_z_trk_charge (0);

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (mbPbPbTree) {
    const int nMBEvts = mbPbPbTree->GetEntries ();
    mbPbPbTree->SetBranchAddress ("run_number",   &run_number);
    mbPbPbTree->SetBranchAddress ("event_number", &event_number);
    mbPbPbTree->SetBranchAddress ("lumi_block",   &lumi_block);
    //mbPbPbTree->SetBranchAddress ("passes_toroid",&passes_toroid);
    mbPbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    mbPbPbTree->SetBranchAddress ("zdcEnergy",    &zdcEnergy);
    mbPbPbTree->SetBranchAddress ("q2",           &q2);
    mbPbPbTree->SetBranchAddress ("psi2",         &psi2);
    mbPbPbTree->SetBranchAddress ("vz",           &vz);
    mbPbPbTree->SetBranchAddress ("ntrk",         &ntrk);
    mbPbPbTree->SetBranchAddress ("trk_pt",       trk_pt);
    mbPbPbTree->SetBranchAddress ("trk_eta",      trk_eta);
    mbPbPbTree->SetBranchAddress ("trk_phi",      trk_phi);

    mbPbPbTree->LoadBaskets (8000000000); //2000000000 = 2GB


    int iMBEvt = 0;
    mbPbPbTree->GetEntry (iMBEvt);

    std::vector <int> mbEventOrder = {};
    std::vector <bool> mbEventsUsed = {};
    for (int i = 0; i < nMBEvts; i++) {
      mbEventOrder.push_back (i);
      mbEventsUsed.push_back (false);
    }
    //std::srand (std::time (0));
    std::random_shuffle (mbEventOrder.begin (), mbEventOrder.end ());

    TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
    //TTree* zTree = LoadEventMixingTree (inFileName, "PbPbZTrackTree");
    TTree* zTree = (TTree*) inFile->Get ("PbPbZTrackTree");
    if (!zTree)
      cout << "Got a null Z boson tree!" << endl;
    int nZEvts = zTree->GetEntries ();
    zTree->LoadBaskets (4000000000);
    zTree->SetBranchAddress ("isEE",          &isEE);
    zTree->SetBranchAddress ("z_pt",          &z_pt);
    zTree->SetBranchAddress ("z_y",           &z_y);
    zTree->SetBranchAddress ("z_phi",         &z_phi);
    zTree->SetBranchAddress ("z_m",           &z_m);
    zTree->SetBranchAddress ("ntrk",          &z_ntrk);
    zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
    zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
    zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    zTree->SetBranchAddress ("trk_charge",    &z_trk_charge);
    zTree->SetBranchAddress ("l1_pt",         &l1_pt);
    zTree->SetBranchAddress ("l1_eta",        &l1_eta);
    zTree->SetBranchAddress ("l1_phi",        &l1_phi);
    zTree->SetBranchAddress ("l1_trk_pt",     &l1_trk_pt);
    zTree->SetBranchAddress ("l1_trk_eta",    &l1_trk_eta);
    zTree->SetBranchAddress ("l1_trk_phi",    &l1_trk_phi);
    zTree->SetBranchAddress ("l1_charge",     &l1_charge);
    zTree->SetBranchAddress ("l2_pt",         &l2_pt);
    zTree->SetBranchAddress ("l2_eta",        &l2_eta);
    zTree->SetBranchAddress ("l2_phi",        &l2_phi);
    zTree->SetBranchAddress ("l2_trk_pt",     &l2_trk_pt);
    zTree->SetBranchAddress ("l2_trk_eta",    &l2_trk_eta);
    zTree->SetBranchAddress ("l2_trk_phi",    &l2_trk_phi);
    zTree->SetBranchAddress ("l2_charge",     &l2_charge);
    zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
    zTree->SetBranchAddress ("zdcEnergy",     &z_zdcEnergy);
    zTree->SetBranchAddress ("q2",            &z_q2);
    zTree->SetBranchAddress ("psi2",          &z_psi2);
    zTree->SetBranchAddress ("vz",            &z_vz);
    zTree->SetBranchAddress ("event_number",  &z_event_number);
    zTree->SetBranchAddress ("run_number",    &z_run_number);
    //zTree->SetBranchAddress ("lumi_block",    &z_lumi_block);
    zTree->SetBranchAddress ("event_weight",  &z_event_weight);


    TFile* mixedEventsFile = new TFile (Form ("%s/%s_mixed.root", rootPath.Data (), inFileName), "recreate");
    TTree* mixedEventsTree = new TTree ("PbPbMixedTree", "PbPbMixedTree");
    mixedEventsTree->Branch ("run_number",    &run_number,    "run_number/i");
    mixedEventsTree->Branch ("event_number",  &event_number,  "event_number/i");
    mixedEventsTree->Branch ("lumi_block",    &lumi_block,    "lumi_block/i");
    mixedEventsTree->Branch ("isEE",          &isEE,          "isEE/O");

    mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/i");
    mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/i");
    //mixedEventsTree->Branch ("z_lumi_block",    &z_lumi_block,    "z_lumi_block/i");
    mixedEventsTree->Branch ("z_fcal_et",       &z_fcal_et,       "z_fcal_et/F");
    mixedEventsTree->Branch ("z_zdcEnergy",     &z_zdcEnergy,     "z_zdcEnergy/F");
    mixedEventsTree->Branch ("z_q2",            &z_q2,            "z_q2/F");
    mixedEventsTree->Branch ("z_psi2",          &z_psi2,          "z_psi2/F");
    mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
    mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
    mixedEventsTree->Branch ("z_trk_pt",        &out_z_trk_pt);
    mixedEventsTree->Branch ("z_trk_eta",       &out_z_trk_eta);
    mixedEventsTree->Branch ("z_trk_phi",       &out_z_trk_phi);
    mixedEventsTree->Branch ("z_trk_charge",    &out_z_trk_charge);

    mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
    mixedEventsTree->Branch ("z_pt",            &z_pt,            "z_pt/F");
    mixedEventsTree->Branch ("z_y",             &z_y,             "z_y/F");
    mixedEventsTree->Branch ("z_phi",           &z_phi,           "z_phi/F");
    mixedEventsTree->Branch ("z_m",             &z_m,             "z_m/F");
    mixedEventsTree->Branch ("l1_pt",           &l1_pt,           "l1_pt/F");
    mixedEventsTree->Branch ("l1_eta",          &l1_eta,          "l1_eta/F");
    mixedEventsTree->Branch ("l1_phi",          &l1_phi,          "l1_phi/F");
    mixedEventsTree->Branch ("l1_trk_pt",       &l1_trk_pt,       "l1_trk_pt/F");
    mixedEventsTree->Branch ("l1_trk_eta",      &l1_trk_eta,      "l1_trk_eta/F");
    mixedEventsTree->Branch ("l1_trk_phi",      &l1_trk_phi,      "l1_trk_phi/F");
    mixedEventsTree->Branch ("l1_charge",       &l1_charge,       "l1_charge/I");
    mixedEventsTree->Branch ("l2_pt",           &l2_pt,           "l2_pt/F");
    mixedEventsTree->Branch ("l2_eta",          &l2_eta,          "l2_eta/F");
    mixedEventsTree->Branch ("l2_phi",          &l2_phi,          "l2_phi/F");
    mixedEventsTree->Branch ("l2_trk_pt",       &l2_trk_pt,       "l2_trk_pt/F");
    mixedEventsTree->Branch ("l2_trk_eta",      &l2_trk_eta,      "l2_trk_eta/F");
    mixedEventsTree->Branch ("l2_trk_phi",      &l2_trk_phi,      "l2_trk_phi/F");
    mixedEventsTree->Branch ("l2_charge",       &l2_charge,       "l2_charge/I");

    mixedEventsTree->Branch ("fcal_et",       &fcal_et,       "fcal_et/F");
    mixedEventsTree->Branch ("zdcEnergy",     &zdcEnergy,     "zdcEnergy/F");
    mixedEventsTree->Branch ("q2",            &q2,            "q2/F");
    mixedEventsTree->Branch ("psi2",          &psi2,          "psi2/F");
    mixedEventsTree->Branch ("vz",            &vz,            "vz/F");
    mixedEventsTree->Branch ("ntrk",          &ntrk,          "ntrk/I");

    mixedEventsTree->Branch ("trk_pt",        &trk_pt,        "trk_pt[ntrk]/F");
    mixedEventsTree->Branch ("trk_eta",       &trk_eta,       "trk_eta[ntrk]/F");
    mixedEventsTree->Branch ("trk_phi",       &trk_phi,       "trk_phi[ntrk]/F");


    if (nZEvts == 0)
      cout << "Warning! No Z's to mix with in this run!" << endl;
    cout << "For this PbPb tree, maximum mixing fraction = " << nMBEvts / nZEvts << endl;

    bool doShuffle = false;
    if (mixingFraction * nZEvts > nMBEvts) {
      cout << "Warning! Mixing fraction too high, will use " << (float)(nMBEvts / mixingFraction) / (float)(nZEvts) * 100. << "% of Z events" << endl;
      nZEvts = nMBEvts / mixingFraction;
      doShuffle = true;
    }

    std::vector <int> zEventOrder = {};
    std::vector <int> zEventsUsed = {};
    for (int i = 0; i < nZEvts; i++) {
      zEventOrder.push_back (i);
      zEventsUsed.push_back (false);
    }
    //std::srand (std::time (0));
    if (doShuffle) std::random_shuffle (zEventOrder.begin (), zEventOrder.end ());

    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (zEventOrder[iZEvt % nZEvts]);
      zEventsUsed[iZEvt % nZEvts]++;

      //if (fabs (z_vz) > 150)
      //  continue;

      if (z_event_weight == 0)
        continue;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;
      if (iPtZ < 2)
        continue; // no one cares about these events anyways

      {
        const short iFCalEt = GetSuperFineCentBin (z_fcal_et);
        if (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1)
          continue;
        //const short iPsi2 = GetPsi2Bin (z_psi2);
        //if (iPsi2 < 0 || iPsi2 > numPsi2Bins-1)
        //  continue;
        //const short iVZ = GetVZBin (z_vz);
        //if (iVZ < 0 || iVZ > numVZBins)
        //  continue;
        
        bool goodMixEvent = false;
        const int _iMBEvt = iMBEvt;
        do {
          iMBEvt = (iMBEvt+1) % nMBEvts;
          mbPbPbTree->GetEntry (mbEventOrder[iMBEvt]);
          goodMixEvent = (!mbEventsUsed[iMBEvt] && iFCalEt == GetSuperFineCentBin (fcal_et));// && fabs (vz) <= 150);// && iVZ == GetVZBin (vz));
        } while (!goodMixEvent && iMBEvt != _iMBEvt);
        if (_iMBEvt == iMBEvt) {
          cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          continue;
        }
        mbEventsUsed[iMBEvt] = true;
      }

      out_z_trk_pt.clear ();
      out_z_trk_eta.clear ();
      out_z_trk_phi.clear ();
      out_z_trk_charge.clear ();
      for (int iTrk = 0; iTrk < z_ntrk; iTrk++) {
        out_z_trk_pt.push_back (z_trk_pt->at (iTrk));
        out_z_trk_eta.push_back (z_trk_eta->at (iTrk));
        out_z_trk_phi.push_back (z_trk_phi->at (iTrk));
        out_z_trk_charge.push_back (z_trk_charge->at (iTrk));
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      mixedEventsTree->Fill ();

      const short iSpc = (isEE ? 0 : 1); // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      const short iFinerCent = GetFinerCentBin (fcal_et);
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      //{
      //  float dphi = DeltaPhi (z_phi, psi2, false);
      //  if (dphi > pi/2)
      //    dphi = pi - dphi;
      //  //q2_weight = h_PbPbQ2_weights[iSpc][iFinerCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFinerCent][iPtZ]->FindBin (q2));
      //  psi2_weight = h_PbPbPsi2_weights[iSpc][iFinerCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFinerCent][iPtZ]->FindBin (dphi));

      //  //z_event_weight = fcal_weight * q2_weight * psi2_weight * vz_weight;
      //  z_event_weight *= psi2_weight;
      //}
      //z_event_weight = 1;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, z_event_weight);
      h_q2[iFinerCent]->Fill (q2);
      h_q2_reweighted[iFinerCent]->Fill (q2, z_event_weight);
      h_psi2[iFinerCent]->Fill (psi2);
      h_psi2_reweighted[iFinerCent]->Fill (psi2, z_event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, z_event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, z_event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;
        //if (trk_eta[iTrk] < 0)
        //  continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkPur == 0 || trkEff == 0)
          continue;
        const double trkWeight = z_event_weight * trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_z_trk_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
          }
        }
      }

    } // end loop over Pb+Pb tree

    mixedEventsTree->SetDirectory (mixedEventsFile);
    mixedEventsTree->Write ("", TObject :: kOverwrite);

    mixedEventsFile->Close ();

    inFile->Close ();

    cout << "Done minbias Pb+Pb loop." << endl;
  }


  if (mbppTree) {
    const int nMBEvts = mbppTree->GetEntries ();
    mbppTree->SetBranchAddress ("run_number",   &run_number);
    mbppTree->SetBranchAddress ("event_number", &event_number);
    mbppTree->SetBranchAddress ("lumi_block",   &lumi_block);
    mbppTree->SetBranchAddress ("vz",           &vz);
    mbppTree->SetBranchAddress ("ntrk",         &ntrk);
    mbppTree->SetBranchAddress ("trk_pt",       trk_pt);
    mbppTree->SetBranchAddress ("trk_eta",      trk_eta);
    mbppTree->SetBranchAddress ("trk_phi",      trk_phi);
    mbppTree->LoadBaskets (3000000000);


    int iMBEvt = 0;
    mbppTree->GetEntry (iMBEvt);

    std::vector <int> mbEventOrder = {};
    std::vector <bool> mbEventsUsed = {};
    for (int i = 0; i < nMBEvts; i++) {
      mbEventOrder.push_back (i);
      mbEventsUsed.push_back (false);
    }
    //std::srand (std::time (0));
    std::random_shuffle (mbEventOrder.begin (), mbEventOrder.end ());


    TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
    //TTree* zTree = LoadEventMixingTree (inFileName, "ppZTrackTree");
    TTree* zTree = (TTree*) inFile->Get ("ppZTrackTree");
    if (!zTree)
      cout << "Got a null Z boson tree!" << endl;
    int nZEvts = zTree->GetEntries ();
    zTree->LoadBaskets (3000000000);
    zTree->SetBranchAddress ("isEE",          &isEE);
    zTree->SetBranchAddress ("z_pt",          &z_pt);
    zTree->SetBranchAddress ("z_y",           &z_y);
    zTree->SetBranchAddress ("z_phi",         &z_phi);
    zTree->SetBranchAddress ("z_m",           &z_m);
    zTree->SetBranchAddress ("ntrk",          &z_ntrk);
    zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
    zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
    zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    zTree->SetBranchAddress ("trk_charge",    &z_trk_charge);
    zTree->SetBranchAddress ("l1_pt",         &l1_pt);
    zTree->SetBranchAddress ("l1_eta",        &l1_eta);
    zTree->SetBranchAddress ("l1_phi",        &l1_phi);
    zTree->SetBranchAddress ("l1_trk_pt",     &l1_trk_pt);
    zTree->SetBranchAddress ("l1_trk_eta",    &l1_trk_eta);
    zTree->SetBranchAddress ("l1_trk_phi",    &l1_trk_phi);
    zTree->SetBranchAddress ("l1_charge",     &l1_charge);
    zTree->SetBranchAddress ("l2_pt",         &l2_pt);
    zTree->SetBranchAddress ("l2_eta",        &l2_eta);
    zTree->SetBranchAddress ("l2_phi",        &l2_phi);
    zTree->SetBranchAddress ("l2_trk_pt",     &l2_trk_pt);
    zTree->SetBranchAddress ("l2_trk_eta",    &l2_trk_eta);
    zTree->SetBranchAddress ("l2_trk_phi",    &l2_trk_phi);
    zTree->SetBranchAddress ("l2_charge",     &l2_charge);
    zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
    zTree->SetBranchAddress ("zdcEnergy",     &z_zdcEnergy);
    zTree->SetBranchAddress ("q2",            &z_q2);
    zTree->SetBranchAddress ("psi2",          &z_psi2);
    zTree->SetBranchAddress ("vz",            &z_vz);
    zTree->SetBranchAddress ("event_number",  &z_event_number);
    zTree->SetBranchAddress ("run_number",    &z_run_number);
    //zTree->SetBranchAddress ("lumi_block",    &z_lumi_block);
    zTree->SetBranchAddress ("event_weight",  &z_event_weight);


    TFile* mixedEventsFile = new TFile (Form ("%s/%s_mixed.root", rootPath.Data (), inFileName), "recreate");
    TTree* mixedEventsTree = new TTree ("ppMixedTree", "ppMixedTree");
    mixedEventsTree->Branch ("run_number",    &run_number,    "run_number/i");
    mixedEventsTree->Branch ("event_number",  &event_number,  "event_number/i");
    mixedEventsTree->Branch ("lumi_block",    &lumi_block,    "lumi_block/i");
    mixedEventsTree->Branch ("isEE",          &isEE,          "isEE/O");
    mixedEventsTree->Branch ("vz",            &vz,            "vz/F");
    mixedEventsTree->Branch ("ntrk",          &ntrk,          "ntrk/I");

    mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/i");
    mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/i");
    //mixedEventsTree->Branch ("z_lumi_block",    &z_lumi_block,    "z_lumi_block/i");
    mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
    mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
    mixedEventsTree->Branch ("z_trk_pt",        &out_z_trk_pt);
    mixedEventsTree->Branch ("z_trk_eta",       &out_z_trk_eta);
    mixedEventsTree->Branch ("z_trk_phi",       &out_z_trk_phi);
    mixedEventsTree->Branch ("z_trk_charge",    &out_z_trk_charge);

    mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
    mixedEventsTree->Branch ("z_pt",            &z_pt,            "z_pt/F");
    mixedEventsTree->Branch ("z_y",             &z_y,             "z_y/F");
    mixedEventsTree->Branch ("z_phi",           &z_phi,           "z_phi/F");
    mixedEventsTree->Branch ("z_m",             &z_m,             "z_m/F");
    mixedEventsTree->Branch ("l1_pt",           &l1_pt,           "l1_pt/F");
    mixedEventsTree->Branch ("l1_eta",          &l1_eta,          "l1_eta/F");
    mixedEventsTree->Branch ("l1_phi",          &l1_phi,          "l1_phi/F");
    mixedEventsTree->Branch ("l1_trk_pt",       &l1_trk_pt,       "l1_trk_pt/F");
    mixedEventsTree->Branch ("l1_trk_eta",      &l1_trk_eta,      "l1_trk_eta/F");
    mixedEventsTree->Branch ("l1_trk_phi",      &l1_trk_phi,      "l1_trk_phi/F");
    mixedEventsTree->Branch ("l1_charge",       &l1_charge,       "l1_charge/I");
    mixedEventsTree->Branch ("l2_pt",           &l2_pt,           "l2_pt/F");
    mixedEventsTree->Branch ("l2_eta",          &l2_eta,          "l2_eta/F");
    mixedEventsTree->Branch ("l2_phi",          &l2_phi,          "l2_phi/F");
    mixedEventsTree->Branch ("l2_trk_pt",       &l2_trk_pt,       "l2_trk_pt/F");
    mixedEventsTree->Branch ("l2_trk_eta",      &l2_trk_eta,      "l2_trk_eta/F");
    mixedEventsTree->Branch ("l2_trk_phi",      &l2_trk_phi,      "l2_trk_phi/F");
    mixedEventsTree->Branch ("l2_charge",       &l2_charge,       "l2_charge/I");

    mixedEventsTree->Branch ("trk_pt",        &trk_pt,        "trk_pt[ntrk]/F");
    mixedEventsTree->Branch ("trk_eta",       &trk_eta,       "trk_eta[ntrk]/F");
    mixedEventsTree->Branch ("trk_phi",       &trk_phi,       "trk_phi[ntrk]/F");

    if (nZEvts == 0)
      cout << "Warning! No Z's to mix with in this run!" << endl;
    cout << "For this pp tree, maximum mixing fraction = " << nMBEvts / nZEvts << endl;

    bool doShuffle = false;
    if (mixingFraction * nZEvts > nMBEvts) {
      cout << "Warning! Mixing fraction too high, will use " << (float)(nMBEvts / mixingFraction) / (float)(nZEvts) * 100. << "% of Z events" << endl;
      nZEvts = nMBEvts / mixingFraction;
      doShuffle = true;
    }

    std::vector <int> zEventOrder = {};
    std::vector <int> zEventsUsed = {};
    for (int i = 0; i < nZEvts; i++) {
      zEventOrder.push_back (i);
      zEventsUsed.push_back (false);
    }
    //std::srand (std::time (0));
    if (doShuffle) std::random_shuffle (zEventOrder.begin (), zEventOrder.end ());

    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (zEventOrder[iZEvt % nZEvts]);
      zEventsUsed[iZEvt % nZEvts]++;

      //if (fabs (z_vz) > 150)
      //  continue;

      if (z_event_weight == 0)
        continue;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;
      if (iPtZ < 2)
        continue; // no one cares about these events anyways

      {
        bool goodMixEvent = false;
        const int _iMBEvt = iMBEvt;
        do {
          iMBEvt = (iMBEvt+1) % nMBEvts;
          mbppTree->GetEntry (mbEventOrder[iMBEvt]);
          goodMixEvent = (!mbEventsUsed[iMBEvt]);
        } while (!goodMixEvent && iMBEvt != _iMBEvt);
        if (_iMBEvt == iMBEvt) {
          cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          continue;
        }
        mbEventsUsed[iMBEvt] = true;
      }

      out_z_trk_pt.clear ();
      out_z_trk_eta.clear ();
      out_z_trk_phi.clear ();
      out_z_trk_charge.clear ();
      for (int iTrk = 0; iTrk < z_ntrk; iTrk++) {
        out_z_trk_pt.push_back (z_trk_pt->at (iTrk));
        out_z_trk_eta.push_back (z_trk_eta->at (iTrk));
        out_z_trk_phi.push_back (z_trk_phi->at (iTrk));
        out_z_trk_charge.push_back (z_trk_charge->at (iTrk));
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      mixedEventsTree->Fill ();

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      //z_event_weight = 1;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, z_event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, z_event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, z_event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;
        //if (trk_eta[iTrk] < 0)
        //  continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const double trkWeight = z_event_weight * trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_z_trk_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
          }
        }
      }

    } // end loop over pp tree

    mixedEventsTree->SetDirectory (mixedEventsFile);
    mixedEventsTree->Write ("", TObject :: kOverwrite);

    mixedEventsFile->Close ();
    SaferDelete (mixedEventsFile);

    inFile->Close ();
    SaferDelete (inFile);

    cout << "Done minbias pp loop." << endl;
  }

  SaveHists (outFileName);

  mbInFile->Close ();
  SaferDelete (mbInFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Generates weights from pre-mixed events.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: GenerateWeights (const char* inFileName, const char* outFileName) {

  const int nQ2Bins = 20;
  const double* q2Bins = linspace (0, 0.3, nQ2Bins);
  const int nPsi2Bins = 8;
  const double* psi2Bins = linspace (0, pi/2, nPsi2Bins);
  
  //const int nNchBins = 160;
  //const double* nchBins = linspace (-0.5, 160.5, nNchBins);

  TH1D* h_z_q2_dist[3][numFinerCentBins][nPtZBins];
  TH1D* h_mixed_q2_dist[3][numFinerCentBins][nPtZBins];
  TH1D* h_z_psi2_dist[3][numFinerCentBins][nPtZBins];
  TH1D* h_mixed_psi2_dist[3][numFinerCentBins][nPtZBins];

  TFile* inFile = new TFile (inFileName, "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbMixedTree");

  for (int iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (int iFinerCent = 0; iFinerCent < numFinerCentBins; iFinerCent++) {
        h_z_q2_dist[iSpc][iFinerCent][iPtZ] = new TH1D (Form ("h_z_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFinerCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_mixed_q2_dist[iSpc][iFinerCent][iPtZ] = new TH1D (Form ("h_mixed_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFinerCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_z_q2_dist[iSpc][iFinerCent][iPtZ]->Sumw2 ();
        h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]->Sumw2 ();

        h_z_psi2_dist[iSpc][iFinerCent][iPtZ] = new TH1D (Form ("h_z_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFinerCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ] = new TH1D (Form ("h_mixed_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFinerCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->Sumw2 ();
        h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]->Sumw2 ();
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    bool isEE = false;
    float fcal_et = 0, z_pt = 0, z_phi = 0, q2 = 0, psi2 = 0, z_q2 = 0, z_psi2 = 0;

    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("z_q2",         &z_q2);
    PbPbTree->SetBranchAddress ("z_psi2",       &z_psi2);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iFinerCent = GetFinerCentBin (fcal_et);
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      const short iSpc = (isEE ? 0 : 1);

      h_z_q2_dist[iSpc][iFinerCent][iPtZ]->Fill (z_q2);
      h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]->Fill (q2);

      float dphi = DeltaPhi (z_phi, z_psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;
      h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->Fill (dphi);

      dphi = DeltaPhi (z_phi, psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;
      h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]->Fill (dphi);
    }
    cout << "Done Pb+Pb loop." << endl;
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iFinerCent = 0; iFinerCent < numFinerCentBins; iFinerCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int ix = 1; ix <= h_z_q2_dist[iSpc][iFinerCent][iPtZ]->GetNbinsX (); ix++)
          h_z_q2_dist[iSpc][iFinerCent][iPtZ]->SetBinError (ix, h_z_q2_dist[iSpc][iFinerCent][iPtZ]->GetBinError (ix) * TMath::Sqrt (40));
        for (int ix = 1; ix <= h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->GetNbinsX (); ix++)
          h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->SetBinError (ix, h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->GetBinError (ix) * TMath::Sqrt (40));

        h_z_q2_dist[2][iFinerCent][iPtZ]->Add (h_z_q2_dist[iSpc][iFinerCent][iPtZ]);
        h_mixed_q2_dist[2][iFinerCent][iPtZ]->Add (h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]);
        h_z_psi2_dist[2][iFinerCent][iPtZ]->Add (h_z_psi2_dist[iSpc][iFinerCent][iPtZ]);
        h_mixed_psi2_dist[2][iFinerCent][iPtZ]->Add (h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]);
      }
    }
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iFinerCent = 0; iFinerCent < numFinerCentBins; iFinerCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_q2_dist[iSpc][iFinerCent][iPtZ]->Scale (1./h_z_q2_dist[iSpc][iFinerCent][iPtZ]->Integral ());
        h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]->Scale (1./h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]->Integral ());
        h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->Scale (1./h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->Integral ());
        h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]->Scale (1./h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]->Integral ());
      }
    }
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iFinerCent = 0; iFinerCent < numFinerCentBins; iFinerCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_q2_dist[iSpc][iFinerCent][iPtZ]->Divide (h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]);
        h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]->Divide (h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]);
        h_z_psi2_dist[iSpc][iFinerCent][iPtZ]->Divide (h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]);
        h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]->Divide (h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]);
      }
    }
  }

  if (eventWeightsFile && eventWeightsFile->IsOpen ()) {
    eventWeightsFile->Close ();
    eventWeightsLoaded = false;
  }

  eventWeightsFile = new TFile (outFileName, "recreate");
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iFinerCent = 0; iFinerCent < numFinerCentBins; iFinerCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_z_q2_dist[iSpc][iFinerCent][iPtZ]);
        SafeWrite (h_mixed_q2_dist[iSpc][iFinerCent][iPtZ]);
        SafeWrite (h_z_psi2_dist[iSpc][iFinerCent][iPtZ]);
        SafeWrite (h_mixed_psi2_dist[iSpc][iFinerCent][iPtZ]);
      }
    }
  }
  eventWeightsFile->Close ();
}

#endif
