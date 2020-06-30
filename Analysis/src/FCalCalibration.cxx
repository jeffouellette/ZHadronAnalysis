#ifndef __FCalCalibration_cxx__
#define __FCalCalibration_cxx__

#include "FCalCalibration.h"
#include "Params.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"

#include <Utilities.h>
#include <AtlasUtils.h>

#include <TH2D.h>
#include <TChain.h>
#include <TSystem.h>

#include <iostream>

using namespace std;

namespace ZHadronAnalysis {

bool FCalCalibration (const char* directory,
                      const int dataSet,
                      const char* inFileName) {

  cout << "Info: In FCalCalibration.cxx: Entered TreeMaker routine." << endl;
  cout << "Info: In FCalCalibration.cxx: Printing systematic onfiguration:";

  SetupDirectories ("FCalCalibration");

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In FCalCalibration.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In FCalCalibration.cxx: Saving output to " << rootPath << endl;

  TTree* tree = nullptr;
  if (isMC) {
    TFile* file = GetFile (directory, dataSet, inFileName);
    tree = (TTree*) file->Get ("bush");
  }
  else {
    TChain* c = new TChain ("bush", "bush");
    TString pattern = "*.root";
    auto dir = gSystem->OpenDirectory (dataPath + directory);
    while (const char* f = gSystem->GetDirEntry (dir)) {
      TString file = TString (f);
      if (!file.Contains (Form ("%i", dataSet)))
        continue;
      cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << endl;
      c->Add (dataPath + directory + "/" + file + "/*.root");
      break;
    }
    cout << "Chain has " << c->GetListOfFiles ()->GetEntries () << " files, " << c->GetEntries () << " entries" << endl;
    tree = c;
  }
  if (tree == nullptr) {
    cout << "Error: In FCalCalibration.cxx: TChain not obtained for given data set. Quitting." << endl;
    return false;
  }

  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In FCalCalibration.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }


  // Input tree
  bool HLT_mb_sptrk_L1ZDC_A_C_VTE50 = false;
  bool HLT_noalg_pc_L1TE50_VTE600_0ETA49 = false;
  bool HLT_noalg_cc_L1TE600_0ETA49 = false;
  float HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale = 0;
  float HLT_noalg_pc_L1TE50_VTE600_0ETA49_prescale = 0;
  float HLT_noalg_cc_L1TE600_0ETA49_prescale = 0;
  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  if (isPbPb) t->SetGetZdc ();
  t->SetBranchAddresses ();
  tree->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",       &HLT_mb_sptrk_L1ZDC_A_C_VTE50);
  tree->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49",  &HLT_noalg_pc_L1TE50_VTE600_0ETA49);
  tree->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49",        &HLT_noalg_cc_L1TE600_0ETA49);
  tree->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale",       &HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale);
  tree->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49_prescale",  &HLT_noalg_pc_L1TE50_VTE600_0ETA49_prescale);
  tree->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49_prescale",        &HLT_noalg_cc_L1TE600_0ETA49_prescale);


  // Load files for output
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
  TH1D* h_PbPb_fcal_et = new TH1D ("h_PbPb_fcal_et", "", 500, 0, 5000);
  TH1D* h_pp_fcal_et = new TH1D ("h_pp_fcal_et", "", 150, -100, 200);
  TH1D* h_centrality_HLT_mb_sptrk_L1ZDC_A_C_VTE50       = new TH1D ("h_centrality_HLT_mb_sptrk_L1ZDC_A_C_VTE50",      "", 100, 0, 100);
  TH1D* h_centrality_HLT_noalg_pc_L1TE50_VTE600_0ETA49  = new TH1D ("h_centrality_HLT_noalg_pc_L1TE50_VTE600_0ETA49", "", 100, 0, 100);
  TH1D* h_centrality_HLT_noalg_cc_L1TE600_0ETA49        = new TH1D ("h_centrality_HLT_noalg_cc_L1TE600_0ETA49",       "", 100, 0, 100);
  TH1D* h_q2x_a = new TH1D ("h_q2x_a", "q_{x}^{a}", 100, -0.4, 0.4);
  TH1D* h_q2y_a = new TH1D ("h_q2y_a", "q_{y}^{a}", 100, -0.4, 0.4);
  TH1D* h_q2x_c = new TH1D ("h_q2x_c", "q_{x}^{c}", 100, -0.4, 0.4);
  TH1D* h_q2y_c = new TH1D ("h_q2y_c", "q_{y}^{c}", 100, -0.4, 0.4);
  TH1D* h_q2xq2x_a = new TH1D ("h_q2xq2x_a", "(q_{x}^{a})^2", 100, 0, 0.4);
  TH1D* h_q2xq2y_a = new TH1D ("h_q2xq2y_a", "q_{x}^{a}q_y^{a}", 100, -0.4, 0.4);
  TH1D* h_q2yq2y_a = new TH1D ("h_q2yq2y_a", "(q_{x}^{a})^2", 100, 0, 0.4);
  TH1D* h_q2xq2x_c = new TH1D ("h_q2xq2x_c", "(q_{x}^{c})^2", 100, 0, 0.4);
  TH1D* h_q2xq2y_c = new TH1D ("h_q2xq2y_c", "q_{x}^{c}q_y^{c}", 100, -0.4, 0.4);
  TH1D* h_q2yq2y_c = new TH1D ("h_q2yq2y_c", "(q_{x}^{c})^2", 100, 0, 0.4);
  h_q2x_a->Sumw2 ();
  h_q2y_a->Sumw2 ();
  h_q2x_c->Sumw2 ();
  h_q2y_c->Sumw2 ();
  h_q2xq2x_a->Sumw2 ();
  h_q2xq2y_a->Sumw2 ();
  h_q2yq2y_a->Sumw2 ();
  h_q2xq2x_c->Sumw2 ();
  h_q2xq2y_c->Sumw2 ();
  h_q2yq2y_c->Sumw2 ();


  const int nEvts = tree->GetEntries ();

  // First loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In FCalCalibration.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    tree->GetEntry (iEvt);

    if (!isMC) {
      if (collisionSystem == PbPb18 && t->BlayerDesyn)
        continue;
      if (isPbPb && t->isOOTPU)
        continue;
    }

    bool hasPrimaryVert = false;
    float vz = -999;
    for (int iVert = 0; !hasPrimaryVert && iVert < t->nvert; iVert++) {
      hasPrimaryVert = t->vert_type[iVert] == 1;
      vz = t->vert_z[iVert];
    }
    if (!hasPrimaryVert || fabs (vz) > 150)
      continue;

    const float fcal_et = t->fcalA_et + t->fcalC_et;

    if (isPbPb && !isMC) {
      const float zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

    if (isPbPb) {
      h_PbPb_fcal_et->Fill (fcal_et);

      const int centrality = GetPercentileCentrality (fcal_et);
      if (HLT_mb_sptrk_L1ZDC_A_C_VTE50)
        h_centrality_HLT_mb_sptrk_L1ZDC_A_C_VTE50->Fill (centrality);//, HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale);
      if (HLT_noalg_pc_L1TE50_VTE600_0ETA49)
        h_centrality_HLT_noalg_pc_L1TE50_VTE600_0ETA49->Fill (centrality);//, HLT_noalg_pc_L1TE50_VTE600_0ETA49_prescale);
      if (HLT_noalg_cc_L1TE600_0ETA49)
        h_centrality_HLT_noalg_cc_L1TE600_0ETA49->Fill (centrality);//, HLT_noalg_cc_L1TE600_0ETA49_prescale);
      h_q2x_a->Fill (t->fcalA_et_Cos2);
      h_q2y_a->Fill (t->fcalA_et_Sin2);
      h_q2x_c->Fill (t->fcalC_et_Cos2);
      h_q2y_c->Fill (t->fcalC_et_Sin2);
      h_q2xq2x_a->Fill (t->fcalA_et_Cos2 * t->fcalA_et_Cos2);
      h_q2xq2y_a->Fill (t->fcalA_et_Cos2 * t->fcalA_et_Sin2);
      h_q2yq2y_a->Fill (t->fcalA_et_Sin2 * t->fcalA_et_Sin2);
      h_q2xq2x_c->Fill (t->fcalC_et_Cos2 * t->fcalC_et_Cos2);
      h_q2xq2y_c->Fill (t->fcalC_et_Cos2 * t->fcalC_et_Sin2);
      h_q2yq2y_c->Fill (t->fcalC_et_Sin2 * t->fcalC_et_Sin2);
    }
    else
      h_pp_fcal_et->Fill (fcal_et);

  } // end event loop
  cout << endl << "Info: In FCalCalibration.cxx: Finished processing events." << endl;

  h_PbPb_fcal_et->Write ();
  h_pp_fcal_et->Write ();
  h_centrality_HLT_mb_sptrk_L1ZDC_A_C_VTE50->Write ();
  h_centrality_HLT_noalg_pc_L1TE50_VTE600_0ETA49->Write ();
  h_centrality_HLT_noalg_cc_L1TE600_0ETA49->Write ();
  h_q2x_a->Write ();
  h_q2y_a->Write ();
  h_q2x_c->Write ();
  h_q2y_c->Write ();
  h_q2xq2x_a->Write ();
  h_q2xq2y_a->Write ();
  h_q2yq2y_a->Write ();
  h_q2xq2x_c->Write ();
  h_q2xq2y_c->Write ();
  h_q2yq2y_c->Write ();

  outFile->Close ();
  SaferDelete (&outFile);

  SaferDelete (&h_zdcCuts);

  return true;
}

} // end namespace

#endif
