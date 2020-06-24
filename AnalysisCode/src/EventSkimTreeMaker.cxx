#ifndef __EventSkimTreeMaker_cxx__
#define __EventSkimTreeMaker_cxx__

#include "EventSkimTreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"
#include "ZTrackUtilities.h"
#include "OutTree.h"

#include <Utilities.h>
#include <AtlasUtils.h>

#include <TH2D.h>
#include <TChain.h>
#include <TSystem.h>
#include <TRandom3.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool EventSkimTreeMaker (const char* directory,
                       const int dataSet,
                       const char* inFileName) {

  cout << "Info: In EventSkimTreeMaker.cxx: Entered TreeMaker routine." << endl;
  cout << "Info: In EventSkimTreeMaker.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar << endl;

  SetupDirectories ("EventSkims");

  const bool isPCStream = (strstr (directory, "pc") != NULL);
  if (isPCStream)
    cout << "Info: In EventSkimTreeMaker.cxx: detected PC stream data" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In EventSkimTreeMaker.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In EventSkimTreeMaker.cxx: Saving output to " << rootPath << endl;

  TChain* tree = new TChain ("bush", "bush");
  TString pattern = "*.root";
  auto dir = gSystem->OpenDirectory (dataPath + directory);
  while (const char* f = gSystem->GetDirEntry (dir)) {
    TString file = TString (f);
    if (!file.Contains (Form ("%i", dataSet)))
      continue;
    cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << endl;
    tree->Add (dataPath + directory + "/" + file + "/*.root");
    break;
  }

  if (tree == nullptr) {
    cout << "Error: In EventSkimTreeMaker.cxx: TChain not obtained for given data set. Quitting." << endl;
    return false;
  }

  cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << endl;

  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In EventSkimTreeMaker.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }


  // Input tree
  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  if (isPbPb && !isMC) t->SetGetZdc ();
  t->SetBranchAddresses ();

  bool HLT_mb_sptrk = false;
  bool HLT_mb_sptrk_L1ZDC_A_C_VTE50 = false;
  bool HLT_noalg_pc_L1TE50_VTE600_0ETA49 = false;
  bool HLT_noalg_cc_L1TE600_0ETA49 = false;
  if (isPbPb && !isMC) {
    tree->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",       &HLT_mb_sptrk_L1ZDC_A_C_VTE50);
    tree->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49",  &HLT_noalg_pc_L1TE50_VTE600_0ETA49);
    tree->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49",        &HLT_noalg_cc_L1TE600_0ETA49);
  }
  else if (!isPbPb && !isMC) {
    tree->SetBranchAddress ("HLT_mb_sptrk",                       &HLT_mb_sptrk);
  }


  // Load files for output
  TFile* outFiles;
  OutTree* outTrees;
  const char* outTreeName = (isPbPb ? "PbPbTrackTree" : "ppTrackTree");
  if (isPbPb) {
    const TString runGroup = (isMC ? "Hijing" : GetRunGroupTString (dataSet));
    outFiles = new TFile (Form ("%s/%s/%s.root", rootPath.Data (), runGroup.Data (), identifier.Data ()), "recreate");
    outFiles->Delete (Form ("%s;*", outTreeName));
    outTrees = new OutTree (outTreeName, outFiles);
    outTrees->SetBranchEventInfo ();
    outTrees->SetBranches ();
    outTrees->tree->Branch ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",      &HLT_mb_sptrk_L1ZDC_A_C_VTE50,      "HLT_mb_sptrk_L1ZDC_A_C_VTE50/O");
    outTrees->tree->Branch ("HLT_noalg_pc_L1TE50_VTE600.0ETA49", &HLT_noalg_pc_L1TE50_VTE600_0ETA49, "HLT_noalg_pc_L1TE50_VTE600.0ETA49/O");
    outTrees->tree->Branch ("HLT_noalg_cc_L1TE600_0ETA49",       &HLT_noalg_cc_L1TE600_0ETA49,       "HLT_noalg_cc_L1TE600_0ETA49/O");
  }
  else {
    outFiles = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
    outFiles->Delete (Form ("%s;*", outTreeName));
    outTrees = new OutTree (outTreeName, outFiles);
    outTrees->SetBranchEventInfo ();
    outTrees->SetBranches ();
    outTrees->tree->Branch ("HLT_mb_sptrk",  &HLT_mb_sptrk,  "HLT_mb_sptrk/O");
  }


  const int nEvts = tree->GetEntries ();

  // First loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In EventSkimTreeMaker.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    tree->GetEntry (iEvt);

    event_number = t->event_number;
    run_number = t->run_number;

    event_weight = 1;
    if (!isMC) {
      if (collisionSystem == PbPb18 && t->BlayerDesyn)
        continue;
      if (isPbPb && t->isOOTPU)
        continue;
    }

    bool hasPrimaryVert = false;
    for (int iVert = 0; !hasPrimaryVert && iVert < t->nvert; iVert++) {
      hasPrimaryVert = t->vert_type[iVert] == 1;
      vz = t->vert_z[iVert];
    }
    if (!hasPrimaryVert || fabs (vz) > 150)
      continue;

    fcal_et = t->fcalA_et + t->fcalC_et;
    q2x_a = t->fcalA_et_Cos2;
    q2y_a = t->fcalA_et_Sin2;
    q2x_c = t->fcalC_et_Cos2;
    q2y_c = t->fcalC_et_Sin2;
    const double q2x = q2x_a + q2x_c;
    const double q2y = q2y_a + q2y_c;
    const double q3x = t->fcalA_et_Cos3 + t->fcalC_et_Cos3;
    const double q3y = t->fcalA_et_Sin3 + t->fcalC_et_Sin3;
    const double q4x = t->fcalA_et_Cos4 + t->fcalC_et_Cos4;
    const double q4y = t->fcalA_et_Sin4 + t->fcalC_et_Sin4;
    if (fcal_et > 0) {
      q2 = sqrt (q2x*q2x + q2y*q2y) / fcal_et;
      q3 = sqrt (q3x*q3x + q3y*q3y) / fcal_et;
      q4 = sqrt (q4x*q4x + q4y*q4y) / fcal_et;
    }
    else {
      q2 = 0.;
      q3 = 0.;
      q4 = 0.;
    }
    psi2 = atan2 (q2y, q2x) / 2.;
    psi3 = atan2 (q3y, q3x) / 3.;
    psi4 = atan2 (q4y, q4x) / 4.;

    if (isPbPb && !isMC) {
      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

    outTrees->Fill ();
  } // end event loop
  cout << endl << "Info: In EventSkimTreeMaker.cxx: Finished processing events." << endl;

  outFiles->Write (0, TObject::kOverwrite);
  outFiles->Close ();
  SaferDelete (&outFiles);

  SaferDelete (&h_zdcCuts);

  return true;
}

} // end namespace

#endif
