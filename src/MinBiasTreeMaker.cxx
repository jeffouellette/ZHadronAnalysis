#ifndef __MinBiasTreeMaker_cxx__
#define __MinBiasTreeMaker_cxx__

#include "MinBiasTreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"
#include "ZTrackUtilities.h"
#include "OutTree.h"

#include <TH2D.h>
#include <TChain.h>
#include <TSystem.h>
#include <TRandom3.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool MinBiasTreeMaker (const char* directory,
                       const int dataSet,
                       const char* inFileName) {

  cout << "Info: In MinBiasTreeMaker.cxx: Entered TreeMaker routine." << endl;
  cout << "Info: In MinBiasTreeMaker.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar << endl;

  const bool isPCStream = (strstr (directory, "pc") != NULL);
  if (isPCStream)
    cout << "Info: In MinBiasTreeMaker.cxx: detected PC stream data, applying PC-derived prescale" << endl;
  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  if (isHijing)
    cout << "Info: In TreeMaker.cxx: File detected as Hijing overlay" << endl;

  if (isMC && isHijing)
    SetupDirectories ("MixingAnalysis/Hijing", false);
  else
    SetupDirectories ("MixingAnalysis");

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In MinBiasTreeMaker.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In MinBiasTreeMaker.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In TreeMaker.C: Cannot identify this MC file! Quitting." << endl;
      return false;
    }
  }
  else
    fileIdentifier = inFileName;

  TChain* tree = new TChain ("bush", "bush");
  TString pattern = "*.root";
  auto dir = gSystem->OpenDirectory (dataPath + directory);
  while (const char* f = gSystem->GetDirEntry (dir)) {
    TString file = TString (f);
    if (!file.Contains (fileIdentifier))
      continue;
    cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << endl;
    tree->Add (dataPath + directory + "/" + file + "/*.root");
    break;
  }

  if (tree == nullptr) {
    cout << "Error: In MinBiasTreeMaker.cxx: TChain not obtained for given data set. Quitting." << endl;
    return false;
  }

  cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << endl;

  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In MinBiasTreeMaker.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }


  // Input tree
  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetTracks ();
  if (isPbPb && !isMC) t->SetGetZdc ();
  t->SetBranchAddresses ();

  if (isHijing) {
    tree->SetBranchAddress ("nTruthEvt",          &(t->nTruthEvt));
    tree->SetBranchAddress ("nPart1",             t->nPart1);
    tree->SetBranchAddress ("nPart2",             t->nPart2);
    tree->SetBranchAddress ("impactParameter",    t->impactParameter);
    tree->SetBranchAddress ("nColl",              t->nColl);
    tree->SetBranchAddress ("nSpectatorNeutrons", t->nSpectatorNeutrons);
    tree->SetBranchAddress ("nSpectatorProtons",  t->nSpectatorProtons);
    tree->SetBranchAddress ("eccentricity",       t->eccentricity);
    tree->SetBranchAddress ("eventPlaneAngle",    t->eventPlaneAngle);
  }

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
  const int numFileBins = (isHijing ? numFileIPBins : numFileCentBins);
  TFile* outFiles[numFileBins];
  OutTree* outTrees[numFileBins];
  const char* outTreeName = (isPbPb ? "PbPbTrackTree" : "ppTrackTree");
  if (isPbPb) {
    const TString runGroup = (isMC ? identifier : GetRunGroupTString (dataSet));
    for (int iCent = 1; iCent < numFileBins; iCent++) {
      outFiles[iCent] = new TFile (Form ("%s/%s/%s_iCent%i.root", rootPath.Data (), runGroup.Data (), identifier.Data (), iCent), "recreate");
      outFiles[iCent]->Delete (Form ("%s;*", outTreeName));
      outTrees[iCent] = new OutTree (outTreeName, outFiles[iCent]);
      outTrees[iCent]->SetBranchEventInfo ();
      outTrees[iCent]->SetBranchTracks ();
      outTrees[iCent]->SetBranches ();
      outTrees[iCent]->tree->Branch ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",      &HLT_mb_sptrk_L1ZDC_A_C_VTE50,      "HLT_mb_sptrk_L1ZDC_A_C_VTE50/O");
      outTrees[iCent]->tree->Branch ("HLT_noalg_pc_L1TE50_VTE600.0ETA49", &HLT_noalg_pc_L1TE50_VTE600_0ETA49, "HLT_noalg_pc_L1TE50_VTE600.0ETA49/O");
      outTrees[iCent]->tree->Branch ("HLT_noalg_cc_L1TE600_0ETA49",       &HLT_noalg_cc_L1TE600_0ETA49,       "HLT_noalg_cc_L1TE600_0ETA49/O");
      if (isHijing) {
        outTrees[iCent]->tree->Branch ("impactParameter", &ip, "impactParameter/F");
        outTrees[iCent]->tree->Branch ("eventPlane", &eventPlane, "eventPlane/F");
      }
    }
  }
  else {
    outFiles[0] = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
    outFiles[0]->Delete (Form ("%s;*", outTreeName));
    outTrees[0] = new OutTree (outTreeName, outFiles[0]);
    outTrees[0]->SetBranchEventInfo ();
    outTrees[0]->SetBranchTracks ();
    outTrees[0]->SetBranches ();
    outTrees[0]->tree->Branch ("HLT_mb_sptrk",  &HLT_mb_sptrk,  "HLT_mb_sptrk/O");
  }


  const int nEvts = tree->GetEntries ();

  TRandom3* rndm = new TRandom3 ();
  rndm->SetSeed ();

  // First loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In MinBiasTreeMaker.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    if (!isMC && is2015data && rndm->Rndm () > 0.3)
      continue;

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

    if (!isHijing) {
      bool hasPrimary = false;
      bool hasPileup = false;
      vz = -999;
      for (int iVert = 0; iVert < t->nvert; iVert++) {
        const bool isPrimary = (t->vert_type[iVert] == 1);
        hasPrimary = hasPrimary || isPrimary;
        hasPileup = hasPileup || (t->vert_type[iVert] == 3);
        if (isPrimary)
          vz = t->vert_z[iVert];
      }
      if (!hasPrimary || (!isPbPb && hasPileup) || fabs (vz) > 150)
        continue;
    }
    else {
      if (t->nvert > 0) vz = t->vert_z[0];
      else vz = 0;
    }

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
    else if (isHijing) {
      assert (t->nTruthEvt > 0);
      ip = t->impactParameter[0];
      eventPlane = t->eventPlaneAngle[0];
    }

    short iCent = 0;
    if (isPbPb && !isHijing) {
      iCent = GetFileCentBin (fcal_et);
      if (!isMC && !is2015data && !isPCStream && fcal_et < 1378.92 && rndm->Rndm () > 0.48647)
        continue;
      if (iCent < 1 || iCent > numFileBins-1)
        continue;
    }
    else if (isPbPb && isHijing) {
      iCent = GetFileIPBin (ip);
      if (isPbPb && (iCent < 1 || iCent > numFileBins-1))
        continue;
    }

    ntrk_perp = 0;
    ntrk = 0;
    for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      if (doHITightVar && !t->trk_HItight[iTrk])
        continue;
      else if (!doHITightVar && !t->trk_HIloose[iTrk])
        continue;

      if (isPbPb) {
        if (fabs (t->trk_d0sig[iTrk]) > 3.0)
          continue; // d0 significance cut in PbPb
        if (fabs (t->trk_z0sig[iTrk]) > 3.0)
          continue; //z0 significance cut in PbPb
      }

      if (t->trk_pt[iTrk] < trk_pt_cut)
        continue; // track minimum pT
      if (fabs (t->trk_eta[iTrk]) > 2.5)
        continue; // track maximum eta

      ntrk_perp++;

      trk_pt[ntrk] = t->trk_pt[iTrk];
      trk_eta[ntrk] = t->trk_eta[iTrk];
      trk_phi[ntrk] = t->trk_phi[iTrk];
      trk_charge[ntrk] = t->trk_charge[iTrk];
      trk_d0[ntrk] = t->trk_d0[iTrk];
      trk_z0[ntrk] = t->trk_z0[iTrk];

      if (isMC)
        trk_truth_matched[ntrk] = (t->trk_prob_truth[iTrk] > 0.5);

      ntrk++;
    }

    //if (isPbPb && isHijing) {
    //  iCent = GetFileNtrkBin (ntrk_perp);
    //  if (isPbPb && (iCent < 1 || iCent > numFileBins-1))
    //    continue;
    //}

    outTrees[iCent]->Fill ();
  } // end event loop
  cout << endl << "Info: In MinBiasTreeMaker.cxx: Finished processing events." << endl;

  SaferDelete (&rndm);

  if (isPbPb) {
    for (int iCent = 1; iCent < numFileBins; iCent++) {
      outFiles[iCent]->Write (0, TObject::kOverwrite);
      outFiles[iCent]->Close ();
      SaferDelete (&(outFiles[iCent]));
    }
  } else {
    outFiles[0]->Write (0, TObject::kOverwrite);
    outFiles[0]->Close ();
    SaferDelete (&(outFiles[0]));
  }

  SaferDelete (&h_zdcCuts);

  return true;
}

} // end namespace

#endif
