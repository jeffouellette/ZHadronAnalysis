#include "MinbiasTreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"
#include "ZTrackUtilities.h"

#include <TH2D.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TSystem.h>
#include <TRandom3.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool MinbiasTreeMaker (const char* directory,
                       const int dataSet,
                       const char* inFileName) {

  cout << "Info: In MinbiasTreeMaker.cxx: Entered TreeMaker routine." << endl;
  cout << "Info: In MinbiasTreeMaker.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar << endl;

  SetupDirectories ("MinbiasAnalysis");

  const bool isPCStream = (strstr (directory, "pc") != NULL);
  if (isPCStream)
    cout << "Info: In MinbiasTreeMaker.cxx: detected PC stream data, applying PC-derived prescale" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In MinbiasTreeMaker.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In MinbiasTreeMaker.cxx: Saving output to " << rootPath << endl;

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
    cout << "Error: In MinbiasTreeMaker.cxx: TChain not obtained for given data set. Quitting." << endl;
    return false;
  }

  TH1D* h_zdcCuts = GetZdcCuts ();
  if (h_zdcCuts == nullptr) {
    cout << "Error: In TreeMaker.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
    return false;
  }

  // Incoming branches
  unsigned int lumi_block = 0, event_number = 0;
  bool  passes_toroid = false, isOOTPU = false, BlayerDesyn = false;
  float fcalA_et = 0, fcalC_et = 0, fcalA_et_Cos = 0, fcalC_et_Cos = 0, fcalA_et_Sin = 0, fcalC_et_Sin = 0, ZdcCalibEnergy_A = 0, ZdcCalibEnergy_C = 0;
  int   nvert = 0;
  float vert_x[30];
  float vert_y[30];
  float vert_z[30];
  int   vert_ntrk[30];
  int   vert_type[30];

  int   ntrk = 0;
  float trk_pt[10000];
  float trk_eta[10000];
  float trk_phi[10000];
  float trk_charge[10000];
  bool  trk_idcut[10000];
  float trk_d0[10000];
  float trk_z0[10000];
  float trk_vz[10000];
  float trk_theta[10000];

  // Outgoing branches
  int run_number = dataSet; 
  int   out_ntrk_all = 0, out_ntrk = 0;
  float vz = 0, fcal_et = 0, zdcEnergy = 0, q2 = 0, psi2 = 0, event_weight = 1;
  float out_trk_pt[10000];
  float out_trk_eta[10000];
  float out_trk_phi[10000];

  tree->SetBranchAddress ("event_number", &event_number);
  tree->SetBranchAddress ("lumi_block", &lumi_block);
  if (isPbPb) {
    tree->SetBranchAddress ("passes_toroid",  &passes_toroid);
    tree->SetBranchAddress ("isOOTPU",        &isOOTPU);
    tree->SetBranchAddress ("BlayerDesyn",    &BlayerDesyn);
  } else {
    tree->SetBranchStatus ("passes_toroid", 0);
    tree->SetBranchStatus ("isOOTPU",       0);
    tree->SetBranchStatus ("BlayerDesyn",   0);
  }

  tree->SetBranchAddress ("nvert",          &nvert);
  tree->SetBranchAddress ("vert_x",         vert_x);
  tree->SetBranchAddress ("vert_y",         vert_y);
  tree->SetBranchAddress ("vert_z",         vert_z);
  tree->SetBranchAddress ("vert_ntrk",      vert_ntrk);
  tree->SetBranchAddress ("vert_type",      vert_type);

  if (isPbPb) {
    tree->SetBranchAddress ("fcalA_et",         &fcalA_et);
    tree->SetBranchAddress ("fcalC_et",         &fcalC_et);
    tree->SetBranchAddress ("fcalA_et_Cos",     &fcalA_et_Cos);
    tree->SetBranchAddress ("fcalC_et_Cos",     &fcalC_et_Cos);
    tree->SetBranchAddress ("fcalA_et_Sin",     &fcalA_et_Sin);
    tree->SetBranchAddress ("fcalC_et_Sin",     &fcalC_et_Sin);
    if (!isMC) {
      tree->SetBranchAddress ("ZdcCalibEnergy_A", &ZdcCalibEnergy_A);
      tree->SetBranchAddress ("ZdcCalibEnergy_C", &ZdcCalibEnergy_C);
    }
  } else {
    tree->SetBranchStatus ("fcalA_et",      0); 
    tree->SetBranchStatus ("fcalC_et",      0); 
    tree->SetBranchStatus ("fcalA_et_Cos",  0); 
    tree->SetBranchStatus ("fcalC_et_Cos",  0); 
    tree->SetBranchStatus ("fcalA_et_Sin",  0); 
    tree->SetBranchStatus ("fcalC_et_Sin",  0); 
  }

  tree->SetBranchAddress ("ntrk",         &ntrk);
  tree->SetBranchAddress ("trk_pt",       trk_pt);
  tree->SetBranchAddress ("trk_eta",      trk_eta);
  tree->SetBranchAddress ("trk_phi",      trk_phi);
  tree->SetBranchAddress ("trk_charge",   trk_charge);
  if (doHITightVar) {
    tree->SetBranchAddress ("trk_HItight",  trk_idcut);
    tree->SetBranchStatus  ("trk_HIloose",  0);
  } else {
    tree->SetBranchAddress ("trk_HIloose",  trk_idcut);
    tree->SetBranchStatus  ("trk_HItight",  0);
  }
  tree->SetBranchAddress ("trk_d0",       trk_d0);
  tree->SetBranchAddress ("trk_z0",       trk_z0);
  tree->SetBranchAddress ("trk_vz",       trk_vz);
  tree->SetBranchAddress ("trk_theta",    trk_theta);

  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  const char* outTreeName = isPbPb ? "PbPbZTrackTree" : "ppZTrackTree";
  outFile->Delete (Form ("%s;*", outTreeName));
  TTree* outTree = new TTree (outTreeName, outTreeName);
  outTree->SetDirectory (outFile);
  outTree->Branch ("run_number",    &run_number);
  outTree->Branch ("event_number",  &event_number);
  outTree->Branch ("lumi_block",    &lumi_block);
  outTree->Branch ("passes_toroid", &passes_toroid);
  outTree->Branch ("vz",            &vz);
  if (isPbPb) {
    outTree->Branch ("fcal_et",   &fcal_et);
    outTree->Branch ("zdcEnergy", &zdcEnergy);
    outTree->Branch ("q2",        &q2);
    outTree->Branch ("psi2",      &psi2);
  }
  outTree->Branch ("event_weight", &event_weight);

  outTree->Branch ("ntrk_all", &out_ntrk_all);
  outTree->Branch ("ntrk",     &out_ntrk);
  outTree->Branch ("trk_pt",   &out_trk_pt,  "trk_pt[ntrk]/F");
  outTree->Branch ("trk_eta",  &out_trk_eta, "trk_eta[ntrk]/F");
  outTree->Branch ("trk_phi",  &out_trk_phi, "trk_phi[ntrk]/F");

  const int nEvts = tree->GetEntries ();

  TRandom3* rndm = new TRandom3 ();
  rndm->SetSeed ();

  // First loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In MinbiasTreeMaker.cxx: Events " << iEvt / (nEvts / 100) << "\% done...\r" << flush;

    if (!isMC && is2015data && rndm->Rndm () > 0.3)
      continue;

    tree->GetEntry (iEvt);

    run_number = dataSet; 

    if (!isMC) {
      if (collisionSystem == PbPb18 && BlayerDesyn)
        continue;
      if (isPbPb && isOOTPU)
        continue;
      if (isPCStream)
        event_weight = 0.213632;
      else
        event_weight = 1;
    }
    else {
      event_weight = 1;
    }

    bool hasPrimaryVert = false;
    for (int iVert = 0; !hasPrimaryVert && iVert < nvert; iVert++) {
      hasPrimaryVert = vert_type[iVert] == 1;
      vz = vert_z[iVert];
    }
    if (!hasPrimaryVert)
      continue;

    if (isPbPb) {
      fcal_et = fcalA_et + fcalC_et;
      if (fcal_et < 50)
        continue; // loose cut on "noise" events

      const double qx = fcalA_et_Cos + fcalC_et_Cos;
      const double qy = fcalA_et_Sin + fcalC_et_Sin;
      if (fcal_et > 0)
        q2 = sqrt (qx*qx + qy*qy) / fcal_et;
      else
        q2 = 0;
      psi2 = 0.5 * atan2 (qy, qx);

      if (!isMC) {
        zdcEnergy = ZdcCalibEnergy_A + ZdcCalibEnergy_C; // gets zdc energy in TeV
        const float nNeutrons = (zdcEnergy) / (2.51);
        const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
        if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
          continue; // Zdc-based in-time pile-up cut
      }
    }

    out_ntrk_all = 0;
    out_ntrk = 0;
    for (int iTrk = 0; iTrk < ntrk; iTrk++) {
      if (!trk_idcut[iTrk])
        continue;

      out_ntrk_all++;

      if (trk_pt[iTrk] < trk_pt_cut)
        continue;

      out_trk_pt[out_ntrk] = trk_pt[iTrk];
      out_trk_eta[out_ntrk] = trk_eta[iTrk];
      out_trk_phi[out_ntrk] = trk_phi[iTrk];
      out_ntrk++;
    }

    outTree->Fill ();
  } // end electron selection
  cout << endl << "Info: In MinbiasTreeMaker.cxx: Finished processing events." << endl;

  SaferDelete (rndm);

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (outFile);

  SaferDelete (h_zdcCuts);

  return true;
}

} // end namespace
