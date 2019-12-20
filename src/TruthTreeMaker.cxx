#ifndef __TruthTreeMaker_cxx__
#define __TruthTreeMaker_cxx__

#include "TruthTreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"
#include "OutTree.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool TruthTreeMaker (const char* directory,
                     const int dataSet,
                     const char* inFileName) {
 
  SetupDirectories ("TruthAnalysis");

  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TreeMaker.cxx: Running over data overlay, will check data conditions" << endl;

  if (!isMC) {
    cout << "Error: In TruthTreeMaker.cxx: Cannot run truth-level analysis over data. Quitting." << endl;
    return false;
  }
  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TruthTreeMaker.cxx: File Identifier: " << identifier << endl;

  TChain* tree = new TChain ("bush", "bush");
  TString pattern = "*.root";
  auto dir = gSystem->OpenDirectory (dataPath + directory);

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In TruthTreeMaker.C: Cannot identify this MC file! Quitting." << endl;
      return false;
    }
  }
  else
    fileIdentifier = inFileName;

  while (const char* f = gSystem->GetDirEntry (dir)) {
    TString file = TString (f);
    if (!file.Contains (fileIdentifier))
      continue;
    cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << endl;
    tree->Add (dataPath + directory + "/" + file + "/*.root");
    break;
  }
  cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << endl;

  if (tree == nullptr) {
    cout << "Error: In TagAndProbe.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }

  //TFile* file = GetFile (directory, dataSet, inFileName);
  //TTree* tree = nullptr;
  //if (file != nullptr) tree = (TTree*)file->Get ("bush");
  //if (file == nullptr || tree == nullptr) {
  //  cout << "Error: In TruthTreeMaker.cxx: TTree not obtained for given data set. Quitting." << endl;
  //  return false;
  //}

  //TH1D* h_zdcCuts = GetZdcCuts ();
  //if (h_zdcCuts == nullptr) {
  //  cout << "Error: In TruthTreeMaker.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
  //  return false;
  //}

  //First sort jets & tracks into many, smaller TTrees.
  //This is where the sorting based on event information (e.g. centrality, Ntrk, jet pT) will go.
  //Event mixing will take place based on these categories so that total memory usage at any point in time is minimized.

  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetTruthElectrons ();
  t->SetGetTruthMuons ();
  t->SetGetTruthTracks ();
  t->SetGetTruthJets ();
  t->SetBranchAddresses ();

  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  const char* outTreeName = isPbPb ? "PbPbZTrackTree" : "ppZTrackTree";
  outFile->Delete (Form ("%s;*", outTreeName));
  OutTree* outTree = new OutTree (outTreeName, outFile);
  outTree->SetBranchEventInfo ();
  outTree->SetBranchLeptons ();
  outTree->SetBranchZs ();
  outTree->SetBranchTracks ();
  outTree->SetBranchTruthJets ();
  outTree->SetBranches ();

  const int nEvts = tree->GetEntries ();
  TLorentzVector l1, l2;

  // First loop over events looking for Z->ee candidates
  isEE = true;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TruthTreeMaker.cxx: Electrons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    //if (collisionSystem == PbPb18 && t->BlayerDesyn)
    //  continue;
    //if (isPbPb && t->isOOTPU)
    //  continue;

    bool hasPrimary = false;
    for (int iVert = 0; !hasPrimary && iVert < t->nvert; iVert++) {
      hasPrimary = t->vert_type[iVert] == 1;
      vz = t->vert_z[iVert];
    }
    if (!hasPrimary)
      continue; // check for primary vertex

    event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency / mcNumberEvents;

    if (isPbPb) {
      fcal_et = t->fcalA_et + t->fcalC_et;
      if (fcal_et < 50)
        continue; // loose cut on "noise" events
      const double qx = t->fcalA_et_Cos + t->fcalC_et_Cos;
      const double qy = t->fcalA_et_Sin + t->fcalC_et_Sin;
      if (fcal_et > 0)
        q2 = sqrt (qx*qx + qy*qy) / fcal_et;
      else
        q2 = 0;
      psi2 = 0.5 * atan2 (qy, qx);

      //zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C;
      //const float nNeutrons = (zdcEnergy) / (2.51);
      //const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3);
      //if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || nNeutrons > h_zdcCuts->GetBinContent (bin))
      //  continue; // Zdc-based in-time pile-up cut
    }

    for (int iE1 = 0; iE1 < t->truth_electron_n; iE1++) {
      if (t->truth_electron_pt[iE1] < electron_pt_cut)
        continue; // basic electron pT cut
      //if (!InEMCal (t->truth_electron_eta[iE1]))
      //  continue; // reject electrons reconstructed outside the EMCal
      if (t->truth_electron_barcode[iE1] > 200000)
        continue; // non-primary electrons (e.g. from G4)

      l1_pt = t->truth_electron_pt[iE1];
      l1_eta = t->truth_electron_eta[iE1];
      l1_phi = t->truth_electron_phi[iE1];
      l1_charge = t->truth_electron_charge[iE1];

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, electron_mass);

      for (int iE2 = 0; iE2 < iE1; iE2++) {
        if (t->truth_electron_pt[iE2] < electron_pt_cut)
          continue; // basic electron pT cut
        //if (!InEMCal (t->truth_electron_eta[iE2]))
        //  continue; // reject electrons reconstructed outside the EMCal
        if (t->truth_electron_barcode[iE2] > 200000)
          continue; // non-primary electrons (e.g. from G4)

        l2_pt = t->truth_electron_pt[iE2];
        l2_eta = t->truth_electron_eta[iE2];
        l2_phi = t->truth_electron_phi[iE2];
        l2_charge = t->truth_electron_charge[iE2];

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, electron_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged electrons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        ntrk = 0;
        for (int iTrk = 0; iTrk < t->truth_trk_n; iTrk++) {
          if (t->truth_trk_pt[iTrk] < trk_pt_cut)
            continue; // track minimum pT

          trk_pt[ntrk] = t->truth_trk_pt[iTrk];
          trk_eta[ntrk] = t->truth_trk_eta[iTrk];
          trk_phi[ntrk] = t->truth_trk_phi[iTrk];
          trk_charge[ntrk] = t->truth_trk_charge[iTrk];
          ntrk++;
        } // end loop over tracks

        truth_jet_n = 0;
        for (int iJet = 0; iJet < t->truth_jet_n; iJet++) {
          if (t->truth_jet_pt[iJet] < jet_pt_cut)
            continue; // jet minimum pT

          truth_jet_pt[truth_jet_n] = t->truth_jet_pt[iJet];
          truth_jet_eta[truth_jet_n] = t->truth_jet_eta[iJet];
          truth_jet_phi[truth_jet_n] = t->truth_jet_phi[iJet];
          truth_jet_e[truth_jet_n] = t->truth_jet_e[iJet];
          truth_jet_n++;
        } // end loop over truth jets

        outTree->Fill ();
        
      } // end iE2 loop
    } // end iE1 loop

  } // end electron selection
  cout << endl << "Info: In TruthTreeMaker.cxx: Finished processing electrons." << endl;


  // Now loop over events looking for Z->mumu candidates
  isEE = false;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TruthTreeMaker.cxx: Muons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    //if (collisionSystem == PbPb18 && t->BlayerDesyn)
    //  continue;
    //if (isPbPb && t->isOOTPU)
    //  continue;
    //if ((!isMC || isOverlayMC) && !t->passes_toroid)
    //  continue; // additional check for muon quality

    bool hasPrimary = false;
    for (int iVert = 0; !hasPrimary && iVert < t->nvert; iVert++) {
      hasPrimary = (t->vert_type[iVert] == 1);
      vz = t->vert_z[iVert];
    }
    if (!hasPrimary)
      continue; // check for primary vertex

    event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency / mcNumberEvents;

    if (isPbPb) {
      fcal_et = t->fcalA_et + t->fcalC_et;
      if (fcal_et < 50)
        continue; // loose cut on "noise" events
      const double qx = t->fcalA_et_Cos + t->fcalC_et_Cos;
      const double qy = t->fcalA_et_Sin + t->fcalC_et_Sin;
      if (fcal_et > 0)
        q2 = sqrt (qx*qx + qy*qy) / fcal_et;
      else
        q2 = 0;
      psi2 = 0.5 * atan2 (qy, qx);

      //zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C;
      //const float nNeutrons = (zdcEnergy) / (2.51);
      //const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3);
      //if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || nNeutrons > h_zdcCuts->GetBinContent (bin))
      //  continue; // Zdc-based in-time pile-up cut
    }

    for (int iM1 = 0; iM1 < t->truth_muon_n; iM1++) {
      if (t->truth_muon_pt[iM1] < muon_pt_cut)
        continue; // basic muon pT cut
      //if (fabs (t->truth_muon_eta[iM1]) > 2.5)
      //  continue; // reject muons reconstructed outside the EMCal
      if (t->truth_muon_barcode[iM1] > 200000)
        continue; // non-primary muons (e.g. from G4)

      l1_pt = t->truth_muon_pt[iM1];
      l1_eta = t->truth_muon_eta[iM1];
      l1_phi = t->truth_muon_phi[iM1];
      l1_charge = t->truth_muon_charge[iM1];

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, muon_mass);

      for (int iM2 = 0; iM2 < iM1; iM2++) {
        if (t->truth_muon_pt[iM2] < muon_pt_cut)
          continue; // basic muon pT cut
        //if (fabs (t->truth_muon_eta[iM2]) > 2.5)
        //  continue; // reject muons reconstructed outside the EMCal
        if (t->truth_muon_barcode[iM2] > 200000)
          continue; // non-primary muons (e.g. from G4)

        l2_pt = t->truth_muon_pt[iM2];
        l2_eta = t->truth_muon_eta[iM2];
        l2_phi = t->truth_muon_phi[iM2];
        l2_charge = t->truth_muon_charge[iM2];

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, muon_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged muons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        ntrk = 0;
        for (int iTrk = 0; iTrk < t->truth_trk_n; iTrk++) {
          if (t->truth_trk_pt[iTrk] < trk_pt_cut)
            continue; // track minimum pT

          trk_pt[ntrk] = t->truth_trk_pt[iTrk];
          trk_eta[ntrk] = t->truth_trk_eta[iTrk];
          trk_phi[ntrk] = t->truth_trk_phi[iTrk];
          trk_charge[ntrk] = t->truth_trk_charge[iTrk];
          ntrk++;
        } // end loop over tracks

        truth_jet_n = 0;
        for (int iJet = 0; iJet < t->truth_jet_n; iJet++) {
          if (t->truth_jet_pt[iJet] < jet_pt_cut)
            continue; // jet minimum pT

          truth_jet_pt[truth_jet_n] = t->truth_jet_pt[iJet];
          truth_jet_eta[truth_jet_n] = t->truth_jet_eta[iJet];
          truth_jet_phi[truth_jet_n] = t->truth_jet_phi[iJet];
          truth_jet_e[truth_jet_n] = t->truth_jet_e[iJet];
          truth_jet_n++;
        } // end loop over truth jets

        outTree->Fill ();
        
      } // end iM2 loop
    } // end iM1 loop

  } // end muon selection
  cout << endl << "Info: In TruthTreeMaker.cxx: Finished processing muons." << endl;

  SaferDelete (t);

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (outFile);

  return true;
}

} // end namespace

#endif
