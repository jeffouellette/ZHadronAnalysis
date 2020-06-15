#ifndef __TruthSelectedTreeMaker_cxx__
#define __TruthSelectedTreeMaker_cxx__

#include "TruthSelectedTreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"
#include "OutTree.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"
#include "Cuts.h"

#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool TruthSelectedTreeMaker (const char* directory,
                             const int dataSet,
                             const char* inFileName) {
 
  cout << "Info: In TruthSelectedTreeMaker.cxx: Entered TruthSelectedTreeMaker routine." << endl;
  cout << "Info: In TruthSelectedTreeMaker.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar;
  cout << endl;

  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  if (isHijing)
    cout << "Info: In TruthSelectedTreeMaker.cxx: File detected as Hijing overlay" << endl;
  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TruthSelectedTreeMaker.cxx: File detected as data overlay, will check data conditions" << endl;

  if (isMC && isHijing)
    SetupDirectories ("MCAnalysis/Hijing", false);
  else if (isMC)
    SetupDirectories ("MCAnalysis");
  else
    SetupDirectories ("DataAnalysis");
    

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TruthSelectedTreeMaker.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TruthSelectedTreeMaker.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In TruthSelectedTreeMaker.C: Cannot identify this MC file! Quitting." << endl;
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
  cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << endl;

  if (tree == nullptr) {
    cout << "Error: In TruthSelectedTreeMaker.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }


  // Input tree
  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetMuons ();
  t->SetGetTruthElectrons ();
  t->SetGetTruthMuons ();
  t->SetGetTracks ();
  t->SetGetTruthTracks ();
  if (!isMC && isPbPb)
    t->SetGetZdc ();
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


  // Load files for output
  const int numFileBins = (isHijing ? numFileIPBins : numFileCentBins);
  TFile* outFiles[numFileBins];
  OutTree* outTrees[numFileBins];
  const char* outTreeName = (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree");
  if (isPbPb) {
    for (int iCent = 1; iCent < numFileBins; iCent++) {
      if (!isHijing) {
        outFiles[iCent] = new TFile (Form ("%s/%s/%s_iCent%i.root", rootPath.Data (), GetRunGroupTString (dataSet).Data (), identifier.Data (), iCent), "recreate");
        outFiles[iCent]->Delete (Form ("%s;*", outTreeName));
        outTrees[iCent] = new OutTree (outTreeName, outFiles[iCent]);
        outTrees[iCent]->SetBranchEventInfo ();
        outTrees[iCent]->SetBranchLeptons ();
        outTrees[iCent]->SetBranchZs ();
        outTrees[iCent]->SetBranchTracks ();
        outTrees[iCent]->SetBranchTruthTracks ();
        outTrees[iCent]->SetBranches ();
      }
      else {
        outFiles[iCent] = new TFile (Form ("%s/%s/%s_iCent%i.root", rootPath.Data (), identifier.Data (), identifier.Data (), iCent), "recreate");
        outFiles[iCent]->Delete (Form ("%s;*", outTreeName));
        outTrees[iCent] = new OutTree (outTreeName, outFiles[iCent]);
        outTrees[iCent]->SetBranchEventInfo ();
        outTrees[iCent]->SetBranchLeptons ();
        outTrees[iCent]->SetBranchZs ();
        outTrees[iCent]->SetBranchTracks ();
        outTrees[iCent]->SetBranchTruthTracks ();
        outTrees[iCent]->SetBranches ();
        outTrees[iCent]->tree->Branch ("impactParameter", &ip, "impactParameter/F");
        outTrees[iCent]->tree->Branch ("eventPlane", &eventPlane, "eventPlane/F");
      }
    } // end loop over iCent
  }
  else {
    outFiles[0] = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
    outFiles[0]->Delete (Form ("%s;*", outTreeName));
    outTrees[0] = new OutTree (outTreeName, outFiles[0]);
    outTrees[0]->SetBranchEventInfo ();
    outTrees[0]->SetBranchLeptons ();
    outTrees[0]->SetBranchZs ();
    outTrees[0]->SetBranchTracks ();
    outTrees[0]->SetBranchTruthTracks ();
    outTrees[0]->SetBranches ();
  }


  const int nEvts = tree->GetEntries ();
  TLorentzVector l1, l2;

  // First loop over events looking for Z->ee candidates
  isEE = true;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TruthSelectedTreeMaker.cxx: Electrons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    if (!isHijing) {
      if (collisionSystem == PbPb18 && t->BlayerDesyn)
        continue;
      if (isPbPb && t->isOOTPU)
        continue;
    }

    event_number = t->event_number;
    run_number = t->run_number;

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

    event_weight = crossSectionPicoBarns * mcFilterEfficiency * (isPbPb ? 0.00171786 : 258.4); // sigma * f * L_int
    //event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency * 0.00171786; // sigma * f * L_int
    if (isHijing) {
      assert (t->nTruthEvt > 0);
      ip = t->impactParameter[0];
      eventPlane = t->eventPlaneAngle[0];
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

    if (isPbPb)
      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV

    short iCent = 0;
    if (isPbPb && !isHijing) {
      iCent = GetFileCentBin (fcal_et);
      if (iCent < 1 || iCent > numFileBins-1)
        continue;
    }
    else if (isPbPb && isHijing) {
      iCent = GetFileIPBin (ip);
      if (iCent < 1 || iCent > numFileBins-1)
        continue;
    }

    for (int iE1 = 0; iE1 < t->truth_electron_n; iE1++) {
      if (t->truth_electron_pt[iE1] < electron_pt_cut)
        continue; // basic electron pT cut
      if (!InEMCal (t->truth_electron_eta[iE1]))
        continue; // reject electrons reconstructed outside the EMCal
      if (t->truth_electron_barcode[iE1] >= 200000)
        continue; // non-primary electrons (e.g. from G4)

      l1_pt = t->truth_electron_pt[iE1];
      l1_eta = t->truth_electron_eta[iE1];
      l1_phi = t->truth_electron_phi[iE1];
      l1_charge = t->truth_electron_charge[iE1];

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, electron_mass);

      int iRE1 = -1;
      for (int iRE = 0; iRE < t->electron_n; iRE++) {
        if (DeltaR (l1_eta, t->electron_eta[iRE], l1_phi, t->electron_phi[iRE]) < 0.08) {
          iRE1 = iRE;
          break;
        }
      }

      for (int iE2 = 0; iE2 < iE1; iE2++) {
        if (t->truth_electron_pt[iE2] < electron_pt_cut)
          continue; // basic electron pT cut
        if (!InEMCal (t->truth_electron_eta[iE2]))
          continue; // reject electrons reconstructed outside the EMCal
        if (t->truth_electron_barcode[iE2] >= 200000)
          continue; // non-primary electrons (e.g. from G4)

        l2_pt = t->truth_electron_pt[iE2];
        l2_eta = t->truth_electron_eta[iE2];
        l2_phi = t->truth_electron_phi[iE2];
        l2_charge = t->truth_electron_charge[iE2];

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, electron_mass);

        int iRE2 = -1;
        for (int iRE = 0; iRE < t->electron_n; iRE++) {
          if (DeltaR (l2_eta, t->electron_eta[iRE], l2_phi, t->electron_phi[iRE]) < 0.08) {
            iRE2 = iRE;
            break;
          }
        }

        const float nom_z_pt = (l1+l2).Pt ();
        const float nom_z_y = (l1+l2).Rapidity ();

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged electrons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        if (isMC && !isHijing) {
          if (isPbPb) {
            if (nom_z_pt < 25 && fabs (nom_z_y) < 1.75)       event_weight = event_weight / (415887+638346+638321+976912);
            else if (nom_z_pt >= 25 && fabs (nom_z_y) < 1.75) event_weight = event_weight / (415887+638346+638321+976912+415862+638366+638366+976816);
            else if (nom_z_pt < 25 && fabs (nom_z_y) >= 1.75) event_weight = event_weight / (415887+638346+638321+976912+415862+638366+637401+976816);
            else {
              assert (nom_z_pt >= 25 && fabs (nom_z_y) >= 1.75);
              event_weight = event_weight / (415887+638346+638321+976912+415862+638366+638366+976816+415862+638366+637401+976816);
            }
          }
          else event_weight = event_weight / 1000000;
        }

        ntrk_perp = 0;
        ntrk = 0;

        float sumptcw = 0;
        float sumptccw = 0;
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

          if (IsElectronTrack (t, iTrk, iRE1, iRE2))
            continue;

          float dphi = DeltaPhi (t->trk_phi[ntrk], z_phi, true);
          if (pi/3. < fabs (dphi) && fabs (dphi) < 2.*pi/3.) {
            if (dphi < 0) // then trackphi is ccw of zphi
              sumptccw += t->trk_pt[iTrk];
            else
              sumptcw += t->trk_pt[iTrk];
          }
          dphi = DeltaPhi (t->trk_phi[ntrk], z_phi, false);
          if (pi/4. < dphi && dphi < 3.*pi/4.)
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
        } // end loop over tracks

        truth_ntrk = 0;
        for (int iTrk = 0; iTrk < t->truth_trk_n; iTrk++) {
          if (t->truth_trk_pt[iTrk] < trk_pt_cut)
            continue; // track minimum pT
          if (fabs (t->truth_trk_eta[iTrk]) > 2.5)
            continue;

          if (t->truth_trk_barcode[iTrk] <= 0 || 200000 <= t->truth_trk_barcode[iTrk] || abs (t->truth_trk_pdgid[iTrk]) == 3112 || abs (t->truth_trk_pdgid[iTrk]) == 3222 || abs (t->truth_trk_pdgid[iTrk]) == 3312 || abs (t->truth_trk_pdgid[iTrk]) == 3334)
            continue; // select primary truth particles
    
          //if (fabs (t->truth_trk_pdgid[iTrk]) == 11 || fabs (t->truth_trk_pdgid[iTrk]) == 13)
          //  continue;
          if (!(t->truth_trk_isHadron[iTrk]))
            continue;

          truth_trk_pt[truth_ntrk] = t->truth_trk_pt[iTrk];
          truth_trk_eta[truth_ntrk] = t->truth_trk_eta[iTrk];
          truth_trk_phi[truth_ntrk] = t->truth_trk_phi[iTrk];
          truth_trk_charge[truth_ntrk] = t->truth_trk_charge[iTrk];
          truth_ntrk++;
        } // end loop over tracks

        if (sumptccw > sumptcw) {
          phi_transmax = z_phi + pi/2.;
          phi_transmin = z_phi - pi/2.;
        }
        else {
          phi_transmax = z_phi - pi/2.;
          phi_transmin = z_phi + pi/2.;
        }

        ntrk_perp *= 2;
        //if (isPbPb && isHijing) {
        //  iCent = GetFileNtrkBin (ntrk_perp);
        //  if (isPbPb && (iCent < 1 || iCent > numFileBins-1))
        //    continue;
        //}

        outTrees[iCent]->Fill ();
        
      } // end iE2 loop
    } // end iE1 loop

  } // end electron selection
  cout << endl << "Info: In TruthSelectedTreeMaker.cxx: Finished processing electrons." << endl;


  // Now loop over events looking for Z->mumu candidates
  isEE = false;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TruthSelectedTreeMaker.cxx: Muons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    if (!isHijing) {
      if (collisionSystem == PbPb18 && t->BlayerDesyn)
        continue;
      if (isPbPb && t->isOOTPU)
        continue;
      if ((!isMC || isOverlayMC) && !t->passes_toroid)
        continue; // additional check for muon quality
    }

    event_number = t->event_number;
    run_number = t->run_number;

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

    event_weight = crossSectionPicoBarns * mcFilterEfficiency * (isPbPb ? 0.00143844 : 258.4); // sigma * f * L_int
    //event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency * (isPbPb ? 0.00143844 : 258.4); // sigma * f * L_int
    if (isHijing) {
      assert (t->nTruthEvt > 0);
      ip = t->impactParameter[0];
      eventPlane = t->eventPlaneAngle[0];
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

    if (isPbPb)
      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV

    short iCent = 0;
    if (isPbPb && !isHijing) {
      iCent = GetFileCentBin (fcal_et);
      if (iCent < 1 || iCent > numFileBins-1)
        continue;
    }
    else if (isPbPb && isHijing) {
      iCent = GetFileIPBin (ip);
      if (iCent < 1 || iCent > numFileBins-1)
        continue;
    }

    for (int iM1 = 0; iM1 < t->truth_muon_n; iM1++) {
      if (t->truth_muon_pt[iM1] < muon_pt_cut)
        continue; // basic muon pT cut
      if (fabs (t->truth_muon_eta[iM1]) > 2.5)
        continue; // reject muons reconstructed outside the EMCal
      if (t->truth_muon_barcode[iM1] >= 200000)
        continue; // non-primary muons (e.g. from G4)

      l1_pt = t->truth_muon_pt[iM1];
      l1_eta = t->truth_muon_eta[iM1];
      l1_phi = t->truth_muon_phi[iM1];
      l1_charge = t->truth_muon_charge[iM1];

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, muon_mass);

      int iRM1 = -1;
      for (int iRM = 0; iRM < t->muon_n; iRM++) {
        if (DeltaR (l1_eta, t->muon_eta[iRM], l1_phi, t->muon_phi[iRM]) < 0.08) {
          iRM1 = iRM;
          break;
        }
      }

      for (int iM2 = 0; iM2 < iM1; iM2++) {
        if (t->truth_muon_pt[iM2] < muon_pt_cut)
          continue; // basic muon pT cut
        if (fabs (t->truth_muon_eta[iM2]) > 2.5)
          continue; // reject muons reconstructed outside the EMCal
        if (t->truth_muon_barcode[iM2] >= 200000)
          continue; // non-primary muons (e.g. from G4)

        l2_pt = t->truth_muon_pt[iM2];
        l2_eta = t->truth_muon_eta[iM2];
        l2_phi = t->truth_muon_phi[iM2];
        l2_charge = t->truth_muon_charge[iM2];

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, muon_mass);

        int iRM2 = -1;
        for (int iRM = 0; iRM < t->muon_n; iRM++) {
          if (DeltaR (l2_eta, t->muon_eta[iRM], l2_phi, t->muon_phi[iRM]) < 0.08) {
            iRM2 = iRM;
            break;
          }
        }

        const float nom_z_pt = (l1+l2).Pt ();
        const float nom_z_y = (l1+l2).Rapidity ();

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged muons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        if (isMC && !isHijing) {
          if (isPbPb) {
            if (nom_z_pt < 25 && fabs (nom_z_y) < 1.75)       event_weight = event_weight / (357861+541661+561006+840470);
            else if (nom_z_pt >= 25 && fabs (nom_z_y) < 1.75) event_weight = event_weight / (357861+541661+561006+840470+357894+560953+560953+841487);
            else if (nom_z_pt < 25 && fabs (nom_z_y) >= 1.75) event_weight = event_weight / (357861+541661+561006+840470+357900+560051+561091+839951);
            else {
              assert (nom_z_pt >= 25 && fabs (nom_z_y) >= 1.75);
              event_weight = event_weight / (357861+541661+561006+840470+357894+560953+560953+841487+357900+560051+561091+839951);
            }
          }
          else event_weight = event_weight / 3360000;
        }

        ntrk_perp = 0;
        ntrk = 0;

        float sumptcw = 0;
        float sumptccw = 0;
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

          if (IsMuonTrack (t, iTrk, iRM1, iRM2))
            continue;

          float dphi = DeltaPhi (t->trk_phi[ntrk], z_phi, true);
          if (pi/3. < fabs (dphi) && fabs (dphi) < 2.*pi/3.) {
            if (dphi < 0) // then trackphi is ccw of zphi
              sumptccw += t->trk_pt[iTrk];
            else
              sumptcw += t->trk_pt[iTrk];
          }
          dphi = DeltaPhi (t->trk_phi[ntrk], z_phi, false);
          if (pi/4. < dphi && dphi < 3.*pi/4.)
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
        } // end loop over tracks

        truth_ntrk = 0;
        for (int iTrk = 0; iTrk < t->truth_trk_n; iTrk++) {
          if (t->truth_trk_pt[iTrk] < trk_pt_cut)
            continue; // track minimum pT
          if (fabs (t->truth_trk_eta[iTrk]) > 2.5)
            continue;

          if (t->truth_trk_barcode[iTrk] <= 0 || 200000 <= t->truth_trk_barcode[iTrk] || abs (t->truth_trk_pdgid[iTrk]) == 3112 || abs (t->truth_trk_pdgid[iTrk]) == 3222 || abs (t->truth_trk_pdgid[iTrk]) == 3312 || abs (t->truth_trk_pdgid[iTrk]) == 3334)
            continue; // select primary truth particles
    
          //if (fabs (t->truth_trk_pdgid[iTrk]) == 11 || fabs (t->truth_trk_pdgid[iTrk]) == 13)
          //  continue;
          if (!(t->truth_trk_isHadron[iTrk]))
            continue;

          truth_trk_pt[truth_ntrk] = t->truth_trk_pt[iTrk];
          truth_trk_eta[truth_ntrk] = t->truth_trk_eta[iTrk];
          truth_trk_phi[truth_ntrk] = t->truth_trk_phi[iTrk];
          truth_trk_charge[truth_ntrk] = t->truth_trk_charge[iTrk];
          truth_ntrk++;
        } // end loop over tracks

        if (sumptccw > sumptcw) {
          phi_transmax = z_phi + pi/2.;
          phi_transmin = z_phi - pi/2.;
        }
        else {
          phi_transmax = z_phi - pi/2.;
          phi_transmin = z_phi + pi/2.;
        }

        ntrk_perp *= 2;
        //if (isPbPb && isHijing) {
        //  iCent = GetFileNtrkBin (ntrk_perp);
        //  if (isPbPb && (iCent < 1 || iCent > numFileBins-1))
        //    continue;
        //}

        outTrees[iCent]->Fill ();
        
      } // end iM2 loop
    } // end iM1 loop

  } // end muon selection
  cout << endl << "Info: In TruthSelectedTreeMaker.cxx: Finished processing muons." << endl;

  SaferDelete (&t);

  if (isPbPb) {
    for (int iCent = 1; iCent < numFileBins; iCent++) {
      outFiles[iCent]->Write (0, TObject::kOverwrite);
      outFiles[iCent]->Close ();
      SaferDelete (&(outFiles[iCent]));
    }
  }
  else {
    outFiles[0]->Write ();//0, TObject::kOverwrite);
    outFiles[0]->Close ();
    SaferDelete (&(outFiles[0]));
  }

  return true;
}

} // end namespace

#endif
