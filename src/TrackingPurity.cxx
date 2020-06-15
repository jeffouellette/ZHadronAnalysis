#ifndef __TrackingPurity_cxx__
#define __TrackingPurity_cxx__

#include "TrackingPurity.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool TrackingPurity (const char* directory,
                     const int dataSet,
                     const char* inFileName,
                     const char* eventWeightsFileName) {
 
  cout << "Info: In TrackingPurity.cxx: Entered TrackingPurity routine." << endl;
  cout << "Info: In TrackingPurity.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar << endl;

  SetupDirectories ("TrackingPurities");

  if (!isMC) {
    cout << "Error: In TrackingPurity.cxx: Trying to calculate tracking efficiency in data! Quitting." << endl;
    return false;
  }

  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TrackingPurity.cxx: Running over data overlay, will check data conditions" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TrackingPurity.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TrackingPurity.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    cout << "Error: In TrackingPurity.C: Cannot identify this MC file! Quitting." << endl;
    return false;
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
    cout << "Error: In TrackingPurity.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }

  TFile* eventWeightsFile = nullptr;
  TH1D* h_weights = nullptr;

  //if (isHijing) {
    eventWeightsFile = new TFile (eventWeightsFileName, "read");
    h_weights = (TH1D*) eventWeightsFile->Get (Form ("h_PbPb%s_weights_%s", doNchWeighting ? "Nch" : "FCal", isHijing ? "hijing" : "mc"));
  //}
  //else {
  //  eventWeightsFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/pythiaNchWeights.root", "read");
  //  h_weights = (TH1D*) eventWeightsFile->Get ("h_PythiaNchWeights");
  //}

  //First sort jets & tracks into many, smaller TTrees.
  //This is where the sorting based on event information (e.g. centrality, Ntrk, jet pT) will go.
  //Event mixing will take place based on these categories so that total memory usage at any point in time is minimized.

  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetTracks ();
  //t->SetGetElectrons ();
  //t->SetGetMuons ();
  t->SetGetTruthElectrons ();
  t->SetGetTruthMuons ();
  //t->SetGetTruthTracks ();
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

  cout << "Info : In TrackingPurity.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  const double centBins[4] = {66.402, 1378.92, 2995.94, 5000};
  ////const int centCuts[4] = {80, 30, 10, 0};
  //const double centBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
  //const int centCuts[10] = {80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  const int numCentBins = sizeof (centBins) / sizeof (centBins[0]);

  const double ipBins[4] = {14.031, 8.582, 4.952, 0};
  const int numIPBins = sizeof (ipBins) / sizeof (ipBins[0]);

  const int numFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, numFinerEtaTrkBins);

  const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  const int numEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

  const double ptTrkBinsPP[26] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 60};
  const int numPtTrkBinsPP = sizeof (ptTrkBinsPP) / sizeof (ptTrkBinsPP[0]) - 1;
  const double ptTrkBinsPbPb[25] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 15, 60};
  const int numPtTrkBinsPbPb = sizeof (ptTrkBinsPbPb) / sizeof (ptTrkBinsPbPb[0]) - 1;

  //const int numPtTrkBins = 20;
  //const double* ptTrkBins = logspace (0.5, 60, numPtTrkBins);

  //TEfficiency*** h_trk_pur = new TEfficiency**[numCentBins];
  TH2D** h2_primary_reco_tracks = new TH2D*[numCentBins];
  TH2D** h2_reco_tracks = new TH2D*[numCentBins];
  TH1D*** h_primary_reco_tracks = new TH1D**[numCentBins];
  TH1D*** h_reco_tracks = new TH1D**[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h2_primary_reco_tracks[iCent] = new TH2D (Form ("h2_primary_reco_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, (iCent == 0 ? numPtTrkBinsPP : numPtTrkBinsPbPb), (iCent == 0 ? ptTrkBinsPP : ptTrkBinsPbPb));
    h2_primary_reco_tracks[iCent]->Sumw2 ();

    h2_reco_tracks[iCent] = new TH2D (Form ("h2_reco_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, (iCent == 0 ? numPtTrkBinsPP : numPtTrkBinsPbPb), (iCent == 0 ? ptTrkBinsPP : ptTrkBinsPbPb));
    h2_reco_tracks[iCent]->Sumw2 ();

    h_primary_reco_tracks[iCent] = new TH1D*[numEtaTrkBins];
    h_reco_tracks[iCent] = new TH1D*[numEtaTrkBins];

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h_primary_reco_tracks[iCent][iEta] = new TH1D (Form ("h_primary_reco_tracks_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T} [GeV]", (iCent == 0 ? numPtTrkBinsPP : numPtTrkBinsPbPb), (iCent == 0 ? ptTrkBinsPP : ptTrkBinsPbPb));
      h_primary_reco_tracks[iCent][iEta]->Sumw2 ();

      h_reco_tracks[iCent][iEta] = new TH1D (Form ("h_reco_tracks_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T} [GeV]", (iCent == 0 ? numPtTrkBinsPP : numPtTrkBinsPbPb), (iCent == 0 ? ptTrkBinsPP : ptTrkBinsPbPb));
      h_reco_tracks[iCent][iEta]->Sumw2 ();
    }
  }

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TrackingPurity.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    if (isOverlayMC) {
      if (collisionSystem == PbPb18 && (t->BlayerDesyn || !t->passes_toroid))
        continue;
      if (isPbPb && t->isOOTPU)
        continue;
    }

    bool hasPrimary = false;
    bool hasPileup = false;
    float vz = -999;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      const bool isPrimary = (t->vert_type[iVert] == 1);
      hasPrimary = hasPrimary || isPrimary;
      hasPileup = hasPileup || (t->vert_type[iVert] == 3);
      if (isPrimary)
        vz = t->vert_z[iVert];
    }
    if (!hasPrimary || hasPileup || fabs (vz) > 150)
      continue;

    short iCent = 0;
    if (isPbPb && !isHijing) { // for data overlay centrality is defined via FCal Et
      const float fcal_et = t->fcalA_et + t->fcalC_et;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;
    }
    else if (isPbPb && isHijing) { // for Hijing centrality is defined via impact parameter
      assert (t->nTruthEvt > 0);
      const float ip = t->impactParameter[0];
      while (iCent < numIPBins) { 
        if (ip > ipBins[iCent])
          break;
        iCent++;
      }
      if (iCent < 1 || iCent > numIPBins-1)
        continue;
    }

    /*
    int iL1 = -1, iL2 = -1;
    bool hasTruthZee = false, hasTruthZmumu = false;

    if (t->truth_electron_n > 2) {
      for (int iE1 = 0; iE1 < t->truth_electron_n; iE1++) {
        if (t->truth_electron_pt[iE1] < electron_pt_cut)
          continue; // basic electron pT cut
        if (!InEMCal (t->truth_electron_eta[iE1]))
          continue; // reject electrons reconstructed outside the EMCal
        if (t->truth_electron_barcode[iE1] >= 200000)
          continue; // non-primary electrons (e.g. from G4)

        const float l1_pt = t->truth_electron_pt[iE1];
        const float l1_eta = t->truth_electron_eta[iE1];
        const float l1_phi = t->truth_electron_phi[iE1];
        const int l1_charge = t->truth_electron_charge[iE1];

        TLorentzVector l1;
        l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, electron_mass);

        for (int iE2 = 0; iE2 < iE1; iE2++) {
          if (t->truth_electron_pt[iE2] < electron_pt_cut)
            continue; // basic electron pT cut
          if (!InEMCal (t->truth_electron_eta[iE2]))
            continue; // reject electrons reconstructed outside the EMCal
          if (t->truth_electron_barcode[iE2] >= 200000)
            continue; // non-primary electrons (e.g. from G4)

          const float l2_pt = t->truth_electron_pt[iE2];
          const float l2_eta = t->truth_electron_eta[iE2];
          const float l2_phi = t->truth_electron_phi[iE2];
          const int l2_charge = t->truth_electron_charge[iE2];

          TLorentzVector l2;
          l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, electron_mass);

          const float z_m = (l1+l2).M ();

          if (l1_charge == l2_charge) 
            continue; // require oppositely charged electrons
          if (z_m < 76 || 106 < z_m)
            continue; // require Z to be in mass window

          hasTruthZee = true;
          iL1 = iE1;
          iL2 = iE2;
        } // end loop over iE2
      } // end loop over iE1
    }

    if (t->truth_muon_n > 2) {
      for (int iM1 = 0; iM1 < t->truth_muon_n; iM1++) {
        if (t->truth_muon_pt[iM1] < muon_pt_cut)
          continue; // basic muon pT cut
        if (fabs (t->truth_muon_eta[iM1]) > 2.5)
          continue; // reject muons reconstructed outside the EMCal
        if (t->truth_muon_barcode[iM1] >= 200000)
          continue; // non-primary muons (e.g. from G4)

        const float l1_pt = t->truth_muon_pt[iM1];
        const float l1_eta = t->truth_muon_eta[iM1];
        const float l1_phi = t->truth_muon_phi[iM1];
        const int l1_charge = t->truth_muon_charge[iM1];

        TLorentzVector l1;
        l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, muon_mass);

        for (int iM2 = 0; iM2 < iM1; iM2++) {
          if (t->truth_muon_pt[iM2] < muon_pt_cut)
            continue; // basic muon pT cut
          if (fabs (t->truth_muon_eta[iM2]) > 2.5)
            continue; // reject muons reconstructed outside the EMCal
          if (t->truth_muon_barcode[iM2] >= 200000)
            continue; // non-primary muons (e.g. from G4)

          const float l2_pt = t->truth_muon_pt[iM2];
          const float l2_eta = t->truth_muon_eta[iM2];
          const float l2_phi = t->truth_muon_phi[iM2];
          const int l2_charge = t->truth_muon_charge[iM2];

          TLorentzVector l2;
          l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, muon_mass);

          const float z_m = (l1+l2).M ();

          if (l1_charge == l2_charge) 
            continue; // require oppositely charged muons
          if (z_m < 76 || 106 < z_m)
            continue; // require Z to be in mass window

          hasTruthZmumu = true;
          iL1 = iM1;
          iL2 = iM2;
        } // end loop over iM2
      } // end loop over iM1
    }
    */

    const float eventWeight = 1.;
    //const float eventWeight =  crossSectionPicoBarns * (isPbPb ? 1. : 258.4) * mcFilterEfficiency / mcNumberEvents; // sigma * f * L_int
    //const float eventWeight = ((isPbPb && !isHijing) ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? t->ntrk : t->fcalA_et+t->fcalC_et)) : 1); // weight is 1 for pp
    //const float eventWeight = (isPbPb ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? t->ntrk : t->fcalA_et+t->fcalC_et)) : 1); // weight is 1 for pp

    for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      if (doHITightVar && !t->trk_HItight[iTrk])
        continue;
      //else if (!doHITightVar && !t->trk_TightPrimary[iTrk])
      else if (!doHITightVar && !t->trk_HIloose[iTrk])
        continue; // track selection criteria

      if (isPbPb) {
        if (fabs (t->trk_d0sig[iTrk]) > 3.0)
          continue; // d0 significance cut in PbPb
        if (fabs (t->trk_z0sig[iTrk]) > 3.0)
          continue; //z0 significance cut in PbPb
      }
      //else {
      //  if (fabs (t->trk_d0[iTrk]) > 0.47*exp (-0.15 * t->trk_pt[iTrk]) + 0.19*exp (0.00034 * t->trk_pt[iTrk]))
      //    continue; // d0 cut in pp
      //  if (t->trk_nSCTSharedHits[iTrk] > 1)
      //    continue; // maximum no. of shared hits in pp
      //}


      if (t->trk_charge[iTrk] == 0)
        continue;
      if (fabs (t->trk_eta[iTrk]) > 2.5)
        continue; // eta cut on track

      //if (isPbPb && !is2015data) {
      //  if ((isHijing && t->trk_pt[iTrk] > 10) || (!isHijing && t->trk_pt[iTrk] < 10))
      //    continue; // 10 GeV cut: use Hijing below 10 GeV and Pythia+DO above 10 GeV
      //}

      short iEta = 0;
      while (etaTrkBins[iEta+1] < fabs (t->trk_eta[iTrk]))
        iEta++;

      //double mindr = 1.;
      //for (int iE = 0; iE < t->truth_electron_n; iE++)
      //  mindr = fmin (mindr, DeltaR (t->trk_eta[iTrk], t->truth_electron_eta[iE], t->trk_phi[iTrk], t->truth_electron_phi[iE]));
      //for (int iM = 0; iM < t->truth_muon_n; iM++)
      //  mindr = fmin (mindr, DeltaR (t->trk_eta[iTrk], t->truth_muon_eta[iM], t->trk_phi[iTrk], t->truth_muon_phi[iM]));
      //if (mindr < 0.03)
      //  continue;
      //bool isLepton = false;
      //for (int iE = 0; !isLepton && iE < t->electron_n; iE++)
      //  isLepton = isLepton || (t->trk_pt[iTrk] == t->electron_id_track_pt[iE] && t->trk_eta[iTrk] == t->electron_id_track_eta[iE] && t->trk_phi[iTrk] == t->electron_id_track_phi[iE]);
      //for (int iM = 0; !isLepton && iM < t->muon_n; iM++)
      //  isLepton = isLepton || (t->trk_pt[iTrk] == t->muon_id_track_pt[iM] && t->trk_eta[iTrk] == t->muon_id_track_eta[iM] && t->trk_phi[iTrk] == t->muon_id_track_phi[iM]);
      //if (isLepton)
      //  continue;
      //if (hasTruthZee && iL1 != -1 && iL2 != -1) {
      //  if (DeltaR (t->trk_eta[iTrk], t->truth_electron_eta[iL1], t->trk_phi[iTrk], t->truth_electron_phi[iL1]) < 0.10 || DeltaR (t->trk_eta[iTrk], t->truth_electron_eta[iL2], t->trk_phi[iTrk], t->truth_electron_phi[iL2]) < 0.10)
      //    continue;
      //}
      //else if (hasTruthZmumu && iL1 != -1 && iL2 != -1) {
      //  if (DeltaR (t->trk_eta[iTrk], t->truth_muon_eta[iL1], t->trk_phi[iTrk], t->truth_muon_phi[iL1]) < 0.10 || DeltaR (t->trk_eta[iTrk], t->truth_muon_eta[iL2], t->trk_phi[iTrk], t->truth_muon_phi[iL2]) < 0.10)
      //    continue;
      //}

      const bool isTruthMatched = (t->trk_prob_truth[iTrk] > 0.5);

      const bool isFake = !isTruthMatched;
      //if (isPbPb && !isHijing && isFake)
      //  continue;
      //
      if (!isFake && !t->trk_truth_isHadron[iTrk])
        continue;

      const bool isSecondary = isTruthMatched && (t->trk_truth_barcode[iTrk] <= 0 || 200000 <= t->trk_truth_barcode[iTrk]);

      // primary tracks are non-fake, non-secondary tracks. Strange baryons are excluded too.
      const bool isPrimary = !isFake && !isSecondary && (abs (t->trk_truth_pdgid[iTrk]) != 3112 && abs (t->trk_truth_pdgid[iTrk]) != 3222 && abs (t->trk_truth_pdgid[iTrk]) != 3312 && abs (t->trk_truth_pdgid[iTrk]) != 3334);

      if (isPrimary) {
        h2_primary_reco_tracks[iCent]->Fill (t->trk_eta[iTrk], t->trk_pt[iTrk], eventWeight);
        h_primary_reco_tracks[iCent][iEta]->Fill (t->trk_pt[iTrk], eventWeight);
      }
      h2_reco_tracks[iCent]->Fill (t->trk_eta[iTrk], t->trk_pt[iTrk], eventWeight);
      h_reco_tracks[iCent][iEta]->Fill (t->trk_pt[iTrk], eventWeight);
    }

  } // end event loop
  cout << endl << "Info: In TrackingPurity.cxx: Finished event loop." << endl;


  for (int iCent = 0; iCent < numCentBins; iCent++) {
    if ((isPbPb && iCent == 0) || (!isPbPb && iCent != 0)) continue;

    h2_primary_reco_tracks[iCent]->Write ();
    delete h2_primary_reco_tracks[iCent];
    h2_primary_reco_tracks[iCent] = nullptr;

    h2_reco_tracks[iCent]->Write ();
    delete h2_reco_tracks[iCent];
    h2_reco_tracks[iCent] = nullptr;

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h_primary_reco_tracks[iCent][iEta]->Write ();
      delete h_primary_reco_tracks[iCent][iEta];
      h_primary_reco_tracks[iCent][iEta] = nullptr;

      h_reco_tracks[iCent][iEta]->Write ();
      delete h_reco_tracks[iCent][iEta];
      h_primary_reco_tracks[iCent][iEta] = nullptr;
    }

    delete [] h_primary_reco_tracks[iCent];
    delete [] h_reco_tracks[iCent];
  }

  delete [] h_primary_reco_tracks;
  h_primary_reco_tracks = nullptr;

  delete [] h_reco_tracks;
  h_reco_tracks = nullptr;

  delete [] h2_primary_reco_tracks;
  h2_primary_reco_tracks = nullptr;

  delete [] h2_reco_tracks;
  h2_reco_tracks = nullptr;

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
