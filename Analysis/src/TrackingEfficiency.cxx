#ifndef __TrackingEfficiency_cxx__
#define __TrackingEfficiency_cxx__

#include "TrackingEfficiency.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"

#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>

#include <iostream>

using namespace std;

namespace ZHadronAnalysis {

bool TrackingEfficiency (const char* directory,
                         const int dataSet,
                         const char* inFileName,
                         const char* eventWeightsFileName) {
 
  cout << "Info: In TrackingEfficiency.cxx: Entered TrackingEfficiency routine." << endl;
  cout << "Info: In TrackingEfficiency.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar << endl;
  cout << "\n\tdoPionsOnlyVar: " << doPionsOnlyVar << endl;

  SetupDirectories ("TrackingEfficiencies");

  if (!isMC) {
    cout << "Error: In TrackingEfficiency.cxx: Trying to calculate tracking efficiency in data! Quitting." << endl;
    return false;
  }

  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TrackingEfficiency.cxx: Running over data overlay, will check data conditions" << endl;
  if (isHijing)
    cout << "Info: In TrackingEfficiency.cxx: Running over Hijing sample" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TrackingEfficiency.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TrackingEfficiency.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    cout << "Error: In TrackingEfficiency.cxx: Cannot identify this MC file! Quitting." << endl;
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
    cout << "Error: In TrackingEfficiency.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }

  TFile* eventWeightsFile = nullptr;
  TH1D* h_weights = nullptr;

  //if (isHijing) {
    eventWeightsFile = new TFile (eventWeightsFileName, "read");
    h_weights = (TH1D*) eventWeightsFile->Get (Form ("h_PbPb%s_weights_%s", doNchWeighting ? "Nch" : "FCal", isHijing ? "hijing" : "mc"));
  //}
  //else {
  //  eventWeightsFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/rootFiles/MCAnalysis/pythiaNchWeights.root", "read");
  //  h_weights = (TH1D*) eventWeightsFile->Get ("h_PythiaNchWeights");
  //}
  cout << "Info: In TrackingEfficiency.cxx: Found FCal weighting histogram, " << h_weights->GetName () << endl;

  //First sort jets & tracks into many, smaller TTrees.
  //This is where the sorting based on event information (e.g. centrality, Ntrk, jet pT) will go.
  //Event mixing will take place based on these categories so that total memory usage at any point in time is minimized.

  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetZdc ();
  t->SetGetVertices ();
  t->SetGetTracks ();
  t->SetGetTruthTracks ();
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

  cout << "Info : In TrackingEfficiency.cxx: Saving histograms to " << Form ("%s/%s/%s.root", rootPath.Data (), isHijing ? "../Hijing/" : "", identifier.Data ()) << endl;
  TFile* outFile = new TFile (Form ("%s/%s/%s.root", rootPath.Data (), isHijing ? "../Hijing/" : "", identifier.Data ()), "recreate");

  const double centBins[4] = {66.402, 1378.92, 2995.94, 5000};
  //const int centCuts[4] = {80, 30, 10, 0};
  //const double centBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
  ////const int centCuts[10] = {80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  const int numCentBins = sizeof (centBins) / sizeof (centBins[0]);

  const double ipBins[4] = {14.031, 8.582, 4.952, 0};
  const int numIPBins = sizeof (ipBins) / sizeof (ipBins[0]);

  const float multBins[6] = {-0.5, 499.5, 699.5, 1499.5, 1799.5, 3399.5};
  const int numMultBins = sizeof (multBins) / sizeof (multBins[0]) - 1;

  const int numFinerEtachBins = 40;
  const double* finerEtachBins = linspace (-2.5, 2.5, numFinerEtachBins);

  const double etachBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  const int numEtachBins = sizeof (etachBins) / sizeof (etachBins[0]) - 1;

  //const double pTchBins[13] = {0.5, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 10, 15, 30, 60};
  //const int numPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

  const double pTchBinsPP[30] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60};
  const int numPtchBinsPP = sizeof (pTchBinsPP) / sizeof (pTchBinsPP[0]) - 1;
  const double pTchBinsPbPb[28] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 15, 20, 30, 40, 60};
  const int numPtchBinsPbPb = sizeof (pTchBinsPbPb) / sizeof (pTchBinsPbPb[0]) - 1;

  //const int numPtchBins = 20;
  //const double* pTchBins = logspace (0.5, 65, numPtchBins);

  TH1D*  h_truth_matching_prob = new TH1D (Form ("h_truth_matching_prob_%s", isPbPb ? "PbPb" : "pp"), ";Truth matching prob.;N_{ch}^{rec}", 200, 0, 1);
  h_truth_matching_prob->Sumw2 ();
  TH2D** h_truth_matched_reco_tracks = new TH2D*[numCentBins];
  TH2D** h_truth_tracks = new TH2D*[numCentBins];
  TH1D** h_num = new TH1D*[numCentBins];
  TH1D** h_den = new TH1D*[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h_truth_matched_reco_tracks[iCent] = new TH2D (Form ("h_truth_matched_reco_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtachBins, finerEtachBins, (iCent == 0 ? numPtchBinsPP : numPtchBinsPbPb), (iCent == 0 ? pTchBinsPP : pTchBinsPbPb));
    h_truth_matched_reco_tracks[iCent]->Sumw2 ();

    h_truth_tracks[iCent] = new TH2D (Form ("h_truth_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtachBins, finerEtachBins, (iCent == 0 ? numPtchBinsPP : numPtchBinsPbPb), (iCent == 0 ? pTchBinsPP : pTchBinsPbPb));
    h_truth_tracks[iCent]->Sumw2 ();

  }

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TrackingEfficiency.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
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
    //if (!hasPrimary || fabs (vz) > 150)
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
    //short iMult = 0;
    //if (isPbPb) {
    //  while (iMult < numMultBins) {
    //    if (t->ntrk < multBins[iMult+1])
    //      break;
    //    else
    //      iMult++;
    //  }
    //}
    //if (iMult < 0 || iMult > numMultBins-1)
    //  continue;

    //const float eventWeight = 1;
    const float eventWeight = ((isPbPb && !isHijing) ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? t->ntrk : t->fcalA_et+t->fcalC_et)) : 1); // weight is 1 for pp
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

      if (t->trk_pt[iTrk] > 1)
        h_truth_matching_prob->Fill (t->trk_prob_truth[iTrk]);

      const bool isTruthMatched = (t->trk_prob_truth[iTrk] > 0.5);

      const bool isFake = !isTruthMatched;
      const bool isSecondary = isTruthMatched && (t->trk_truth_barcode[iTrk] <= 0 || 200000 <= t->trk_truth_barcode[iTrk]);

      // primary tracks are non-fake, non-secondary tracks. Strange baryons are excluded too.
      const bool isPrimary = !isFake && !isSecondary && (abs (t->trk_truth_pdgid[iTrk]) != 3112 && abs (t->trk_truth_pdgid[iTrk]) != 3222 && abs (t->trk_truth_pdgid[iTrk]) != 3312 && abs (t->trk_truth_pdgid[iTrk]) != 3334);

      if (!isPrimary)
        continue; // restrict to only primary tracks

      if (t->trk_truth_charge[iTrk] == 0)
        continue;
      if (fabs (t->trk_truth_eta[iTrk]) > 2.5)
        continue;
      if (!(t->trk_truth_isHadron[iTrk]))
        continue;
      if (doPionsOnlyVar && abs (t->trk_truth_pdgid[iTrk]) != 211)
        continue; //just pions
      //if (fabs (t->trk_truth_pdgid[iTrk]) != 211)
      //  continue; //just pions

      //h_truth_matched_reco_tracks[iCent]->Fill (t->trk_eta[iTrk], t->trk_truth_pt[iTrk], 1);
      h_truth_matched_reco_tracks[iCent]->Fill (t->trk_truth_eta[iTrk], t->trk_truth_pt[iTrk], eventWeight);
      //h_truth_matched_reco_tracks[iCent]->Fill (t->trk_eta[iTrk], t->trk_pt[iTrk], eventWeight); // do an unfolding
    }


    for (int iTTrk = 0; iTTrk < t->truth_trk_n; iTTrk++) {
      if (t->truth_trk_barcode[iTTrk] <= 0 || 200000 <= t->truth_trk_barcode[iTTrk] || abs (t->truth_trk_pdgid[iTTrk]) == 3112 || abs (t->truth_trk_pdgid[iTTrk]) == 3222 || abs (t->truth_trk_pdgid[iTTrk]) == 3312 || abs (t->truth_trk_pdgid[iTTrk]) == 3334)
        continue; // select primary truth particles

      if (fabs (t->truth_trk_pdgid[iTTrk]) == 11 || fabs (t->truth_trk_pdgid[iTTrk]) == 13)
        continue;
      if (fabs (t->truth_trk_eta[iTTrk]) > 2.5)
        continue;
      if (!(t->truth_trk_isHadron[iTTrk]))
        continue;
      if (doPionsOnlyVar && abs (t->truth_trk_pdgid[iTTrk]) != 211)
        continue; //just pions
      //if (fabs (t->truth_trk_pdgid[iTTrk]) != 211)
      //  continue; //just pions

      //h_truth_tracks[iCent]->Fill (t->truth_trk_eta[iTTrk], t->truth_trk_pt[iTTrk], 1);
      h_truth_tracks[iCent]->Fill (t->truth_trk_eta[iTTrk], t->truth_trk_pt[iTTrk], eventWeight);
    }

  } // end event loop
  cout << endl << "Info: In TrackingEfficiency.cxx: Finished event loop." << endl;


  for (int iCent = 0; iCent < numCentBins; iCent++) {
    TH2D h2_num = *(h_truth_matched_reco_tracks[iCent]);
    TH2D h2_den = *(h_truth_tracks[iCent]);

    h_num[iCent] = new TH1D[numEtachBins];
    h_den[iCent] = new TH1D[numEtachBins];
    for (int iEtach = 0; iEtach < numEtachBins; iEtach++) {
      h_num[iCent][iEtach] = TH1D (Form ("h_trk_eff_num_iCent%i_iEta%i", iCent, iEtach), "", (iCent == 0 ? numPtchBinsPP : numPtchBinsPbPb), (iCent == 0 ? pTchBinsPP : pTchBinsPbPb));
      h_num[iCent][iEtach].Sumw2 ();
      h_den[iCent][iEtach] = TH1D (Form ("h_trk_eff_den_iCent%i_iEta%i", iCent, iEtach), "", (iCent == 0 ? numPtchBinsPP : numPtchBinsPbPb), (iCent == 0 ? pTchBinsPP : pTchBinsPbPb));
      h_den[iCent][iEtach].Sumw2 ();
    }

    for (int iFinerEtach = 0; iFinerEtach < numFinerEtachBins; iFinerEtach++) {
      const float binCenter = 0.5 * fabs (finerEtachBins[iFinerEtach] + finerEtachBins[iFinerEtach+1]);

      int iEtach = 0;
      while (iEtach < numEtachBins) {
        if (etachBins[iEtach] < binCenter && binCenter < etachBins[iEtach+1])
          break;
        iEtach++;
      }

      for (int iPtch = 1; iPtch <= (iCent == 0 ? numPtchBinsPP : numPtchBinsPbPb); iPtch++) {
        h_num[iCent][iEtach].SetBinContent (iPtch, h_num[iCent][iEtach].GetBinContent (iPtch) + h2_num.GetBinContent (iFinerEtach+1, iPtch));
        h_num[iCent][iEtach].SetBinError (iPtch, h_num[iCent][iEtach].GetBinError (iPtch) + pow (h2_num.GetBinError (iFinerEtach+1, iPtch), 2));
        h_den[iCent][iEtach].SetBinContent (iPtch, h_den[iCent][iEtach].GetBinContent (iPtch) + h2_den.GetBinContent (iFinerEtach+1, iPtch));
        h_den[iCent][iEtach].SetBinError (iPtch, h_den[iCent][iEtach].GetBinError (iPtch) + pow (h2_den.GetBinError (iFinerEtach+1, iPtch), 2));
      }
    }

    for (int iEtach = 0; iEtach < numEtachBins; iEtach++) {
      for (int ix = 1; ix <= h_num[iCent][iEtach].GetNbinsX (); ix++) {
        if (h_num[iCent][iEtach].GetBinContent (ix) > h_den[iCent][iEtach].GetBinContent (ix)) {
          cout << "Num > den at (ix, x) =" << ix << ", " << h_num[iCent][iEtach].GetBinCenter (ix) << endl;
        }
      }

      for (int iPtch = 1; iPtch <= (iCent == 0 ? numPtchBinsPP : numPtchBinsPbPb); iPtch++) {
        h_num[iCent][iEtach].SetBinError (iPtch, sqrt (h_num[iCent][iEtach].GetBinError (iPtch)));
        h_den[iCent][iEtach].SetBinError (iPtch, sqrt (h_den[iCent][iEtach].GetBinError (iPtch)));
      }
    }
  } 



  //for (int iCent = 0; iCent < numCentBins; iCent++) {
  //  if ((isPbPb && iCent == 0) || (!isPbPb && iCent != 0))
  //    continue;

  //  h_truth_matched_reco_tracks[iCent]->Write ();
  //  h_truth_tracks[iCent]->Write ();

  //  for (int iEta = 0; iEta < numEtachBins; iEta++) {
  //    h_num[iCent][iEta].Write ();
  //    h_den[iCent][iEta].Write ();
  //  }
  //}

  h_truth_matching_prob->Write ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h_truth_matched_reco_tracks[iCent]->Write ();
    delete h_truth_matched_reco_tracks[iCent];
    h_truth_matched_reco_tracks[iCent] = nullptr;

    h_truth_tracks[iCent]->Write ();
    delete h_truth_tracks[iCent];
    h_truth_tracks[iCent] = nullptr;

    for (int iEta = 0; iEta < numEtachBins; iEta++) {
      h_num[iCent][iEta].Write ();
      h_den[iCent][iEta].Write ();
    }
    delete [] h_num[iCent];
    delete [] h_den[iCent];
  }
  delete [] h_num;
  delete [] h_den;

  delete [] h_truth_matched_reco_tracks;
  h_truth_matched_reco_tracks = nullptr;

  delete [] h_truth_tracks;
  h_truth_tracks = nullptr;

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
