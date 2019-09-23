#include "TrackingEfficiency.h"
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
    cout << "Info: In TreeMaker.cxx: Running over data overlay, will check data conditions" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TrackingEfficiency.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TrackingEfficiency.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    cout << "Error: In TreeMaker.C: Cannot identify this MC file! Quitting." << endl;
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
    cout << "Error: In TagAndProbe.cxx: TTree not obtained for given data set. Quitting." << endl;
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

  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  const double centBins[4] = {66.402, 1378.92, 2995.94, 5000};
  //const int centCuts[4] = {80, 30, 10, 0};
  //const double centBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
  ////const int centCuts[10] = {80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  const int numCentBins = sizeof (centBins) / sizeof (centBins[0]);

  const float multBins[6] = {-0.5, 499.5, 699.5, 1499.5, 1799.5, 3399.5};
  const int numMultBins = sizeof (multBins) / sizeof (multBins[0]) - 1;

  const int numFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, numFinerEtaTrkBins);

  const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  const int numEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

  const double ptTrkBins[13] = {0.5, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 10, 15, 30, 60};
  const int numPtTrkBins = sizeof (ptTrkBins) / sizeof (ptTrkBins[0]) - 1;

  //const int numPtTrkBins = 20;
  //const double* ptTrkBins = logspace (0.5, 65, numPtTrkBins);

  TH2D** h_truth_matched_reco_tracks = new TH2D*[numCentBins];
  TH2D** h_truth_tracks = new TH2D*[numCentBins];
  TH1D** h_num = new TH1D*[numCentBins];
  TH1D** h_den = new TH1D*[numCentBins];

  TLorentzVector l1, l2;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h_truth_matched_reco_tracks[iCent] = new TH2D (Form ("h_truth_matched_reco_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtTrkBins, ptTrkBins);
    h_truth_matched_reco_tracks[iCent]->Sumw2 ();

    h_truth_tracks[iCent] = new TH2D (Form ("h_truth_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtTrkBins, ptTrkBins);
    h_truth_tracks[iCent]->Sumw2 ();

  }

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TrackingEfficiency.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    if (isOverlayMC) {
      if (collisionSystem == PbPb18 && (t->BlayerDesyn || !t->passes_toroid))
        continue;
      if (isPbPb && t->isOOTPU)
        continue;
    }

    bool hasPrimary = false;
    for (int iVert = 0; !hasPrimary && iVert < t->nvert; iVert++) {
      hasPrimary = (t->vert_type[iVert] == 1);
    }
    if (!hasPrimary)
      continue;

    short iCent = 0;
    if (isPbPb) {
      const float fcal_et = t->fcalA_et + t->fcalC_et;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
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
    //const float eventWeight = ((isPbPb && !isHijing) ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? t->ntrk : t->fcalA_et+t->fcalC_et)) : 1); // weight is 1 for pp
    const float eventWeight = (isPbPb ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? t->ntrk : t->fcalA_et+t->fcalC_et)) : 1); // weight is 1 for pp

    for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      if (doHITightVar && !t->trk_HItight[iTrk])
        continue;
      else if (!doHITightVar && !t->trk_HIloose[iTrk])
        continue; // track selection criteria

      if (t->trk_prob_truth[iTrk] < 0.5)
        continue; // select non-fake tracks

      if (t->trk_truth_barcode[iTrk] <= 0 || 200000 < t->trk_truth_barcode[iTrk])
        continue; // select primary tracks

      if (0 == t->trk_truth_charge[iTrk])
        continue;
      if (fabs (t->trk_truth_eta[iTrk]) > 2.5)
        continue;
      if (!(t->trk_truth_isHadron[iTrk]))
        continue;
      if (doPionsOnlyVar && fabs (t->trk_truth_pdgid[iTrk]) != 211)
        continue; //just pions

      //h_truth_matched_reco_tracks[iCent]->Fill (t->trk_eta[iTrk], t->trk_truth_pt[iTrk], 1);
      h_truth_matched_reco_tracks[iCent]->Fill (t->trk_truth_eta[iTrk], t->trk_truth_pt[iTrk], eventWeight);
      //h_truth_matched_reco_tracks[iCent]->Fill (t->trk_eta[iTrk], t->trk_pt[iTrk], eventWeight); // do an unfolding
    }


    for (int iTTrk = 0; iTTrk < t->truth_trk_n; iTTrk++) {
      if (t->truth_trk_barcode[iTTrk] <= 0 || 200000 < t->truth_trk_barcode[iTTrk])
        continue; // select primary truth particles

      if (fabs (t->truth_trk_pdgid[iTTrk]) == 11 || fabs (t->truth_trk_pdgid[iTTrk]) == 13)
        continue;
      if (fabs (t->truth_trk_eta[iTTrk]) > 2.5)
        continue;
      if (!(t->truth_trk_isHadron[iTTrk]))
        continue;
      if (doPionsOnlyVar && fabs (t->truth_trk_pdgid[iTTrk]) != 211)
        continue; //just pions

      //h_truth_tracks[iCent]->Fill (t->truth_trk_eta[iTTrk], t->truth_trk_pt[iTTrk], 1);
      h_truth_tracks[iCent]->Fill (t->truth_trk_eta[iTTrk], t->truth_trk_pt[iTTrk], eventWeight);
    }

  } // end event loop
  cout << endl << "Info: In TrackingEfficiency.cxx: Finished event loop." << endl;


  for (int iCent = 0; iCent < numCentBins; iCent++) {
    TH2D h2_num = *(h_truth_matched_reco_tracks[iCent]);
    TH2D h2_den = *(h_truth_tracks[iCent]);

    h_num[iCent] = new TH1D[numEtaTrkBins];
    h_den[iCent] = new TH1D[numEtaTrkBins];
    for (int iEtaTrk = 0; iEtaTrk < numEtaTrkBins; iEtaTrk++) {
      h_num[iCent][iEtaTrk] = TH1D (Form ("h_trk_eff_num_iCent%i_iEta%i", iCent, iEtaTrk), "", numPtTrkBins, ptTrkBins);
      h_num[iCent][iEtaTrk].Sumw2 ();
      h_den[iCent][iEtaTrk] = TH1D (Form ("h_trk_eff_den_iCent%i_iEta%i", iCent, iEtaTrk), "", numPtTrkBins, ptTrkBins);
      h_den[iCent][iEtaTrk].Sumw2 ();
    }

    for (int iFinerEtaTrk = 0; iFinerEtaTrk < numFinerEtaTrkBins; iFinerEtaTrk++) {
      const float binCenter = 0.5 * fabs (finerEtaTrkBins[iFinerEtaTrk] + finerEtaTrkBins[iFinerEtaTrk+1]);

      int iEtaTrk = 0;
      while (iEtaTrk < numEtaTrkBins) {
        if (etaTrkBins[iEtaTrk] < binCenter && binCenter < etaTrkBins[iEtaTrk+1])
          break;
        iEtaTrk++;
      }

      for (int iPtTrk = 1; iPtTrk <= numPtTrkBins; iPtTrk++) {
        h_num[iCent][iEtaTrk].SetBinContent (iPtTrk, h_num[iCent][iEtaTrk].GetBinContent (iPtTrk) + h2_num.GetBinContent (iFinerEtaTrk+1, iPtTrk));
        h_num[iCent][iEtaTrk].SetBinError (iPtTrk, h_num[iCent][iEtaTrk].GetBinError (iPtTrk) + pow (h2_num.GetBinError (iFinerEtaTrk+1, iPtTrk), 2));
        h_den[iCent][iEtaTrk].SetBinContent (iPtTrk, h_den[iCent][iEtaTrk].GetBinContent (iPtTrk) + h2_den.GetBinContent (iFinerEtaTrk+1, iPtTrk));
        h_den[iCent][iEtaTrk].SetBinError (iPtTrk, h_den[iCent][iEtaTrk].GetBinError (iPtTrk) + pow (h2_den.GetBinError (iFinerEtaTrk+1, iPtTrk), 2));
      }
    }

    for (int iEtaTrk = 0; iEtaTrk < numEtaTrkBins; iEtaTrk++) {
      for (int ix = 1; ix <= h_num[iCent][iEtaTrk].GetNbinsX (); ix++) {
        if (h_num[iCent][iEtaTrk].GetBinContent (ix) > h_den[iCent][iEtaTrk].GetBinContent (ix)) {
          cout << "Num > den at (ix, x) =" << ix << ", " << h_num[iCent][iEtaTrk].GetBinCenter (ix) << endl;
        }
      }

      for (int iPtTrk = 1; iPtTrk <= numPtTrkBins; iPtTrk++) {
        h_num[iCent][iEtaTrk].SetBinError (iPtTrk, sqrt (h_num[iCent][iEtaTrk].GetBinError (iPtTrk)));
        h_den[iCent][iEtaTrk].SetBinError (iPtTrk, sqrt (h_den[iCent][iEtaTrk].GetBinError (iPtTrk)));
      }
    }
  } 



  for (int iCent = 0; iCent < numCentBins; iCent++) {
    if ((isPbPb && iCent == 0) || (!isPbPb && iCent != 0))
      continue;

    h_truth_matched_reco_tracks[iCent]->Write ();
    h_truth_tracks[iCent]->Write ();

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h_num[iCent][iEta].Write ();
      h_den[iCent][iEta].Write ();
    }
  }

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (outFile);


  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h_truth_matched_reco_tracks[iCent]->Write ();
    delete h_truth_matched_reco_tracks[iCent];
    h_truth_matched_reco_tracks[iCent] = nullptr;

    h_truth_tracks[iCent]->Write ();
    delete h_truth_tracks[iCent];
    h_truth_tracks[iCent] = nullptr;

    delete [] h_num[iCent];
    delete [] h_den[iCent];
  }
  delete [] h_num;
  delete [] h_den;

  delete [] h_truth_matched_reco_tracks;
  h_truth_matched_reco_tracks = nullptr;

  delete [] h_truth_tracks;
  h_truth_tracks = nullptr;

  return true;
}

} // end namespace
