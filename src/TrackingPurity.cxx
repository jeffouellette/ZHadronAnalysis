#include "TrackingPurity.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
//#include <TEfficiency.h>
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
    cout << "Info: In TreeMaker.cxx: Running over data overlay, will check data conditions" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TrackingPurity.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TrackingPurity.cxx: Saving output to " << rootPath << endl;

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

  //First sort jets & tracks into many, smaller TTrees.
  //This is where the sorting based on event information (e.g. centrality, Ntrk, jet pT) will go.
  //Event mixing will take place based on these categories so that total memory usage at any point in time is minimized.

  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetTracks ();
  t->SetGetTruthTracks ();
  t->SetBranchAddresses ();

  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  const double centBins[4] = {66.402, 1378.92, 2995.94, 5000};
  ////const int centCuts[4] = {80, 30, 10, 0};
  //const double centBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
  //const int centCuts[10] = {80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  const int numCentBins = sizeof (centBins) / sizeof (centBins[0]);

  const int numFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, numFinerEtaTrkBins);

  const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  const int numEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

  const double ptTrkBins[13] = {0.5, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 10, 15, 30, 60};
  const int numPtTrkBins = sizeof (ptTrkBins) / sizeof (ptTrkBins[0]) - 1;

  //const int numPtTrkBins = 20;
  //const double* ptTrkBins = logspace (0.5, 60, numPtTrkBins);

  //TEfficiency*** h_trk_pur = new TEfficiency**[numCentBins];
  TH2D** h2_primary_reco_tracks = new TH2D*[numCentBins];
  TH2D** h2_reco_tracks = new TH2D*[numCentBins];
  TH1D*** h_primary_reco_tracks = new TH1D**[numCentBins];
  TH1D*** h_reco_tracks = new TH1D**[numCentBins];

  TLorentzVector l1, l2;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h2_primary_reco_tracks[iCent] = new TH2D (Form ("h2_primary_reco_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtTrkBins, ptTrkBins);
    h2_primary_reco_tracks[iCent]->Sumw2 ();

    h2_reco_tracks[iCent] = new TH2D (Form ("h2_reco_tracks_iCent%i", iCent), ";#eta;#it{p}_{T} [GeV]", numFinerEtaTrkBins, finerEtaTrkBins, numPtTrkBins, ptTrkBins);
    h2_reco_tracks[iCent]->Sumw2 ();

    h_primary_reco_tracks[iCent] = new TH1D*[numEtaTrkBins];
    h_reco_tracks[iCent] = new TH1D*[numEtaTrkBins];

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h_primary_reco_tracks[iCent][iEta] = new TH1D (Form ("h_primary_reco_tracks_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T} [GeV]", numPtTrkBins, ptTrkBins);
      h_primary_reco_tracks[iCent][iEta]->Sumw2 ();

      h_reco_tracks[iCent][iEta] = new TH1D (Form ("h_reco_tracks_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T} [GeV]", numPtTrkBins, ptTrkBins);
      h_reco_tracks[iCent][iEta]->Sumw2 ();
    }
  }

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TrackingPurity.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
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

    const float eventWeight = (isPbPb ? h_weights->GetBinContent (h_weights->FindFixBin (doNchWeighting ? t->ntrk : t->fcalA_et+t->fcalC_et)) : 1); // weight is 1 for pp

    for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
      if (doHITightVar && !t->trk_HItight[iTrk])
        continue;
      else if (!doHITightVar && !t->trk_HIloose[iTrk])
        continue; // track selection criteria

      if (0 == t->trk_truth_charge[iTrk])
        continue;
      if (fabs (t->trk_truth_eta[iTrk]) > 2.5)
        continue; // eta cut on truth track
      if (!(t->trk_truth_isHadron[iTrk]))
        continue; // require hadrons

      if (isPbPb && !is2015data) {
        if ((isHijing && t->trk_pt[iTrk] > 10) || (!isHijing && t->trk_pt[iTrk] < 10))
          continue; // 10 GeV cut: use Hijing below 10 GeV and Pythia+DO above 10 GeV
      }

      short iEta = 0;
      while (etaTrkBins[iEta+1] < fabs (t->trk_eta[iTrk]))
        iEta++;

      const bool isFake = (t->trk_prob_truth[iTrk] < 0.5);
      const bool isPrimary = !isFake && (0 <= t->trk_truth_barcode[iTrk]) && (t->trk_truth_barcode[iTrk] < 200000);
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
    if ((isPbPb && iCent == 0) || (!isPbPb && iCent != 0))
      continue;

    h2_primary_reco_tracks[iCent]->Write ();
    h2_reco_tracks[iCent]->Write ();

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h_primary_reco_tracks[iCent][iEta]->Write ();
      h_reco_tracks[iCent][iEta]->Write ();
    }
  }

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (outFile);


  for (int iCent = 0; iCent < numCentBins; iCent++) {
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

  return true;
}

} // end namespace
