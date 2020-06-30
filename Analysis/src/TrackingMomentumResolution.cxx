#ifndef __TrackingMomentumResolution_cxx__
#define __TrackingMomentumResolution_cxx__

#include "TrackingMomentumResolution.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "LocalUtilities.h"

#include <Utilities.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TF1.h>

#include <iostream>

using namespace std;

namespace ZHadronAnalysis {

bool TrackingMomentumResolution (const char* directory,
                                 const int dataSet,
                                 const char* inFileName,
                                 const char* eventWeightsFileName) {
 
  cout << "Info: In TrackingMomentumResolution.cxx: Entered TrackingMomentumResolution routine." << endl;
  cout << "Info: In TrackingMomentumResolution.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoHITightVar: " << doHITightVar << endl;
  cout << "\n\tdoPionsOnlyVar: " << doPionsOnlyVar << endl;

  SetupDirectories ("TrackingMomentumResolution");

  if (!isMC) {
    cout << "Error: In TrackingMomentumResolution.cxx: Trying to calculate tracking momentum resolution in data! Quitting." << endl;
    return false;
  }

  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TrackingMomentumResolution.cxx: Running over data overlay, will check data conditions" << endl;
  if (isHijing)
    cout << "Info: In TrackingMomentumResolution.cxx: Running over Hijing sample" << endl;

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TrackingMomentumResolution.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TrackingMomentumResolution.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    cout << "Error: In TrackingMomentumResolution.cxx: Cannot identify this MC file! Quitting." << endl;
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
    cout << "Error: In TrackingMomentumResolution.cxx: TTree not obtained for given data set. Quitting." << endl;
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
  cout << "Info: In TrackingMomentumResolution.cxx: Found FCal weighting histogram, " << h_weights->GetName () << endl;

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

  cout << "Info : In TrackingMomentumResolution.cxx: Saving histograms to " << Form ("%s/%s.root", rootPath.Data (), identifier.Data ()) << endl;
  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

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

  //const double etachBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  //const int numEtachBins = sizeof (etachBins) / sizeof (etachBins[0]) - 1;

  const double pTchBinsPP[33] = {0.80, 0.84, 0.88, 0.92, 0.96, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 60, 100};
  const int numPtchBinsPP = sizeof (pTchBinsPP) / sizeof (pTchBinsPP[0]) - 1;
  const double pTchBinsPbPb[31] = {0.80, 0.84, 0.88, 0.92, 0.96, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 15, 20, 30, 60, 100};
  const int numPtchBinsPbPb = sizeof (pTchBinsPbPb) / sizeof (pTchBinsPbPb[0]) - 1;

  const double* pTchBins = (isPbPb ? pTchBinsPP : pTchBinsPbPb);
  const int numPtchBins = (isPbPb ? numPtchBinsPP : numPtchBinsPbPb);

  TH1D**** h_tms = new TH1D***[numCentBins];
  //TH2D** h2_avg_tms = new TH2D*[numCentBins];
  //TH2D** h2_avg_tmr = new TH2D*[numCentBins];
  //TH1D*** h_avg_tms = new TH1D**[numCentBins];
  //TH1D*** h_avg_tmr = new TH1D**[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    //h2_avg_tms[iCent] = new TH2D (Form ("h2_avg_tms_iCent%i", iCent), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMS [%]", numPtchBins, pTchBins, numFinerEtachBins, finerEtachBins);
    //h2_avg_tmr[iCent] = new TH2D (Form ("h2_avg_tmr_iCent%i", iCent), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMR [%]", numPtchBins, pTchBins, numFinerEtachBins, finerEtachBins);
    //h_avg_tms[iCent] = new TH1D*[numEtachBins];
    //h_avg_tmr[iCent] = new TH1D*[numEtachBins];

    //for (int iEta = 0; iEta < numEtachBins; iEta++) {
    //  h_avg_tms[iCent][iEta] = new TH1D (Form ("h_avg_tms_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T}^{truth} [GeV];TMS [%]", numPtchBins, pTchBins);
    //  h_avg_tmr[iCent][iEta] = new TH1D (Form ("h_avg_tmr_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T}^{truth} [GeV];TMR [%]", numPtchBins, pTchBins);
    //}
    

    h_tms[iCent] = new TH1D**[numPtchBins];
    for (int iPtch = 0; iPtch < numPtchBins; iPtch++) {
      h_tms[iCent][iPtch] = new TH1D*[numFinerEtachBins];
      for (int iEta = 0; iEta < numFinerEtachBins; iEta++) {
        h_tms[iCent][iPtch][iEta] = new TH1D (Form ("h_tms_iCent%i_iPtch%i_iEta%i", iCent, iPtch, iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 200, 0.5, 1.5);
        h_tms[iCent][iPtch][iEta]->Sumw2 ();
      }
    }

  }

  const int nEvts = tree->GetEntries ();

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TrackingMomentumResolution.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
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

    const float eventWeight = 1;
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

      short iPtch = -1;
      if (pTchBins[0] <= t->trk_truth_pt[iTrk]) {
        iPtch = 0;
        while (iPtch < numPtchBins && pTchBins[iPtch+1] < t->trk_truth_pt[iTrk]) iPtch++;
      }

      short iEta = -1;
      if (finerEtachBins[0] <= t->trk_truth_eta[iTrk]) {
        iEta = 0;
        while (iEta < numFinerEtachBins && finerEtachBins[iEta+1] < t->trk_truth_eta[iTrk]) iEta++;
      }

      if (iPtch >= 0 && iPtch < numPtchBins && iEta >= 0 && iEta < numFinerEtachBins) {
        h_tms[iCent][iPtch][iEta]->Fill (t->trk_pt[iTrk] / t->trk_truth_pt[iTrk]);
      }
    }
  } // end event loop
  cout << endl << "Info: In TrackingMomentumResolution.cxx: Finished event loop." << endl;



  //for (int iCent = 0; iCent < numCentBins; iCent++) {

  //  for (int iPtch = 0; iPtch < numPtchBins; iPtch++) {

  //    TH1D** h_tms_integratedEta = new TH1D*[numEtachBins];
  //    for (int iEta = 0; iEta < numEtachBins; iEta++) {
  //      h_tms_integratedEta[iEta] = new TH1D (Form ("h_tms_integratedEta_iEta%i", iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 200, 0.5, 1.5);
  //      h_tms_integratedEta[iEta]->Sumw2 ();
  //    }

  //    for (int iEta = 0; iEta < numFinerEtachBins; iEta++) {

  //      if (h_tms[iCent][iPtch][iEta]->Integral () == 0) continue;

  //      // first add to eta-integrated plot
  //      const float binCenter = 0.5 * fabs (finerEtachBins[iEta] + finerEtachBins[iEta+1]);
  //      int iEtach = 0;
  //      while (iEtach < numEtachBins) {
  //        if (etachBins[iEtach] < binCenter && binCenter < etachBins[iEtach+1])
  //          break;
  //        iEtach++;
  //      }

  //      h_tms_integratedEta[iEta]->Add (h_tms[iCent][iPtch][iEta]);

  //      TF1* fit = new TF1 ("fit", "gaus(0)", 0.5, 1.5);
  //      fit->SetParameter (0, h_tms[iCent][iPtch][iEta]->Integral ());
  //      fit->SetParameter (1, 1);
  //      fit->SetParameter (2, 1);

  //      h_tms[iCent][iPtch][iEta]->Fit (fit, "RN0Q");

  //      double mean = fit->GetParameter (1);
  //      double sigma = fit->GetParameter (2);

  //      delete fit;
  //      fit = new TF1 ("fit", "gaus(0)", mean-2*sigma, mean+2*sigma);

  //      fit->SetParameter (0, h_tms[iCent][iPtch][iEta]->Integral ());
  //      fit->SetParameter (1, mean);
  //      fit->SetParameter (2, sigma);

  //      h_tms[iCent][iPtch][iEta]->Fit (fit, "RN0Q");

  //      mean = fit->GetParameter (1);
  //      sigma = fit->GetParameter (2);

  //      delete fit;

  //      const double tms = mean;
  //      const double tmr = sigma / mean;

  //      h2_avg_tms[iCent]->SetBinContent (iPtch+1, iEta+1, tms);
  //      h2_avg_tmr[iCent]->SetBinContent (iPtch+1, iEta+1, tmr);
  //    }

  //    for (int iEta = 0; iEta < numEtachBins; iEta++) {
  //      TF1* fit = new TF1 ("fit", "gaus(0)", 0.5, 1.5);
  //      fit->SetParameter (0, h_tms_integratedEta[iEta]->Integral ());
  //      fit->SetParameter (1, 1);
  //      fit->SetParameter (2, 1);

  //      h_tms_integratedEta[iEta]->Fit (fit, "RN0Q");

  //      double mean = fit->GetParameter (1);
  //      double sigma = fit->GetParameter (2);

  //      delete fit;
  //      fit = new TF1 ("fit", "gaus(0)", mean-2*sigma, mean+2*sigma);

  //      fit->SetParameter (0, h_tms_integratedEta[iEta]->Integral ());
  //      fit->SetParameter (1, mean);
  //      fit->SetParameter (2, sigma);

  //      h_tms_integratedEta[iEta]->Fit (fit, "RN0Q");

  //      mean = fit->GetParameter (1);
  //      sigma = fit->GetParameter (2);

  //      delete fit;

  //      const double tms = mean;
  //      const double tmr = sigma / mean;

  //      h_avg_tms[iCent][iEta]->SetBinContent (iPtch+1, tms);
  //      h_avg_tmr[iCent][iEta]->SetBinContent (iPtch+1, tmr);

  //      delete h_tms_integratedEta[iEta];
  //    } // end loop over iEta
  //
  //    delete[] h_tms_integratedEta;
  //  } // end loop over iPtch
  //} // end loop over iCent



  for (int iCent = 0; iCent < numCentBins; iCent++) {
    for (int iPtch = 0; iPtch < numPtchBins; iPtch++) {
      for (int iEta = 0; iEta < numFinerEtachBins; iEta++) {
        h_tms[iCent][iPtch][iEta]->Write ();
      } // end loop over iEta
    } // end loop over iPtch
  } // end loop over iCent


  //  h2_avg_tms[iCent]->Write ();
  //  h2_avg_tmr[iCent]->Write ();
  //  for (int iEta = 0; iEta < numEtachBins; iEta++) {
  //    h_avg_tms[iCent][iEta]->Write ();
  //    h_avg_tmr[iCent][iEta]->Write ();
  //  } // end loop over iEta
  //} // end loop over iCent



  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  //SaferDelete (&outFile);

  return true;
}

} // end namespace

#endif
