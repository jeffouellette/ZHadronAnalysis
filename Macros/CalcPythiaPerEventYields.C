#ifndef __CalcPythiaPerEventYields_C__
#define __CalcPythiaPerEventYields_C__

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TF1.h>
#include <TRandom3.h>

#include <Utilities.h>

#include "Params.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <math.h>

void DoStuff (TString inFileName, TString outFileName) {

  //int code = 0;
  //int id1 = 0;
  //int id2 = 0;
  //float x1pdf = 0;
  //float x2pdf = 0;
  //float Q = 0;
  //bool isValence1 = false;
  //bool isValence2 = false;

  float z_pt = 0;
  //float z_eta = 0;
  float z_phi = 0;
  //float z_m = 0;
  float z_pt_mixed = 0;
  //float z_eta_mixed = 0;
  float z_phi_mixed = 0;
  //float z_m_mixed = 0;
  int part_n = 0;
  float part_pt[10000];
  float part_eta[10000];
  float part_phi[10000];
  float phi_transmin = 0;
  float phi_transmax = 0;

  TH1D* h_trk_pTch_yield[3];
  TH2D* h2_trk_pTch_cov[3];
  TH1D* h_trk_FinePtch_yield[3];
  TH2D* h2_trk_FinePtch_cov[3];
  TH1D* h_trk_xhZ_yield[3];
  TH2D* h2_trk_xhZ_cov[3];
  TH1D* h_trk_FineXhZ_yield[3];
  TH2D* h2_trk_FineXhZ_cov[3];
  TH1D* h_trk_pTch_yield_bkg[3];
  TH2D* h2_trk_pTch_cov_bkg[3];
  TH1D* h_trk_FinePtch_yield_bkg[3];
  TH2D* h2_trk_FinePtch_cov_bkg[3];
  TH1D* h_trk_xhZ_yield_bkg[3];
  TH2D* h2_trk_xhZ_cov_bkg[3];
  TH1D* h_trk_FineXhZ_yield_bkg[3];
  TH2D* h2_trk_FineXhZ_cov_bkg[3];

  const int nFinePtchBins = 40;
  double* FinePtchBins[3] = {};
  FinePtchBins[0] = logspace (1., 15., nFinePtchBins);
  FinePtchBins[1] = logspace (1., 30., nFinePtchBins);
  FinePtchBins[2] = logspace (1., 60., nFinePtchBins);
  const int nFineXhZBins = 40;
  double* FineXhZBins[3] = {};
  FineXhZBins[0] = logspace (1./15., 1., nFineXhZBins);
  FineXhZBins[1] = logspace (1./30., 1., nFineXhZBins);
  FineXhZBins[2] = logspace (1./60., 1., nFineXhZBins);

  float trk_counts_pTch[maxNPtchBins] = {};
  float trk_counts_FinePtch[nFinePtchBins] = {};
  float trk_counts_xhZ[maxNXhZBins] = {};
  float trk_counts_FineXhZ[nFineXhZBins] = {};

  TFile* outFile = new TFile (Form ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TruthAnalysis/Nominal/%s", outFileName.Data ()), "recreate");

  TH1D* h_z_counts = new TH1D ("h_z_counts", "iPtZ", nPtZBins, 0, nPtZBins);

  for (int iPtZ = 0; iPtZ < 3; iPtZ++) {
    h_trk_pTch_yield[iPtZ] = new TH1D (Form ("h_trk_pTch_yield_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV]", nPtchBins[iPtZ+2], pTchBins[iPtZ+2]);
    h_trk_pTch_yield[iPtZ]->Sumw2 ();
    h_trk_xhZ_yield[iPtZ] = new TH1D (Form ("h_trk_xhZ_yield_iPtZ%i", iPtZ+2), ";#it{x}_{hZ}", nXhZBins[iPtZ+2], xhZBins[iPtZ+2]);
    h_trk_xhZ_yield[iPtZ]->Sumw2 ();
    h2_trk_pTch_cov[iPtZ] = new TH2D (Form ("h2_trk_pTch_cov_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV]", nPtchBins[iPtZ+2], pTchBins[iPtZ+2], nPtchBins[iPtZ+2], pTchBins[iPtZ+2]);
    h2_trk_pTch_cov[iPtZ]->Sumw2 ();
    h2_trk_xhZ_cov[iPtZ] = new TH2D (Form ("h2_trk_xhZ_cov_iPtZ%i", iPtZ+2), ";#it{x}_{hZ};#it{x}_{hZ}", nXhZBins[iPtZ+2], xhZBins[iPtZ+2], nXhZBins[iPtZ+2], xhZBins[iPtZ+2]);
    h2_trk_xhZ_cov[iPtZ]->Sumw2 ();

    h_trk_FinePtch_yield[iPtZ] = new TH1D (Form ("h_trk_FinePtch_yield_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV]", nFinePtchBins, FinePtchBins[iPtZ]);
    h_trk_FinePtch_yield[iPtZ]->Sumw2 ();
    h_trk_FineXhZ_yield[iPtZ] = new TH1D (Form ("h_trk_FineXhZ_yield_iPtZ%i", iPtZ+2), ";#it{x}_{hZ}", nFineXhZBins, FineXhZBins[iPtZ]);
    h_trk_FineXhZ_yield[iPtZ]->Sumw2 ();
    h2_trk_FinePtch_cov[iPtZ] = new TH2D (Form ("h2_trk_FinePtch_cov_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV]", nFinePtchBins, FinePtchBins[iPtZ], nFinePtchBins, FinePtchBins[iPtZ]);
    h2_trk_FinePtch_cov[iPtZ]->Sumw2 ();
    h2_trk_FineXhZ_cov[iPtZ] = new TH2D (Form ("h2_trk_FineXhZ_cov_iPtZ%i", iPtZ+2), ";#it{x}_{hZ};#it{x}_{hZ}", nFineXhZBins, FineXhZBins[iPtZ], nFineXhZBins, FineXhZBins[iPtZ]);
    h2_trk_FineXhZ_cov[iPtZ]->Sumw2 ();

    h_trk_pTch_yield_bkg[iPtZ] = new TH1D (Form ("h_trk_pTch_yield_bkg_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV]", nPtchBins[iPtZ+2], pTchBins[iPtZ+2]);
    h_trk_pTch_yield_bkg[iPtZ]->Sumw2 ();
    h_trk_xhZ_yield_bkg[iPtZ] = new TH1D (Form ("h_trk_xhZ_yield_bkg_iPtZ%i", iPtZ+2), ";#it{x}_{hZ}", nXhZBins[iPtZ+2], xhZBins[iPtZ+2]);
    h_trk_xhZ_yield_bkg[iPtZ]->Sumw2 ();
    h2_trk_pTch_cov_bkg[iPtZ] = new TH2D (Form ("h2_trk_pTch_cov_bkg_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV]", nPtchBins[iPtZ+2], pTchBins[iPtZ+2], nPtchBins[iPtZ+2], pTchBins[iPtZ+2]);
    h2_trk_pTch_cov_bkg[iPtZ]->Sumw2 ();
    h2_trk_xhZ_cov_bkg[iPtZ] = new TH2D (Form ("h2_trk_xhZ_cov_bkg_iPtZ%i", iPtZ+2), ";#it{x}_{hZ};#it{x}_{hZ}", nXhZBins[iPtZ+2], xhZBins[iPtZ+2], nXhZBins[iPtZ+2], xhZBins[iPtZ+2]);
    h2_trk_xhZ_cov_bkg[iPtZ]->Sumw2 ();

    h_trk_FinePtch_yield_bkg[iPtZ] = new TH1D (Form ("h_trk_FinePtch_yield_bkg_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV]", nFinePtchBins, FinePtchBins[iPtZ]);
    h_trk_FinePtch_yield_bkg[iPtZ]->Sumw2 ();
    h_trk_FineXhZ_yield_bkg[iPtZ] = new TH1D (Form ("h_trk_FineXhZ_yield_bkg_iPtZ%i", iPtZ+2), ";#it{x}_{hZ}", nFineXhZBins, FineXhZBins[iPtZ]);
    h_trk_FineXhZ_yield_bkg[iPtZ]->Sumw2 ();
    h2_trk_FinePtch_cov_bkg[iPtZ] = new TH2D (Form ("h2_trk_FinePtch_cov_bkg_iPtZ%i", iPtZ+2), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV]", nFinePtchBins, FinePtchBins[iPtZ], nFinePtchBins, FinePtchBins[iPtZ]);
    h2_trk_FinePtch_cov_bkg[iPtZ]->Sumw2 ();
    h2_trk_FineXhZ_cov_bkg[iPtZ] = new TH2D (Form ("h2_trk_FineXhZ_cov_bkg_iPtZ%i", iPtZ+2), ";#it{x}_{hZ};#it{x}_{hZ}", nFineXhZBins, FineXhZBins[iPtZ], nFineXhZBins, FineXhZBins[iPtZ]);
    h2_trk_FineXhZ_cov_bkg[iPtZ]->Sumw2 ();
  }


  TFile* inFile = new TFile (Form ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TruthAnalysis/Nominal/%s", inFileName.Data ()), "read");
  TTree* inTree = (TTree*) inFile->Get ("ppZTrackTree");

  inTree->LoadBaskets (2000000000);

  //inTree->SetBranchAddress ("code",       &code);
  //inTree->SetBranchAddress ("id1",        &id1);
  //inTree->SetBranchAddress ("id2",        &id2);
  //inTree->SetBranchAddress ("x1pdf",      &x1pdf);
  //inTree->SetBranchAddress ("x2pdf",      &x2pdf);
  //inTree->SetBranchAddress ("Q",          &Q);
  //inTree->SetBranchAddress ("isValence1", &isValence1);
  //inTree->SetBranchAddress ("isValence2", &isValence2);
  inTree->SetBranchAddress ("z_pt",         &z_pt);
  //inTree->SetBranchAddress ("z_eta",        &z_eta);
  inTree->SetBranchAddress ("z_phi",        &z_phi);
  //inTree->SetBranchAddress ("z_m",          &z_m);
  inTree->SetBranchAddress ("ntrk",         &part_n);
  inTree->SetBranchAddress ("trk_pt",       &part_pt);
  inTree->SetBranchAddress ("trk_eta",      &part_eta);
  inTree->SetBranchAddress ("trk_phi",      &part_phi);
  inTree->SetBranchAddress ("phi_transmin", &phi_transmin);
  inTree->SetBranchAddress ("phi_transmax", &phi_transmax);

  const int nEvents = inTree->GetEntries ();

  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    inTree->GetEntry (iEvent);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

    if (iPtZ < 2) continue; // skip unneeded events

    h_z_counts->SetBinContent (iPtZ+1, h_z_counts->GetBinContent (iPtZ+1) + 1);

    // first loop over the particles in the recorded event
    for (int iPart = 0; iPart < part_n; iPart++) {
      if (fabs (part_eta[iPart]) > 2.5)
        continue;
      if (DeltaPhi (part_phi[iPart], z_phi) < 3*TMath::Pi ()/4)
        continue;

      const float pptch = part_pt[iPart];
      const float pxhz = pptch / z_pt;

      if (pptch >= pTchBins[iPtZ][0]) {
        int iPtch = 0;
        while (iPtch < nPtchBins[iPtZ] && pptch >= pTchBins[iPtZ][iPtch+1]) iPtch++;
        if (0 <= iPtch && iPtch < nPtchBins[iPtZ])
          trk_counts_pTch[iPtch]++;
      }
      if (pptch >= FinePtchBins[iPtZ-2][0]) {
        int iPtch = 0;
        while (iPtch < nFinePtchBins && pptch >= FinePtchBins[iPtZ-2][iPtch+1]) iPtch++;
        if (0 <= iPtch && iPtch < nFinePtchBins)
          trk_counts_FinePtch[iPtch]++;
      }
      if (pxhz >= xhZBins[iPtZ][0]) {
        int iXhZ = 0;
        while (iXhZ < nXhZBins[iPtZ] && pxhz >= xhZBins[iPtZ][iXhZ+1]) iXhZ++;
        if (0 <= iXhZ && iXhZ < nXhZBins[iPtZ])
          trk_counts_xhZ[iXhZ]++;
      }
      if (pxhz >= FineXhZBins[iPtZ-2][0]) {
        int iXhZ = 0;
        while (iXhZ < nFineXhZBins && pxhz >= FineXhZBins[iPtZ-2][iXhZ+1]) iXhZ++;
        if (0 <= iXhZ && iXhZ < nFineXhZBins)
          trk_counts_FineXhZ[iXhZ]++;
      }
    } // end loop over iPart

    TH1D* h = h_trk_pTch_yield[iPtZ-2];
    TH2D* h2 = h2_trk_pTch_cov[iPtZ-2];
    for (short iX = 0; iX < nPtchBins[iPtZ]; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_pTch[iX]);
      for (short iY = 0; iY < nPtchBins[iPtZ]; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_pTch[iX])*(trk_counts_pTch[iY]));
    }
    h = h_trk_FinePtch_yield[iPtZ-2];
    h2 = h2_trk_FinePtch_cov[iPtZ-2];
    for (short iX = 0; iX < nFinePtchBins; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_FinePtch[iX]);
      for (short iY = 0; iY < nFinePtchBins; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_FinePtch[iX])*(trk_counts_FinePtch[iY]));
    }
    h = h_trk_xhZ_yield[iPtZ-2];
    h2 = h2_trk_xhZ_cov[iPtZ-2];
    for (short iX = 0; iX < nXhZBins[iPtZ]; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_xhZ[iX]);
      for (short iY = 0; iY < nXhZBins[iPtZ]; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_xhZ[iX])*(trk_counts_xhZ[iY]));
    }
    h = h_trk_FineXhZ_yield[iPtZ-2];
    h2 = h2_trk_FineXhZ_cov[iPtZ-2];
    for (short iX = 0; iX < nFineXhZBins; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_FineXhZ[iX]);
      for (short iY = 0; iY < nFineXhZBins; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_FineXhZ[iX])*(trk_counts_FineXhZ[iY]));
    }

    for (int i = 0; i < maxNPtchBins; i++)
      trk_counts_pTch[i] = 0;
    for (int i = 0; i < nFinePtchBins; i++)
      trk_counts_FinePtch[i] = 0;
    for (int i = 0; i < maxNXhZBins; i++)
      trk_counts_xhZ[i] = 0;
    for (int i = 0; i < nFineXhZBins; i++)
      trk_counts_FineXhZ[i] = 0;

    // Find the next unused minimum bias event
    z_pt_mixed = z_pt;
    //z_eta_mixed = z_eta;
    z_phi_mixed = z_phi;
    //z_m_mixed = z_m;
    //
    // disable big branches during event mixing
    inTree->SetBranchStatus ("trk_pt", 0);
    inTree->SetBranchStatus ("trk_eta", 0);
    inTree->SetBranchStatus ("trk_phi", 0);
 
    bool goodMixEvent = false;
    int iMixedEvent = iEvent;
    do {
      iMixedEvent = (iMixedEvent+1) % nEvents;
      inTree->GetEntry (iMixedEvent);
      goodMixEvent = (1 < z_pt && z_pt < 12 && DeltaPhi (phi_transmin, z_phi_mixed) >= 7.*pi/8.);
    } while (!goodMixEvent && iMixedEvent != iEvent); // only check each event once
    if (iMixedEvent == iEvent || !goodMixEvent) {
      cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
      continue;
    }
 
    // reenable big branches after finding a suitable event
    inTree->SetBranchStatus ("trk_pt", 1);
    inTree->SetBranchStatus ("trk_eta", 1);
    inTree->SetBranchStatus ("trk_phi", 1);
    inTree->GetEntry (iMixedEvent);
    

    // now loop over the particles in the mixed event
    for (int iPart = 0; iPart < part_n; iPart++) {
      if (fabs (part_eta[iPart]) > 2.5)
        continue;
      if (DeltaPhi (part_phi[iPart], z_phi_mixed) < 7*TMath::Pi ()/8)
        continue;

      const float pptch = part_pt[iPart];
      const float pxhz = pptch / z_pt_mixed;

      if (pptch >= pTchBins[iPtZ][0]) {
        int iPtch = 0;
        while (iPtch < nPtchBins[iPtZ] && pptch >= pTchBins[iPtZ][iPtch+1]) iPtch++;
        if (0 <= iPtch && iPtch < nPtchBins[iPtZ])
          trk_counts_pTch[iPtch]++;
      }
      if (pptch >= FinePtchBins[iPtZ-2][0]) {
        int iPtch = 0;
        while (iPtch < nFinePtchBins && pptch >= FinePtchBins[iPtZ-2][iPtch+1]) iPtch++;
        if (0 <= iPtch && iPtch < nFinePtchBins)
          trk_counts_FinePtch[iPtch]++;
      }
      if (pxhz >= xhZBins[iPtZ][0]) {
        int iXhZ = 0;
        while (iXhZ < nXhZBins[iPtZ] && pxhz >= xhZBins[iPtZ][iXhZ+1]) iXhZ++;
        if (0 <= iXhZ && iXhZ < nXhZBins[iPtZ])
          trk_counts_xhZ[iXhZ]++;
      }
      if (pxhz >= FineXhZBins[iPtZ-2][0]) {
        int iXhZ = 0;
        while (iXhZ < nFineXhZBins && pxhz >= FineXhZBins[iPtZ-2][iXhZ+1]) iXhZ++;
        if (0 <= iXhZ && iXhZ < nFineXhZBins)
          trk_counts_FineXhZ[iXhZ]++;
      }
    } // end loop over iPart

    h = h_trk_pTch_yield_bkg[iPtZ-2];
    h2 = h2_trk_pTch_cov_bkg[iPtZ-2];
    for (short iX = 0; iX < nPtchBins[iPtZ]; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_pTch[iX]);
      for (short iY = 0; iY < nPtchBins[iPtZ]; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_pTch[iX])*(trk_counts_pTch[iY]));
    }
    h = h_trk_FinePtch_yield_bkg[iPtZ-2];
    h2 = h2_trk_FinePtch_cov_bkg[iPtZ-2];
    for (short iX = 0; iX < nFinePtchBins; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_FinePtch[iX]);
      for (short iY = 0; iY < nFinePtchBins; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_FinePtch[iX])*(trk_counts_FinePtch[iY]));
    }
    h = h_trk_xhZ_yield_bkg[iPtZ-2];
    h2 = h2_trk_xhZ_cov_bkg[iPtZ-2];
    for (short iX = 0; iX < nXhZBins[iPtZ]; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_xhZ[iX]);
      for (short iY = 0; iY < nXhZBins[iPtZ]; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_xhZ[iX])*(trk_counts_xhZ[iY]));
    }
    h = h_trk_FineXhZ_yield_bkg[iPtZ-2];
    h2 = h2_trk_FineXhZ_cov_bkg[iPtZ-2];
    for (short iX = 0; iX < nFineXhZBins; iX++) {
      h->SetBinContent (iX+1, h->GetBinContent (iX+1) + trk_counts_FineXhZ[iX]);
      for (short iY = 0; iY < nFineXhZBins; iY++)
        h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) + (trk_counts_FineXhZ[iX])*(trk_counts_FineXhZ[iY]));
    }

    for (int i = 0; i < maxNPtchBins; i++)
      trk_counts_pTch[i] = 0;
    for (int i = 0; i < nFinePtchBins; i++)
      trk_counts_FinePtch[i] = 0;
    for (int i = 0; i < maxNXhZBins; i++)
      trk_counts_xhZ[i] = 0;
    for (int i = 0; i < nFineXhZBins; i++)
      trk_counts_FineXhZ[i] = 0;
 
  } // end loop over iEvent

  inFile->Close ();
  SaferDelete (&inFile);


  /*// Post processing
  for (int iPtZ = 0; iPtZ < 3; iPtZ++) {
    const float nTotalEvents = h_z_counts->GetBinContent (iPtZ+2+1);

    h_trk_FinePtch_yield[iPtZ]->Scale (1./nTotalEvents, "width");
    h_trk_FineXhZ_yield[iPtZ]->Scale (1./nTotalEvents, "width");
    h_trk_FinePtch_yield_bkg[iPtZ]->Scale (1./nTotalEvents, "width");
    h_trk_FineXhZ_yield_bkg[iPtZ]->Scale (1./nTotalEvents, "width");

    h2_trk_FinePtch_cov[iPtZ]->Scale (1., "width");
    h2_trk_FineXhZ_cov[iPtZ]->Scale (1., "width");
    h2_trk_FinePtch_cov_bkg[iPtZ]->Scale (1., "width");
    h2_trk_FineXhZ_cov_bkg[iPtZ]->Scale (1., "width");

    TH2D* h2 = h2_trk_FinePtch_cov[iPtZ];
    TH1D* h = h_trk_FinePtch_yield[iPtZ];
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
        h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
    h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));

    h2 = h2_trk_FinePtch_cov_bkg[iPtZ];
    h = h_trk_FinePtch_yield_bkg[iPtZ];
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
        h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
    h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));

    h2 = h2_trk_FineXhZ_cov[iPtZ];
    h = h_trk_FineXhZ_yield[iPtZ];
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
        h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
    h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));

    h2 = h2_trk_FineXhZ_cov_bkg[iPtZ];
    h = h_trk_FineXhZ_yield_bkg[iPtZ];
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
        h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
    h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
    for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
      h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
  }
  */


  // now save histograms to a rootfile
  outFile->cd ();
  for (int iPtZ = 0; iPtZ < 3; iPtZ++) {
    h_trk_pTch_yield[iPtZ]->Write ();
    h_trk_xhZ_yield[iPtZ]->Write ();
    h2_trk_pTch_cov[iPtZ]->Write ();
    h2_trk_xhZ_cov[iPtZ]->Write ();
    h_trk_pTch_yield_bkg[iPtZ]->Write ();
    h_trk_xhZ_yield_bkg[iPtZ]->Write ();
    h2_trk_pTch_cov_bkg[iPtZ]->Write ();
    h2_trk_xhZ_cov_bkg[iPtZ]->Write ();
    h_trk_FinePtch_yield[iPtZ]->Write ();
    h_trk_FineXhZ_yield[iPtZ]->Write ();
    h2_trk_FinePtch_cov[iPtZ]->Write ();
    h2_trk_FineXhZ_cov[iPtZ]->Write ();
    h_trk_FinePtch_yield_bkg[iPtZ]->Write ();
    h_trk_FineXhZ_yield_bkg[iPtZ]->Write ();
    h2_trk_FinePtch_cov_bkg[iPtZ]->Write ();
    h2_trk_FineXhZ_cov_bkg[iPtZ]->Write ();
  }
  h_z_counts->Write ();
  
  outFile->Close ();

}


void CalcPythiaPerEventYields () {
  DoStuff ("pp_Zee.root", "pp_Zee_hists.root");
  DoStuff ("pp_Zmumu.root", "pp_Zmumu_hists.root");
}

#endif // __CalcPythiaPerEventYields_C__
