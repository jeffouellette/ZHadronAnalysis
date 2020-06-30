#ifndef __TagAndProbe_cxx__
#define __TagAndProbe_cxx__

#include "TagAndProbe.h"
#include "Params.h"
#include "TreeVariables.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"

#include <Utilities.h>
#include <AtlasUtils.h>

#include <TChain.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace std;

namespace ZHadronAnalysis {

bool TagAndProbe (const char* directory,
                            const int dataSet,
                            const char* inFileName) {
 
  SetupDirectories ("TagAndProbe");

  //if (isMC) {
  //  cout << "Error: In TagAndProbe.cxx: Trying to calculate trigger efficiency in MC! Quitting." << endl;
  //  return false;
  //}

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TagAndProbe.cxx: File Identifier: " << identifier << endl;

  //TChain* tree = new TChain ("bush", "bush");
  //TString pattern = "*.root";
  //auto dir = gSystem->OpenDirectory (dataPath + directory);
  //while (const char* f = gSystem->GetDirEntry (dir)) {
  //  TString file = TString (f);
  //  if (!file.Contains (Form ("%i", dataSet)))
  //    continue;
  //  cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << endl;
  //  tree->Add (dataPath + directory + "/" + file + "/*.root");
  //  break;
  //}
  //cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In TagAndProbe.C: Cannot identify this MC file! Quitting." << endl;
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
    cout << "Error: In TagAndProbe.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }

  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In TagAndProbe.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }


  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetClusters ();
  t->SetGetMuons ();
  if (isMC) {
    t->SetGetTruthElectrons ();
    t->SetGetTruthMuons ();
  }
  //t->SetGetTracks ();
  if (isPbPb) t->SetGetZdc ();
  t->SetBranchAddresses ();

  Trigger* electronTrigger = nullptr, *muonTrigger = nullptr;
  std::string electronTrigName, muonTrigName;
  if (collisionSystem == PbPb15) {
    electronTrigName = "HLT_e15_loose_L1EM12";
    muonTrigName = "HLT_mu14";
  }
  else if (collisionSystem == pp17) {
    electronTrigName = "HLT_e15_lhloose_L1EM12";
    muonTrigName = "HLT_mu14";
  }
  else if (collisionSystem == PbPb18) {
    electronTrigName = "HLT_e15_lhloose_ion_L1EM12";
    muonTrigName = "HLT_mu14";
  }
  else {
    cout << "Error: In TagAndProbe.cxx: Invalid collision system, quitting." << endl;
    return false;
  }
  if (!isMC) {
    for (int iTrig = 0; iTrig < nElectronTrig; iTrig++) {
      if (electronTrigNames[iTrig] != electronTrigName)
        continue;
      Trigger* trig = new Trigger (electronTrigNames[iTrig].Data (), electronTrigMinPtCuts[iTrig], electronTrigMinEtaCuts[iTrig], electronTrigMaxEtaCuts[iTrig]);
      trig->minPt = electronTrigMinPtCuts[iTrig];
      trig->maxPt = electronTrigMaxPtCuts[iTrig];
      tree->SetBranchAddress (electronTrigNames[iTrig], &(trig->trigBool));
      tree->SetBranchAddress (electronTrigNames[iTrig]+"_prescale", &(trig->trigPrescale));
      electronTrigger = trig;
    }
    for (int iTrig = 0; iTrig < nMuonTrig; iTrig++) {
      if (muonTrigNames[iTrig] != muonTrigName)
        continue;
      Trigger* trig = new Trigger (muonTrigNames[iTrig].Data (), muonTrigMinPtCuts[iTrig], muonTrigMinEtaCuts[iTrig], muonTrigMaxEtaCuts[iTrig]);
      trig->minPt = muonTrigMinPtCuts[iTrig];
      trig->maxPt = muonTrigMaxPtCuts[iTrig];
      tree->SetBranchAddress (muonTrigNames[iTrig], &(trig->trigBool));
      tree->SetBranchAddress (muonTrigNames[iTrig]+"_prescale", &(trig->trigPrescale));
      muonTrigger = trig;
    }
    if (!electronTrigger || !muonTrigger) {
      cout << "Error: In TagAndProbe.cxx: Failed to define triggers, quitting!"<< endl;
      return false;
    }
  }

  double fcal_et = 0, vz = 0;

  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  TH1D* h_muonTrigEffNum_pt[3];
  TH1D* h_muonTrigEffDen_pt[3];
  TH2D* h2_muonTrigEffNum_eta_phi[3];
  TH2D* h2_muonTrigEffDen_eta_phi[3];
  TH1D* h_muonTrigEffNum_fcal[2];
  TH1D* h_muonTrigEffDen_fcal[2];

  TH1D* h_muonIDEffNum_pt[3];
  TH1D* h_muonIDEffDen_pt[3];
  TH2D* h2_muonIDEffNum_eta_phi[3];
  TH2D* h2_muonIDEffDen_eta_phi[3];
  TH1D* h_muonIDEffNum_fcal[2];
  TH1D* h_muonIDEffDen_fcal[2];

  TH2D* h2_zmumuTrigEffNum_pt_y[3];
  TH2D* h2_zmumuTrigEffDen_pt_y[3];
  TH2D* h2_zmumuIDEffNum_pt_y[3];
  TH2D* h2_zmumuIDEffDen_pt_y[3];

  TH1D* h_electronTrigEffNum_pt[3];
  TH1D* h_electronTrigEffDen_pt[3];
  TH1D* h_electronTrigEffNum_eta[3];
  TH1D* h_electronTrigEffDen_eta[3];
  TH2D* h2_electronTrigEffNum_pt_eta[3];
  TH2D* h2_electronTrigEffDen_pt_eta[3];
  TH1D* h_electronTrigEffNum_fcal[2];
  TH1D* h_electronTrigEffDen_fcal[2];

  TH1D* h_electronIDEffNum_pt[3];
  TH1D* h_electronIDEffDen_pt[3];
  TH1D* h_electronIDEffNum_eta[3];
  TH1D* h_electronIDEffDen_eta[3];
  TH2D* h2_electronIDEffNum_pt_eta[3];
  TH2D* h2_electronIDEffDen_pt_eta[3];
  TH1D* h_electronIDEffNum_fcal[2];
  TH1D* h_electronIDEffDen_fcal[2];

  TH2D* h2_zeeTrigEffNum_pt_y[3];
  TH2D* h2_zeeTrigEffDen_pt_y[3];
  TH2D* h2_zeeIDEffNum_pt_y[3];
  TH2D* h2_zeeIDEffDen_pt_y[3];

  double zPtBins[17] = {0, 3, 6, 10, 15, 20, 25, 30, 40, 50, 70, 90, 110, 140, 170, 200, 240};

  h_muonTrigEffNum_pt[0] = new TH1D ("h_muonTrigEffNum_pt_pp", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonTrigEffDen_pt[0] = new TH1D ("h_muonTrigEffDen_pt_pp", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonTrigEffNum_pt[1] = new TH1D ("h_muonTrigEffNum_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonTrigEffDen_pt[1] = new TH1D ("h_muonTrigEffDen_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonTrigEffNum_pt[2] = new TH1D ("h_muonTrigEffNum_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonTrigEffDen_pt[2] = new TH1D ("h_muonTrigEffDen_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonTrigEffNum_pt[0]->Sumw2 ();
  h_muonTrigEffDen_pt[0]->Sumw2 ();
  h_muonTrigEffNum_pt[1]->Sumw2 ();
  h_muonTrigEffDen_pt[1]->Sumw2 ();
  h_muonTrigEffNum_pt[2]->Sumw2 ();
  h_muonTrigEffDen_pt[2]->Sumw2 ();
  h2_muonTrigEffNum_eta_phi[0] = new TH2D ("h2_muonTrigEffNum_eta_phi_pp", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonTrigEffDen_eta_phi[0] = new TH2D ("h2_muonTrigEffDen_eta_phi_pp", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonTrigEffNum_eta_phi[1] = new TH2D ("h2_muonTrigEffNum_eta_phi_PbPb18", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonTrigEffDen_eta_phi[1] = new TH2D ("h2_muonTrigEffDen_eta_phi_PbPb18", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonTrigEffNum_eta_phi[2] = new TH2D ("h2_muonTrigEffNum_eta_phi_PbPb15", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonTrigEffDen_eta_phi[2] = new TH2D ("h2_muonTrigEffDen_eta_phi_PbPb15", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonTrigEffNum_eta_phi[0]->Sumw2 ();
  h2_muonTrigEffDen_eta_phi[0]->Sumw2 ();
  h2_muonTrigEffNum_eta_phi[1]->Sumw2 ();
  h2_muonTrigEffDen_eta_phi[1]->Sumw2 ();
  h2_muonTrigEffNum_eta_phi[2]->Sumw2 ();
  h2_muonTrigEffDen_eta_phi[2]->Sumw2 ();
  h_muonTrigEffNum_fcal[0] = new TH1D ("h_muonTrigEffNum_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonTrigEffDen_fcal[0] = new TH1D ("h_muonTrigEffDen_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonTrigEffNum_fcal[1] = new TH1D ("h_muonTrigEffNum_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonTrigEffDen_fcal[1] = new TH1D ("h_muonTrigEffDen_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonTrigEffNum_fcal[0]->Sumw2 ();
  h_muonTrigEffDen_fcal[0]->Sumw2 ();
  h_muonTrigEffNum_fcal[1]->Sumw2 ();
  h_muonTrigEffDen_fcal[1]->Sumw2 ();

  h_muonIDEffNum_pt[0] = new TH1D ("h_muonIDEffNum_pt_pp", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonIDEffDen_pt[0] = new TH1D ("h_muonIDEffDen_pt_pp", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonIDEffNum_pt[1] = new TH1D ("h_muonIDEffNum_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonIDEffDen_pt[1] = new TH1D ("h_muonIDEffDen_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonIDEffNum_pt[2] = new TH1D ("h_muonIDEffNum_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonIDEffDen_pt[2] = new TH1D ("h_muonIDEffDen_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 25, 10, 85);
  h_muonIDEffNum_pt[0]->Sumw2 ();
  h_muonIDEffDen_pt[0]->Sumw2 ();
  h_muonIDEffNum_pt[1]->Sumw2 ();
  h_muonIDEffDen_pt[1]->Sumw2 ();
  h_muonIDEffNum_pt[2]->Sumw2 ();
  h_muonIDEffDen_pt[2]->Sumw2 ();
  h2_muonIDEffNum_eta_phi[0] = new TH2D ("h2_muonIDEffNum_eta_phi_pp", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonIDEffDen_eta_phi[0] = new TH2D ("h2_muonIDEffDen_eta_phi_pp", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonIDEffNum_eta_phi[1] = new TH2D ("h2_muonIDEffNum_eta_phi_PbPb18", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonIDEffDen_eta_phi[1] = new TH2D ("h2_muonIDEffDen_eta_phi_PbPb18", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonIDEffNum_eta_phi[2] = new TH2D ("h2_muonIDEffNum_eta_phi_PbPb15", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonIDEffDen_eta_phi[2] = new TH2D ("h2_muonIDEffDen_eta_phi_PbPb15", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  h2_muonIDEffNum_eta_phi[0]->Sumw2 ();
  h2_muonIDEffDen_eta_phi[0]->Sumw2 ();
  h2_muonIDEffNum_eta_phi[1]->Sumw2 ();
  h2_muonIDEffDen_eta_phi[1]->Sumw2 ();
  h2_muonIDEffNum_eta_phi[2]->Sumw2 ();
  h2_muonIDEffDen_eta_phi[2]->Sumw2 ();
  h_muonIDEffNum_fcal[0] = new TH1D ("h_muonIDEffNum_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonIDEffDen_fcal[0] = new TH1D ("h_muonIDEffDen_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonIDEffNum_fcal[1] = new TH1D ("h_muonIDEffNum_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonIDEffDen_fcal[1] = new TH1D ("h_muonIDEffDen_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_muonIDEffNum_fcal[0]->Sumw2 ();
  h_muonIDEffDen_fcal[0]->Sumw2 ();
  h_muonIDEffNum_fcal[1]->Sumw2 ();
  h_muonIDEffDen_fcal[1]->Sumw2 ();

  h2_zmumuTrigEffNum_pt_y[0] = new TH2D ("h2_zmumuTrigEffNum_pt_y_pp",     ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuTrigEffDen_pt_y[0] = new TH2D ("h2_zmumuTrigEffDen_pt_y_pp",     ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuTrigEffNum_pt_y[1] = new TH2D ("h2_zmumuTrigEffNum_pt_y_PbPb18", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuTrigEffDen_pt_y[1] = new TH2D ("h2_zmumuTrigEffDen_pt_y_PbPb18", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuTrigEffNum_pt_y[2] = new TH2D ("h2_zmumuTrigEffNum_pt_y_PbPb15", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuTrigEffDen_pt_y[2] = new TH2D ("h2_zmumuTrigEffDen_pt_y_PbPb15", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuIDEffNum_pt_y[0]   = new TH2D ("h2_zmumuIDEffNum_pt_y_pp",       ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuIDEffDen_pt_y[0]   = new TH2D ("h2_zmumuIDEffDen_pt_y_pp",       ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuIDEffNum_pt_y[1]   = new TH2D ("h2_zmumuIDEffNum_pt_y_PbPb18",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuIDEffDen_pt_y[1]   = new TH2D ("h2_zmumuIDEffDen_pt_y_PbPb18",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuIDEffNum_pt_y[2]   = new TH2D ("h2_zmumuIDEffNum_pt_y_PbPb15",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zmumuIDEffDen_pt_y[2]   = new TH2D ("h2_zmumuIDEffDen_pt_y_PbPb15",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);

  const float electronPtBins[17] = {10, 12, 14, 16, 18, 20, 22, 24, 27, 30, 35, 40, 45, 50, 60, 70, 100};
  const float electronEtaBins[21] = {-2.47, -2.37, -2.01, -1.81, -1.52, -1.37, -1.15, -0.8, -0.6, -0.1, 0, 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47};
  h_electronTrigEffNum_pt[0] = new TH1D ("h_electronTrigEffNum_pt_pp", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronTrigEffDen_pt[0] = new TH1D ("h_electronTrigEffDen_pt_pp", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronTrigEffNum_pt[1] = new TH1D ("h_electronTrigEffNum_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronTrigEffDen_pt[1] = new TH1D ("h_electronTrigEffDen_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronTrigEffNum_pt[2] = new TH1D ("h_electronTrigEffNum_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronTrigEffDen_pt[2] = new TH1D ("h_electronTrigEffDen_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronTrigEffNum_pt[0]->Sumw2 ();
  h_electronTrigEffDen_pt[0]->Sumw2 ();
  h_electronTrigEffNum_pt[1]->Sumw2 ();
  h_electronTrigEffDen_pt[1]->Sumw2 ();
  h_electronTrigEffNum_pt[2]->Sumw2 ();
  h_electronTrigEffDen_pt[2]->Sumw2 ();
  h_electronTrigEffNum_eta[0] = new TH1D ("h_electronTrigEffNum_eta_pp", ";#eta;Counts", 20, electronEtaBins);
  h_electronTrigEffDen_eta[0] = new TH1D ("h_electronTrigEffDen_eta_pp", ";#eta;Counts", 20, electronEtaBins);
  h_electronTrigEffNum_eta[1] = new TH1D ("h_electronTrigEffNum_eta_PbPb18", ";#eta;Counts", 20, electronEtaBins);
  h_electronTrigEffDen_eta[1] = new TH1D ("h_electronTrigEffDen_eta_PbPb18", ";#eta;Counts", 20, electronEtaBins);
  h_electronTrigEffNum_eta[2] = new TH1D ("h_electronTrigEffNum_eta_PbPb15", ";#eta;Counts", 20, electronEtaBins);
  h_electronTrigEffDen_eta[2] = new TH1D ("h_electronTrigEffDen_eta_PbPb15", ";#eta;Counts", 20, electronEtaBins);
  h_electronTrigEffNum_eta[0]->Sumw2 ();
  h_electronTrigEffDen_eta[0]->Sumw2 ();
  h_electronTrigEffNum_eta[1]->Sumw2 ();
  h_electronTrigEffDen_eta[1]->Sumw2 ();
  h_electronTrigEffNum_eta[2]->Sumw2 ();
  h_electronTrigEffDen_eta[2]->Sumw2 ();
  h2_electronTrigEffNum_pt_eta[0] = new TH2D ("h2_electronTrigEffNum_pt_eta_pp", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronTrigEffDen_pt_eta[0] = new TH2D ("h2_electronTrigEffDen_pt_eta_pp", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronTrigEffNum_pt_eta[1] = new TH2D ("h2_electronTrigEffNum_pt_eta_PbPb18", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronTrigEffDen_pt_eta[1] = new TH2D ("h2_electronTrigEffDen_pt_eta_PbPb18", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronTrigEffNum_pt_eta[2] = new TH2D ("h2_electronTrigEffNum_pt_eta_PbPb15", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronTrigEffDen_pt_eta[2] = new TH2D ("h2_electronTrigEffDen_pt_eta_PbPb15", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronTrigEffNum_pt_eta[0]->Sumw2 ();
  h2_electronTrigEffDen_pt_eta[0]->Sumw2 ();
  h2_electronTrigEffNum_pt_eta[1]->Sumw2 ();
  h2_electronTrigEffDen_pt_eta[1]->Sumw2 ();
  h2_electronTrigEffNum_pt_eta[2]->Sumw2 ();
  h2_electronTrigEffDen_pt_eta[2]->Sumw2 ();
  h_electronTrigEffNum_fcal[0] = new TH1D ("h_electronTrigEffNum_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronTrigEffDen_fcal[0] = new TH1D ("h_electronTrigEffDen_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronTrigEffNum_fcal[1] = new TH1D ("h_electronTrigEffNum_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronTrigEffDen_fcal[1] = new TH1D ("h_electronTrigEffDen_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronTrigEffNum_fcal[0]->Sumw2 ();
  h_electronTrigEffDen_fcal[0]->Sumw2 ();
  h_electronTrigEffNum_fcal[1]->Sumw2 ();
  h_electronTrigEffDen_fcal[1]->Sumw2 ();

  //h_electronIDEffNum_m[0] = new TH1D ("h_electronIDEffNum_m_pp", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffDen_m[0] = new TH1D ("h_electronIDEffDen_m_pp", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffNum_m[1] = new TH1D ("h_electronIDEffNum_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffDen_m[1] = new TH1D ("h_electronIDEffDen_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffNum_m[2] = new TH1D ("h_electronIDEffNum_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffDen_m[2] = new TH1D ("h_electronIDEffDen_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffBkgNum_m[0] = new TH1D ("h_electronIDEffBkgNum_m_pp", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffBkgDen_m[0] = new TH1D ("h_electronIDEffBkgDen_m_pp", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffBkgNum_m[1] = new TH1D ("h_electronIDEffBkgNum_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffBkgDen_m[1] = new TH1D ("h_electronIDEffBkgDen_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffBkgNum_m[2] = new TH1D ("h_electronIDEffBkgNum_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);
  //h_electronIDEffBkgDen_m[2] = new TH1D ("h_electronIDEffBkgDen_m_PbPb18", ";#it{m}_{ee} [GeV];Counts", 80, 50, 130);

  h_electronIDEffNum_pt[0] = new TH1D ("h_electronIDEffNum_pt_pp", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronIDEffDen_pt[0] = new TH1D ("h_electronIDEffDen_pt_pp", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronIDEffNum_pt[1] = new TH1D ("h_electronIDEffNum_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronIDEffDen_pt[1] = new TH1D ("h_electronIDEffDen_pt_PbPb18", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronIDEffNum_pt[2] = new TH1D ("h_electronIDEffNum_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronIDEffDen_pt[2] = new TH1D ("h_electronIDEffDen_pt_PbPb15", ";#it{p}_{T} [GeV];Counts", 16, electronPtBins);
  h_electronIDEffNum_pt[0]->Sumw2 ();
  h_electronIDEffDen_pt[0]->Sumw2 ();
  h_electronIDEffNum_pt[1]->Sumw2 ();
  h_electronIDEffDen_pt[1]->Sumw2 ();
  h_electronIDEffNum_pt[2]->Sumw2 ();
  h_electronIDEffDen_pt[2]->Sumw2 ();
  h_electronIDEffNum_eta[0] = new TH1D ("h_electronIDEffNum_eta_pp", ";#eta;Counts", 20, electronEtaBins);
  h_electronIDEffDen_eta[0] = new TH1D ("h_electronIDEffDen_eta_pp", ";#eta;Counts", 20, electronEtaBins);
  h_electronIDEffNum_eta[1] = new TH1D ("h_electronIDEffNum_eta_PbPb18", ";#eta;Counts", 20, electronEtaBins);
  h_electronIDEffDen_eta[1] = new TH1D ("h_electronIDEffDen_eta_PbPb18", ";#eta;Counts", 20, electronEtaBins);
  h_electronIDEffNum_eta[2] = new TH1D ("h_electronIDEffNum_eta_PbPb15", ";#eta;Counts", 20, electronEtaBins);
  h_electronIDEffDen_eta[2] = new TH1D ("h_electronIDEffDen_eta_PbPb15", ";#eta;Counts", 20, electronEtaBins);
  h_electronIDEffNum_eta[0]->Sumw2 ();
  h_electronIDEffDen_eta[0]->Sumw2 ();
  h_electronIDEffNum_eta[1]->Sumw2 ();
  h_electronIDEffDen_eta[1]->Sumw2 ();
  h_electronIDEffNum_eta[2]->Sumw2 ();
  h_electronIDEffDen_eta[2]->Sumw2 ();
  h2_electronIDEffNum_pt_eta[0] = new TH2D ("h2_electronIDEffNum_pt_eta_pp", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronIDEffDen_pt_eta[0] = new TH2D ("h2_electronIDEffDen_pt_eta_pp", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronIDEffNum_pt_eta[1] = new TH2D ("h2_electronIDEffNum_pt_eta_PbPb18", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronIDEffDen_pt_eta[1] = new TH2D ("h2_electronIDEffDen_pt_eta_PbPb18", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronIDEffNum_pt_eta[2] = new TH2D ("h2_electronIDEffNum_pt_eta_PbPb15", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronIDEffDen_pt_eta[2] = new TH2D ("h2_electronIDEffDen_pt_eta_PbPb15", ";#eta;#it{p}_{T} [GeV];Counts", 20, electronEtaBins, 16, electronPtBins);
  h2_electronIDEffNum_pt_eta[0]->Sumw2 ();
  h2_electronIDEffDen_pt_eta[0]->Sumw2 ();
  h2_electronIDEffNum_pt_eta[1]->Sumw2 ();
  h2_electronIDEffDen_pt_eta[1]->Sumw2 ();
  h2_electronIDEffNum_pt_eta[2]->Sumw2 ();
  h2_electronIDEffDen_pt_eta[2]->Sumw2 ();
  h_electronIDEffNum_fcal[0] = new TH1D ("h_electronIDEffNum_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronIDEffDen_fcal[0] = new TH1D ("h_electronIDEffDen_fcal_PbPb18", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronIDEffNum_fcal[1] = new TH1D ("h_electronIDEffNum_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronIDEffDen_fcal[1] = new TH1D ("h_electronIDEffDen_fcal_PbPb15", ";#Sigma#it{E}_{T}^{FCal} [TeV];Counts", 11, -0.5, 5);
  h_electronIDEffNum_fcal[0]->Sumw2 ();
  h_electronIDEffDen_fcal[0]->Sumw2 ();
  h_electronIDEffNum_fcal[1]->Sumw2 ();
  h_electronIDEffDen_fcal[1]->Sumw2 ();

  h2_zeeTrigEffNum_pt_y[0] = new TH2D ("h2_zeeTrigEffNum_pt_y_pp",     ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeTrigEffDen_pt_y[0] = new TH2D ("h2_zeeTrigEffDen_pt_y_pp",     ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeTrigEffNum_pt_y[1] = new TH2D ("h2_zeeTrigEffNum_pt_y_PbPb18", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeTrigEffDen_pt_y[1] = new TH2D ("h2_zeeTrigEffDen_pt_y_PbPb18", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeTrigEffNum_pt_y[2] = new TH2D ("h2_zeeTrigEffNum_pt_y_PbPb15", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeTrigEffDen_pt_y[2] = new TH2D ("h2_zeeTrigEffDen_pt_y_PbPb15", ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeIDEffNum_pt_y[0]   = new TH2D ("h2_zeeIDEffNum_pt_y_pp",       ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeIDEffDen_pt_y[0]   = new TH2D ("h2_zeeIDEffDen_pt_y_pp",       ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeIDEffNum_pt_y[1]   = new TH2D ("h2_zeeIDEffNum_pt_y_PbPb18",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeIDEffDen_pt_y[1]   = new TH2D ("h2_zeeIDEffDen_pt_y_PbPb18",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeIDEffNum_pt_y[2]   = new TH2D ("h2_zeeIDEffNum_pt_y_PbPb15",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);
  h2_zeeIDEffDen_pt_y[2]   = new TH2D ("h2_zeeIDEffDen_pt_y_PbPb15",   ";y_{Z};#it{p}_{T}^{Z} [GeV];Counts", 10, -2.5, 2.5, 16, zPtBins);


  const int nEvts = tree->GetEntries ();


  ///////////////////////////////////////////////////////////////////////
  ////**** Calculate background weights for muon Med|ID efficiency ****//
  ///////////////////////////////////////////////////////////////////////

  //// Loop over events
  //for (int iEvt = 0; iEvt < nEvts; iEvt++) {
  //  if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
  //    cout << "Info: In TagAndProbe.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
  //  tree->GetEntry (iEvt);

  //  //**** Event selection ****//
  //  if (collisionSystem == PbPb18 && t->BlayerDesyn)
  //    continue;
  //  if (isPbPb && t->isOOTPU)
  //    continue;

  //  bool hasPrimaryVert = false;
  //  //bool hasPileupVert = false;
  //  for (int iVert = 0; iVert < t->nvert; iVert++) {
  //    if (t->vert_type[iVert] == 1) {
  //      if (hasPrimaryVert) {
  //        hasPrimaryVert = false;
  //        break;
  //      }
  //      hasPrimaryVert = (t->vert_type[iVert] == 1);
  //      vz = t->vert_z[iVert];
  //    }
  //    //else if (t->vert_type[iVert] == 3) {
  //    //  hasPileupVert = true;
  //    //}
  //  }
  //  if (!hasPrimaryVert)// || hasPileupVert)
  //    continue;

  //  if (isPbPb) {
  //    fcal_et = t->fcalA_et + t->fcalC_et;

  //    const float zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
  //    const float nNeutrons = (zdcEnergy) / (2.51);
  //    const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3);
  //    if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
  //      continue; // Zdc-based in-time pile-up cut
  //  }

  //  //**** Loop over tag muons ****//
  //  if (muonTrigger->trigBool && t->passes_toroid) {
  //    for (int iTag = 0; iTag < t->muon_n; iTag++) {
  //      if (!isMC && !t->muon_matched[iTag])
  //        continue;
  //      if (!t->muon_tight[iTag])
  //        continue;
  //      if (fabs (t->muon_eta[iTag]) > 2.5)
  //        continue;
  //      if (t->muon_pt[iTag] < 30)
  //        continue;

  //      TLorentzVector mtag;
  //      mtag.SetPtEtaPhiM (t->muon_pt[iTag], t->muon_eta[iTag], t->muon_phi[iTag], muon_mass);

  //      //**** Calculate background on muon medium efficiency given an ID track hists. ****//
  //      for (int iProbe = 0; iProbe < t->muon_n; iProbe++) {
  //        if (iProbe == iTag)
  //          continue;

  //        if (t->muon_id_track_pt[iProbe] != 0)
  //          continue;
  //        if (t->muon_id_track_d0sig[iProbe] > 3)
  //          continue;
  //        if (fabs (t->muon_eta[iProbe]) > 2.5)
  //          continue;
  //        if (DeltaR (t->muon_eta[iTag], t->muon_eta[iProbe], t->muon_phi[iTag], t->muon_phi[iProbe]) < 0.3)
  //          continue; // avoids ambiguity of muons being too close together
  //        if (!GetMuonTrackCut (t, iProbe, "HILoose"))
  //          continue;

  //        TLorentzVector mprobe;
  //        mprobe.SetPtEtaPhiM (t->muon_pt[iProbe], t->muon_eta[iProbe], t->muon_phi[iProbe], muon_mass);

  //        int iSign = (t->muon_charge[iProbe] == t->muon_charge[iTag] ? 0 : 1);

  //        TLorentzVector z = mtag+mprobe;
  //        if (z.M () < 76 || 106 < z.M ())
  //          continue;

  //        if (!isPbPb) {
  //          h_muonMed_IDEffDen_m[0][iSign]->Fill (z.M ());
  //          if (t->muon_medium[iProbe])
  //            h_muonMed_IDEffNum_m[0][iSign]->Fill (z.M ());
  //        }
  //        else if (isPbPb && !is2015data) {
  //          h_muonMed_IDEffDen_m[1][iSign]->Fill (z.M ());
  //          if (t->muon_medium[iProbe])
  //            h_muonMed_IDEffNum_m[1][iSign]->Fill (z.M ());
  //        }
  //        else if (isPbPb && is2015data) {
  //          h_muonMed_IDEffDen_m[2][iSign]->Fill (z.M ());
  //          if (t->muon_medium[iProbe])
  //            h_muonMed_IDEffNum_m[2][iSign]->Fill (z.M ());
  //        }

  //      } // end loop over probe

  //      //**** Calculate background on muon ID track efficiency given an MS track hists. ****//
  //      for (int iProbe = 0; iProbe < t->muon_n; iProbe++) {
  //        if (iProbe == iTag)
  //          continue;

  //        //if (t->muon_id_track_d0sig[iProbe] > 3)
  //        //  continue;
  //        if (fabs (t->muon_eta[iProbe]) > 2.5)
  //          continue;
  //        if (DeltaR (t->muon_eta[iTag], t->muon_eta[iProbe], t->muon_phi[iTag], t->muon_phi[iProbe]) < 0.3)
  //          continue; // avoids ambiguity of muons being too close together
  //        if (!GetMuonTrackCut (t, iProbe, "HILoose"))
  //          continue;

  //        TLorentzVector mprobe;
  //        mprobe.SetPtEtaPhiM (t->muon_pt[iProbe], t->muon_eta[iProbe], t->muon_phi[iProbe], muon_mass);

  //        int iSign = (t->muon_charge[iProbe] == t->muon_charge[iTag] ? 0 : 1);

  //        TLorentzVector z = mtag+mprobe;
  //        if (z.M () < 76 || 106 < z.M ())
  //          continue;

  //        if (!isPbPb) {
  //          h_muonID_MSEffDen_m[0][iSign]->Fill (z.M ());
  //          if (t->muon_id_track_pt[iProbe] > 0.)
  //            h_muonID_MSEffNum_m[0][iSign]->Fill (z.M ());
  //        }
  //        else if (isPbPb && !is2015data) {
  //          h_muonID_MSEffDen_m[1][iSign]->Fill (z.M ());
  //          if (t->muon_id_track_pt[iProbe] > 0.)
  //            h_muonID_MSEffNum_m[1][iSign]->Fill (z.M ());
  //        }
  //        else if (isPbPb && is2015data) {
  //          h_muonID_MSEffDen_m[2][iSign]->Fill (z.M ());
  //          if (t->muon_id_track_pt[iProbe] > 0.)
  //            h_muonID_MSEffNum_m[2][iSign]->Fill (z.M (0);
  //        }

  //      } // end loop over probe
  //    } // end loop over tag
  //  } // end if muon trigger fired (and passes_toroid)
  //} // end loop over events



  //////////////////////////////////////////////////////////////////////////////
  //**** Calculate single lepton trigger efficiencies with tag-and-probe. ****//
  //////////////////////////////////////////////////////////////////////////////

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TagAndProbe.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    //**** Event selection ****//
    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;

    bool hasPrimaryVert = false;
    //bool hasPileupVert = false;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      if (t->vert_type[iVert] == 1) {
        if (hasPrimaryVert) {
          hasPrimaryVert = false;
          break;
        }
        hasPrimaryVert = (t->vert_type[iVert] == 1);
        vz = t->vert_z[iVert];
      }
      //else if (t->vert_type[iVert] == 3) {
      //  hasPileupVert = true;
      //}
    }
    if (!hasPrimaryVert)// || hasPileupVert)
      continue;

    if (isPbPb) {
      fcal_et = t->fcalA_et + t->fcalC_et;

      const float zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3);
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

    //**** Loop over tag muons ****//
    if ((isMC || muonTrigger->trigBool) && t->passes_toroid) {
      for (int iTag = 0; iTag < t->muon_n; iTag++) {
        if (!isMC && !t->muon_matched[iTag])
          continue;
        if (!t->muon_tight[iTag])
          continue;
        if (fabs (t->muon_eta[iTag]) > 2.5)
          continue;
        if (t->muon_pt[iTag] < 30)
          continue;

        if (isMC) {
          bool hasTruthMatch = false;
          for (int iTM = 0; iTM < t->truth_muon_n; iTM++) {
            if (DeltaR (t->muon_eta[iTag], t->truth_muon_eta[iTM], t->muon_phi[iTag], t->truth_muon_phi[iTM]) < 0.3) {
              hasTruthMatch = true;
              break;
            }
          }
          if (!hasTruthMatch)
            continue;
        }

        TLorentzVector mtag;
        mtag.SetPtEtaPhiM (t->muon_pt[iTag], t->muon_eta[iTag], t->muon_phi[iTag], muon_mass);

        //**** Fill muon trigger efficiency hists. ****//
        for (int iProbe = 0; iProbe < t->muon_n; iProbe++) {
          if (iProbe == iTag)
            continue;

          if (!t->muon_medium[iProbe])
            continue;
          if (fabs (t->muon_eta[iProbe]) > 2.5)
            continue;
          if (DeltaR (t->muon_eta[iTag], t->muon_eta[iProbe], t->muon_phi[iTag], t->muon_phi[iProbe]) < 0.3)
            continue; // avoids ambiguity of muons being too close together

          TLorentzVector mprobe;
          mprobe.SetPtEtaPhiM (t->muon_pt[iProbe], t->muon_eta[iProbe], t->muon_phi[iProbe], muon_mass);

          if (t->muon_charge[iProbe] == t->muon_charge[iTag])
            continue;

          TLorentzVector z = mtag+mprobe;
          if (z.M () < 76 || 106 < z.M ())
            continue;

          if (!isPbPb) {
            h_muonTrigEffDen_pt[0]->Fill (t->muon_pt[iProbe]);
            if (t->muon_pt[iProbe] > 20) {
              h2_muonTrigEffDen_eta_phi[0]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
            }
            if (isMC || t->muon_matched[iProbe]) {
              h_muonTrigEffNum_pt[0]->Fill (t->muon_pt[iProbe]);
              if (t->muon_pt[iProbe] > 20) {
                h2_muonTrigEffNum_eta_phi[0]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
              }
            }
          }
          else if (isPbPb && !is2015data) {
            h_muonTrigEffDen_pt[1]->Fill (t->muon_pt[iProbe]);
            if (t->muon_pt[iProbe] > 20) {
              h2_muonTrigEffDen_eta_phi[1]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
              h_muonTrigEffDen_fcal[0]->Fill (fcal_et*1e-3);
            }
            if (isMC || t->muon_matched[iProbe]) {
              h_muonTrigEffNum_pt[1]->Fill (t->muon_pt[iProbe]);
              if (t->muon_pt[iProbe] > 20) {
                h2_muonTrigEffNum_eta_phi[1]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
                h_muonTrigEffNum_fcal[0]->Fill (fcal_et*1e-3);
              }
            }
          }
          else if (isPbPb && is2015data) {
            h_muonTrigEffDen_pt[2]->Fill (t->muon_pt[iProbe]);
            if (t->muon_pt[iProbe] > 20) {
              h2_muonTrigEffDen_eta_phi[2]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
              h_muonTrigEffDen_fcal[1]->Fill (fcal_et*1e-3);
            }
            if (isMC || t->muon_matched[iProbe]) {
              h_muonTrigEffNum_pt[2]->Fill (t->muon_pt[iProbe]);
              if (t->muon_pt[iProbe] > 20) {
                h2_muonTrigEffNum_eta_phi[2]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
                h_muonTrigEffNum_fcal[1]->Fill (fcal_et*1e-3);
              }
            }
          }
        } // end loop over probe


        ////**** Fill muon medium efficiency given an ID track hists. ****//
        //for (int iProbe = 0; iProbe < t->muon_n; iProbe++) {
        //  if (iProbe == iTag)
        //    continue;

        //  if (t->muon_id_track_pt[iProbe] <= 0)
        //    continue;
        //  if (fabs (t->muon_id_track_d0sig[iProbe]) > 3)
        //    continue;
        //  if (fabs (t->muon_eta[iProbe]) > 2.5)
        //    continue;
        //  if (DeltaR (t->muon_eta[iTag], t->muon_eta[iProbe], t->muon_phi[iTag], t->muon_phi[iProbe]) < 0.3)
        //    continue; // avoids ambiguity of muons being too close together
        //  if (!GetMuonTrackCut (t, iProbe, "HILoose"))
        //    continue;

        //  TLorentzVector mprobe;
        //  mprobe.SetPtEtaPhiM (t->muon_pt[iProbe], t->muon_eta[iProbe], t->muon_phi[iProbe], muon_mass);

        //  int iSign = (t->muon_charge[iProbe] == t->muon_charge[iTag] ? 0 : 1);

        //  TLorentzVector z = mtag+mprobe;
        //  if (z.M () < 76 || 106 < z.M ())
        //    continue;

        //  if (!isPbPb) {
        //    h_muonMed_IDEffDen_pt[0][iSign]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h_muonMed_IDEffDen_eta[0][iSign]->Fill (t->muon_eta[iProbe]);
        //    }
        //    if (t->muon_medium[iProbe]) {
        //      h_muonMed_IDEffNum_pt[0][iSign]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h_muonMed_IDEffNum_eta[0][iSign]->Fill (t->muon_eta[iProbe]);
        //      }
        //    }
        //  }
        //  else if (isPbPb && !is2015data) {
        //    h_muonMed_IDEffDen_pt[1][iSign]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h_muonMed_IDEffDen_eta[1][iSign]->Fill (t->muon_eta[iProbe]);
        //      h_muonMed_IDEffDen_fcal[0][iSign]->Fill (fcal_et);
        //    }
        //    if (t->muon_medium[iProbe]) {
        //      h_muonMed_IDEffNum_pt[1][iSign]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h_muonMed_IDEffNum_eta[1][iSign]->Fill (t->muon_eta[iProbe]);
        //        h_muonMed_IDEffNum_fcal[0][iSign]->Fill (fcal_et);
        //      }
        //    }
        //  }
        //  else if (isPbPb && is2015data) {
        //    h_muonMed_IDEffDen_pt[2][iSign]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h_muonMed_IDEffDen_eta[2][iSign]->Fill (t->muon_eta[iProbe]);
        //      h_muonMed_IDEffDen_fcal[1][iSign]->Fill (fcal_et);
        //    }
        //    if (t->muon_medium[iProbe]) {
        //      h_muonMed_IDEffNum_pt[2][iSign]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h_muonMed_IDEffNum_eta[2][iSign]->Fill (t->muon_eta[iProbe]);
        //        h_muonMed_IDEffNum_fcal[1][iSign]->Fill (fcal_et);
        //      }
        //    }
        //  }
        //} // end loop over probe


        ////**** Fill muon ID track efficiency given an MS track hists. ****//
        //for (int iProbe = 0; iProbe < t->muon_n; iProbe++) {
        //  if (iProbe == iTag)
        //    continue;

        //  if (t->muon_ms_pt[iProbe] <= 0)
        //    continue;
        //  //if (t->muon_id_track_d0sig[iProbe] > 3)
        //  //  continue;
        //  if (fabs (t->muon_eta[iProbe]) > 2.5)
        //    continue;
        //  if (DeltaR (t->muon_eta[iTag], t->muon_eta[iProbe], t->muon_phi[iTag], t->muon_phi[iProbe]) < 0.3)
        //    continue; // avoids ambiguity of muons being too close together
        //  if (!GetMuonTrackCut (t, iProbe, "HILoose"))
        //    continue;

        //  TLorentzVector mprobe;
        //  mprobe.SetPtEtaPhiM (t->muon_pt[iProbe], t->muon_eta[iProbe], t->muon_phi[iProbe], muon_mass);

        //  int iSign = (t->muon_charge[iProbe] == t->muon_charge[iTag] ? 0 : 1);

        //  TLorentzVector z = mtag+mprobe;
        //  if (z.M () < 76 || 106 < z.M ())
        //    continue;

        //  if (!isPbPb) {
        //    h_muonID_MSEffDen_pt[0][iSign]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h_muonID_MSEffDen_eta[0][iSign]->Fill (t->muon_eta[iProbe]);
        //    }
        //    if (t->muon_id_track_pt[iProbe] > 0) {
        //      h_muonID_MSEffNum_pt[0][iSign]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h_muonID_MSEffNum_eta[0][iSign]->Fill (t->muon_eta[iProbe]);
        //      }
        //    }
        //  }
        //  else if (isPbPb && !is2015data) {
        //    h_muonID_MSEffDen_pt[1][iSign]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h_muonID_MSEffDen_eta[1][iSign]->Fill (t->muon_eta[iProbe]);
        //      h_muonID_MSEffDen_fcal[0][iSign]->Fill (fcal_et);
        //    }
        //    if (t->muon_id_track_pt[iProbe] > 0) {
        //      h_muonID_MSEffNum_pt[1][iSign]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h_muonID_MSEffNum_eta[1][iSign]->Fill (t->muon_eta[iProbe]);
        //        h_muonID_MSEffNum_fcal[0][iSign]->Fill (fcal_et);
        //      }
        //    }
        //  }
        //  else if (isPbPb && is2015data) {
        //    h_muonID_MSEffDen_pt[2][iSign]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h_muonID_MSEffDen_eta[2][iSign]->Fill (t->muon_eta[iProbe]);
        //      h_muonID_MSEffDen_fcal[1][iSign]->Fill (fcal_et);
        //    }
        //    if (t->muon_id_track_pt[iProbe] > 0) {
        //      h_muonID_MSEffNum_pt[2][iSign]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h_muonID_MSEffNum_eta[2][iSign]->Fill (t->muon_eta[iProbe]);
        //        h_muonID_MSEffNum_fcal[1][iSign]->Fill (fcal_et);
        //      }
        //    }
        //  }
        //} // end loop over probe


        ////**** Fill muon ID efficiency hists. ****//
        //for (int iProbe = 0; iProbe < t->muon_n; iProbe++) {
        //  if (iProbe == iTag)
        //    continue;

        //  if (!isMC && !t->muon_matched[iProbe])
        //    continue;
        //  if (fabs (t->muon_eta[iProbe]) > 2.5)
        //    continue;

        //  if (t->muon_charge[iProbe] == t->muon_charge[iTag])
        //    continue;

        //  if (isMC) {
        //    bool hasTruthMatch = false;
        //    for (int iTM = 0; iTM < t->truth_muon_n; iTM++) {
        //      if (DeltaR (t->muon_eta[iProbe], t->truth_muon_eta[iTM], t->muon_phi[iProbe], t->truth_muon_phi[iTM]) < 0.3) {
        //        hasTruthMatch = true;
        //        break;
        //      }
        //    }
        //    if (!hasTruthMatch)
        //      continue;
        //  }

        //  TLorentzVector mprobe;
        //  mprobe.SetPtEtaPhiM (t->muon_pt[iProbe], t->muon_eta[iProbe], t->muon_phi[iProbe], muon_mass);

        //  TLorentzVector z = mtag+mprobe;
        //  if (z.M () < 76 || 106 < z.M ())
        //    continue;

        //  if (DeltaR (t->muon_eta[iTag], t->muon_eta[iProbe], t->muon_phi[iTag], t->muon_phi[iProbe]) < 0.3)
        //    continue; // avoids ambiguity of muons being too close together

        //  if (!isPbPb) {
        //    h_muonIDEffDen_pt[0]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h2_muonIDEffDen_eta_phi[0]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
        //    }
        //    if (t->muon_medium[iProbe]) {
        //      h_muonIDEffNum_pt[0]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h2_muonIDEffNum_eta_phi[0]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
        //      }
        //    }
        //  }
        //  else if (isPbPb && !is2015data) {
        //    h_muonIDEffDen_pt[1]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h2_muonIDEffDen_eta_phi[1]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
        //      h_muonIDEffDen_fcal[0]->Fill (fcal_et*1e-3);
        //    }
        //    if (t->muon_medium[iProbe]) {
        //      h_muonIDEffNum_pt[1]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h2_muonIDEffNum_eta_phi[1]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
        //        h_muonIDEffNum_fcal[0]->Fill (fcal_et*1e-3);
        //      }
        //    }
        //  }
        //  else if (isPbPb && is2015data) {
        //    h_muonIDEffDen_pt[2]->Fill (t->muon_pt[iProbe]);
        //    if (t->muon_pt[iProbe] > 20) {
        //      h2_muonIDEffDen_eta_phi[2]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
        //      h_muonIDEffDen_fcal[1]->Fill (fcal_et*1e-3);
        //    }
        //    if (t->muon_medium[iProbe]) {
        //      h_muonIDEffNum_pt[2]->Fill (t->muon_pt[iProbe]);
        //      if (t->muon_pt[iProbe] > 20) {
        //        h2_muonIDEffNum_eta_phi[2]->Fill (t->muon_eta[iProbe], InTwoPi (t->muon_phi[iProbe]));
        //        h_muonIDEffNum_fcal[1]->Fill (fcal_et);
        //      }
        //    }
        //  }
        //} // end loop over probe
      } // end loop over tag
    } // end if muon trigger fired (and passes_toroid)


    //**** Loop over tag electrons ****//
    if (isMC || electronTrigger->trigBool) {
      for (int iTag = 0; iTag < t->electron_n; iTag++) {
        if (!isMC && !t->electron_matched[iTag])
          continue;
        if ((isPbPb && !t->electron_lhmedium_hi[iTag]) || (!isPbPb && !t->electron_lhmedium[iTag]))
          continue;
        if (!InEMCal (t->electron_eta[iTag]))
          continue;
        if (t->electron_pt[iTag] < 30)
          continue;

        if (isMC) {
          bool hasTruthMatch = false;
          for (int iTE = 0; iTE < t->truth_electron_n; iTE++) {
            if (!InEMCal (t->truth_electron_eta[iTE]))
              continue;
            if (t->truth_electron_barcode[iTE] > 10000)
              continue;
            if (DeltaR (t->truth_electron_eta[iTE], t->electron_eta[iTag], t->truth_electron_phi[iTE], t->electron_phi[iTag]) < 0.3) {
              hasTruthMatch = true;
              break;
            }
          }
          if (!hasTruthMatch)
            continue;
        }

        TLorentzVector etag;
        etag.SetPtEtaPhiM (t->electron_pt[iTag], t->electron_eta[iTag], t->electron_phi[iTag], electron_mass);

        //**** Fill electron trigger efficiency hists. ****//
        for (int iProbe = 0; iProbe < t->electron_n; iProbe++) {
          if (iProbe == iTag)
            continue;

          if ((!isPbPb && !t->electron_lhloose[iProbe]) || (isPbPb && !t->electron_lhloose_hi[iProbe]))
            continue;
          if (!InEMCal (t->electron_eta[iProbe]))
            continue;

          if (t->electron_charge[iProbe] == t->electron_charge[iTag])
            continue;

          TLorentzVector eprobe;
          eprobe.SetPtEtaPhiM (t->electron_pt[iProbe], t->electron_eta[iProbe], t->electron_phi[iProbe], electron_mass);

          TLorentzVector z = etag+eprobe;
          if (z.M () < 76 || 106 < z.M ())
            continue;

          if (DeltaR (t->electron_eta[iTag], t->electron_eta[iProbe], t->electron_phi[iTag], t->electron_phi[iProbe]) < 0.3)
            continue; // avoids ambiguity of electrons being too close together

          if (!isPbPb) {
            h_electronTrigEffDen_pt[0]->Fill (t->electron_pt[iProbe]);
            h2_electronTrigEffDen_pt_eta[0]->Fill (t->electron_eta[iProbe], t->electron_pt[iProbe]);
            if (t->electron_pt[iProbe] > 20) {
              h_electronTrigEffDen_eta[0]->Fill (t->electron_eta[iProbe]);
            }
            if (isMC || t->electron_matched[iProbe]) {
              h_electronTrigEffNum_pt[0]->Fill (t->electron_pt[iProbe]);
              h2_electronTrigEffNum_pt_eta[0]->Fill (t->electron_eta[iProbe], t->electron_pt[iProbe]);
              if (t->electron_pt[iProbe] > 20) {
                h_electronTrigEffNum_eta[0]->Fill (t->electron_eta[iProbe]);
              }
            }
          }
          else if (isPbPb && !is2015data) {
            h_electronTrigEffDen_pt[1]->Fill (t->electron_pt[iProbe]);
            h2_electronTrigEffDen_pt_eta[1]->Fill (t->electron_eta[iProbe], t->electron_pt[iProbe]);
            if (t->electron_pt[iProbe] > 20) {
              h_electronTrigEffDen_eta[1]->Fill (t->electron_eta[iProbe]);
              //h_electronTrigEffDen_fcal[0]->Fill (fcal_et*1e-3);
            }
            if (isMC || t->electron_matched[iProbe]) {
              h_electronTrigEffNum_pt[1]->Fill (t->electron_pt[iProbe]);
              h2_electronTrigEffNum_pt_eta[1]->Fill (t->electron_eta[iProbe], t->electron_pt[iProbe]);
              if (t->electron_pt[iProbe] > 20) {
                h_electronTrigEffNum_eta[1]->Fill (t->electron_eta[iProbe]);
                //h_electronTrigEffNum_fcal[0]->Fill (fcal_et*1e-3);
              }
            }
          }
          else if (isPbPb && is2015data) {
            h_electronTrigEffDen_pt[2]->Fill (t->electron_pt[iProbe]);
            h2_electronTrigEffDen_pt_eta[2]->Fill (t->electron_eta[iProbe], t->electron_pt[iProbe]);
            if (t->electron_pt[iProbe] > 20) {
              h_electronTrigEffDen_eta[2]->Fill (t->electron_eta[iProbe]);
              //h_electronTrigEffDen_fcal[1]->Fill (fcal_et*1e-3);
            }
            if (isMC || t->electron_matched[iProbe]) {
              h_electronTrigEffNum_pt[2]->Fill (t->electron_pt[iProbe]);
              h2_electronTrigEffNum_pt_eta[2]->Fill (t->electron_eta[iProbe], t->electron_pt[iProbe]);
              if (t->electron_pt[iProbe] > 20) {
                h_electronTrigEffNum_eta[2]->Fill (t->electron_eta[iProbe]);
                //h_electronTrigEffNum_fcal[1]->Fill (fcal_et*1e-3);
              }
            }
          }
        } // end loop over probe


        ////**** Fill electron ID efficiency hists. ****//
        //for (int iProbe = 0; iProbe < t->Cluster_n; iProbe++) {

        //  if (!InEMCal (t->Cluster_etaBE->at (iProbe)))
        //    continue;
        //  if (DeltaR (t->electron_eta[iTag], t->Cluster_etaBE->at (iProbe), t->electron_phi[iTag], t->Cluster_phi->at (iProbe)) < 0.3)
        //    continue; // avoids tag electron cluster

        //  float plotPt = t->Cluster_pt->at (iProbe);
        //  float plotEta = t->Cluster_etaBE->at (iProbe);

        //  if (isMC) {
        //    bool hasTruthMatch = false;
        //    int iTE = 0;
        //    for (; iTE < t->truth_electron_n; iTE++) {
        //      if (!InEMCal (t->truth_electron_eta[iTE]))
        //        continue;
        //      if (t->truth_electron_barcode[iTE] > 10000)
        //        continue;
        //      if (DeltaR (t->truth_electron_eta[iTE], t->Cluster_etaBE->at (iProbe), t->truth_electron_phi[iTE], t->Cluster_phi->at (iProbe)) < 0.3) {
        //        hasTruthMatch = true;
        //        break;
        //      }
        //    }
        //    if (!hasTruthMatch)
        //      continue;
        //    plotPt = t->truth_electron_pt[iTE];
        //    plotEta = t->truth_electron_eta[iTE];
        //  }

        //  TLorentzVector eprobe;
        //  eprobe.SetPtEtaPhiM (t->Cluster_pt->at (iProbe), t->Cluster_etaBE->at (iProbe), t->Cluster_phi->at (iProbe), electron_mass);

        //  TLorentzVector z = etag+eprobe;
        //  if (z.M () < 76 || 106 < z.M ())
        //    continue;

        //  bool hasMatch = false;
        //  for (int iE = 0; iE < t->electron_n; iE++) {
        //    if (iE == iTag)
        //      continue;
        //    if (!InEMCal (t->electron_eta[iE]))
        //      continue;
        //    if ((isPbPb && !t->electron_lhloose_hi[iE]) || (!isPbPb && !t->electron_lhloose[iE]))
        //      continue;
        //    hasMatch = (DeltaR (t->electron_eta[iE], t->Cluster_etaBE->at (iProbe), t->electron_phi[iE], t->Cluster_phi->at (iProbe)) < 0.3);
        //  }

        //  if (!isPbPb) {
        //    h_electronIDEffDen_pt[0]->Fill (plotPt);
        //    h2_electronIDEffDen_pt_eta[0]->Fill (plotEta, plotPt);
        //    if (plotPt > 20) {
        //      h_electronIDEffDen_eta[0]->Fill (plotEta);
        //    }
        //    if (hasMatch) {
        //      h_electronIDEffNum_pt[0]->Fill (plotPt);
        //      h2_electronIDEffNum_pt_eta[0]->Fill (plotEta, plotPt);
        //      if (plotPt > 20) {
        //        h_electronIDEffNum_eta[0]->Fill (plotEta);
        //      }
        //    }
        //  }
        //  else if (isPbPb && !is2015data) {
        //    h_electronIDEffDen_pt[1]->Fill (plotPt);
        //    h2_electronIDEffDen_pt_eta[1]->Fill (plotEta, plotPt);
        //    if (plotPt > 20) {
        //      h_electronIDEffDen_eta[1]->Fill (plotEta);
        //    }
        //    if (hasMatch) {
        //      h_electronIDEffNum_pt[1]->Fill (plotPt);
        //      h2_electronIDEffNum_pt_eta[1]->Fill (plotEta, plotPt);
        //      if (plotPt > 20) {
        //        h_electronIDEffNum_eta[1]->Fill (plotEta);
        //      }
        //    }
        //  }
        //  else if (isPbPb && is2015data) {
        //    h_electronIDEffDen_pt[2]->Fill (plotPt);
        //    h2_electronIDEffDen_pt_eta[2]->Fill (plotEta, plotPt);
        //    if (plotPt > 20) {
        //      h_electronIDEffDen_eta[2]->Fill (plotEta);
        //    }
        //    if (hasMatch) {
        //      h_electronIDEffNum_pt[2]->Fill (plotPt);
        //      h2_electronIDEffNum_pt_eta[2]->Fill (plotEta, plotPt);
        //      if (plotPt > 20) {
        //        h_electronIDEffNum_eta[2]->Fill (plotEta);
        //      }
        //    }
        //  }
        //} // end loop over probe

      } // end loop over tag
    } // end if electron trigger fired


    //**** Get ID efficiencies from MC ****//
    if (isMC) {
      //**** Loop over truth muons ****//
      for (int iTM = 0; iTM < t->truth_muon_n; iTM++) {
        if (t->truth_muon_barcode[iTM] > 10000)
          continue;
        if (!InEMCal (t->truth_muon_eta[iTM]))
          continue;

        bool hasRecoMatch = false;
        for (int iM = 0; iM < t->muon_n; iM++) {
          if (!t->muon_medium[iM])
            continue;
          if (!InEMCal (t->muon_eta[iM]))
            continue;
          if (DeltaR (t->truth_muon_eta[iTM], t->muon_eta[iM], t->truth_muon_phi[iTM], t->muon_phi[iM]) < 0.3) {
            hasRecoMatch = true;
            break;
          }
        }

        if (!isPbPb) {
          h_muonIDEffDen_pt[0]->Fill (t->truth_muon_pt[iTM]);
          if (t->truth_muon_pt[iTM] > 20) {
            h2_muonIDEffDen_eta_phi[0]->Fill (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM]));
          }
          if (hasRecoMatch) {
            h_muonIDEffNum_pt[0]->Fill (t->truth_muon_pt[iTM]);
            if (t->truth_muon_pt[iTM] > 20) {
              h2_muonIDEffNum_eta_phi[0]->Fill (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM]));
            }
          }
        }
        else if (isPbPb && !is2015data) {
          h_muonIDEffDen_pt[1]->Fill (t->truth_muon_pt[iTM]);
          if (t->truth_muon_pt[iTM] > 20) {
            h2_muonIDEffDen_eta_phi[1]->Fill (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM]));
            //h_muonIDEffDen_fcal[0]->Fill (fcal_et*1e-3);
          }
          if (hasRecoMatch) {
            h_muonIDEffNum_pt[1]->Fill (t->truth_muon_pt[iTM]);
            if (t->truth_muon_pt[iTM] > 20) {
              h2_muonIDEffNum_eta_phi[1]->Fill (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM]));
              //h_muonIDEffNum_fcal[0]->Fill (fcal_et*1e-3);
            }
          }
        }
        else if (isPbPb && is2015data) {
          h_muonIDEffDen_pt[2]->Fill (t->truth_muon_pt[iTM]);
          if (t->truth_muon_pt[iTM] > 20) {
            h2_muonIDEffDen_eta_phi[2]->Fill (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM]));
            //h_muonIDEffDen_fcal[1]->Fill (fcal_et*1e-3);
          }
          if (hasRecoMatch) {
            h_muonIDEffNum_pt[2]->Fill (t->truth_muon_pt[iTM]);
            if (t->truth_muon_pt[iTM] > 20) {
              h2_muonIDEffNum_eta_phi[2]->Fill (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM]));
              //h_muonIDEffNum_fcal[1]->Fill (fcal_et);
            }
          }
        }
      } // end loop over truth electrons

      //**** Loop over truth electrons ****//
      for (int iTE = 0; iTE < t->truth_electron_n; iTE++) {
        if (t->truth_electron_barcode[iTE] > 10000)
          continue;
        if (!InEMCal (t->truth_electron_eta[iTE]))
          continue;

        bool hasRecoMatch = false;
        for (int iE = 0; iE < t->electron_n; iE++) {
          if ((isPbPb && !t->electron_lhloose_hi[iE]) || (!isPbPb && !t->electron_lhloose[iE]))
            continue;
          if (!InEMCal (t->electron_eta[iE]))
            continue;
          if (DeltaR (t->truth_electron_eta[iTE], t->electron_eta[iE], t->truth_electron_phi[iTE], t->electron_phi[iE]) < 0.3) {
            hasRecoMatch = true;
            break;
          }
        }

        if (!isPbPb) {
          h_electronIDEffDen_pt[0]->Fill (t->truth_electron_pt[iTE]);
          h2_electronIDEffDen_pt_eta[0]->Fill (t->truth_electron_eta[iTE], t->truth_electron_pt[iTE]);
          if (t->truth_electron_pt[iTE] > 20) {
            h_electronIDEffDen_eta[0]->Fill (t->truth_electron_eta[iTE]);
          }
          if (hasRecoMatch) {
            h_electronIDEffNum_pt[0]->Fill (t->truth_electron_pt[iTE]);
            h2_electronIDEffNum_pt_eta[0]->Fill (t->truth_electron_eta[iTE], t->truth_electron_pt[iTE]);
            if (t->truth_electron_pt[iTE] > 20) {
              h_electronIDEffNum_eta[0]->Fill (t->truth_electron_eta[iTE]);
            }
          }
        }
        else if (isPbPb && !is2015data) {
          h_electronIDEffDen_pt[1]->Fill (t->truth_electron_pt[iTE]);
          h2_electronIDEffDen_pt_eta[1]->Fill (t->truth_electron_eta[iTE], t->truth_electron_pt[iTE]);
          if (t->truth_electron_pt[iTE] > 20) {
            h_electronIDEffDen_eta[1]->Fill (t->truth_electron_eta[iTE]);
            //h_electronIDEffDen_fcal[0]->Fill (fcal_et*1e-3);
          }
          if (hasRecoMatch) {
            h_electronIDEffNum_pt[1]->Fill (t->truth_electron_pt[iTE]);
            h2_electronIDEffNum_pt_eta[1]->Fill (t->truth_electron_eta[iTE], t->truth_electron_pt[iTE]);
            if (t->truth_electron_pt[iTE] > 20) {
              h_electronIDEffNum_eta[1]->Fill (t->truth_electron_eta[iTE]);
              //h_electronIDEffNum_fcal[0]->Fill (fcal_et*1e-3);
            }
          }
        }
        else if (isPbPb && is2015data) {
          h_electronIDEffDen_pt[2]->Fill (t->truth_electron_pt[iTE]);
          h2_electronIDEffDen_pt_eta[2]->Fill (t->truth_electron_eta[iTE], t->truth_electron_pt[iTE]);
          if (t->truth_electron_pt[iTE] > 20) {
            h_electronIDEffDen_eta[2]->Fill (t->truth_electron_eta[iTE]);
            //h_electronIDEffDen_fcal[1]->Fill (fcal_et*1e-3);
          }
          if (hasRecoMatch) {
            h_electronIDEffNum_pt[2]->Fill (t->truth_electron_pt[iTE]);
            h2_electronIDEffNum_pt_eta[2]->Fill (t->truth_electron_eta[iTE], t->truth_electron_pt[iTE]);
            if (t->truth_electron_pt[iTE] > 20) {
              h_electronIDEffNum_eta[2]->Fill (t->truth_electron_eta[iTE]);
              //h_electronIDEffNum_fcal[1]->Fill (fcal_et*1e-3);
            }
          }
        }
      } // end loop over truth electrons
    } // end if MC


  } // end loop over events



  //////////////////////////////////////////////////////////////
  //**** Calculate single lepton trigger efficiency plots ****//
  //////////////////////////////////////////////////////////////

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TagAndProbe.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);
  }

  TEfficiency* t_muonTrigEff_eta_phi_pp = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[0]), *(h2_muonTrigEffDen_eta_phi[0]));
  TEfficiency* t_muonTrigEff_eta_phi_PbPb18 = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[1]), *(h2_muonTrigEffDen_eta_phi[1]));
  TEfficiency* t_muonTrigEff_eta_phi_PbPb15 = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[2]), *(h2_muonTrigEffDen_eta_phi[2]));
  TEfficiency* t_muonIDEff_eta_phi_pp = new TEfficiency (*(h2_muonIDEffNum_eta_phi[0]), *(h2_muonIDEffDen_eta_phi[0]));
  TEfficiency* t_muonIDEff_eta_phi_PbPb18 = new TEfficiency (*(h2_muonIDEffNum_eta_phi[1]), *(h2_muonIDEffDen_eta_phi[1]));
  TEfficiency* t_muonIDEff_eta_phi_PbPb15 = new TEfficiency (*(h2_muonIDEffNum_eta_phi[2]), *(h2_muonIDEffDen_eta_phi[2]));

  //TEfficiency* t_muonMed_IDEff_eta_pp[2];
  //TEfficiency* t_muonMed_IDEff_eta_PbPb18[2];
  //TEfficiency* t_muonMed_IDEff_eta_PbPb15[2];
  //TEfficiency* t_muonID_MSEff_eta_pp[2];
  //TEfficiency* t_muonID_MSEff_eta_PbPb18[2];
  //TEfficiency* t_muonID_MSEff_eta_PbPb15[2];
  //for (int iSign : {0, 1}) {
  //  t_muonMed_IDEff_eta_pp[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[0][iSign]), *(h_muonMed_IDEffDen_eta[0][iSign]));
  //  t_muonMed_IDEff_eta_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[1][iSign]), *(h_muonMed_IDEffDen_eta[1][iSign]));
  //  t_muonMed_IDEff_eta_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[2][iSign]), *(h_muonMed_IDEffDen_eta[2][iSign]));
  //  t_muonID_MSEff_eta_pp[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[0][iSign]), *(h_muonID_MSEffDen_eta[0][iSign]));
  //  t_muonID_MSEff_eta_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[1][iSign]), *(h_muonID_MSEffDen_eta[1][iSign]));
  //  t_muonID_MSEff_eta_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[2][iSign]), *(h_muonID_MSEffDen_eta[2][iSign]));
  //}

  TEfficiency* t_electronTrigEff_pt_eta_pp = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[0]), *(h2_electronTrigEffDen_pt_eta[0]));
  TEfficiency* t_electronTrigEff_pt_eta_PbPb18 = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[1]), *(h2_electronTrigEffDen_pt_eta[1]));
  TEfficiency* t_electronTrigEff_pt_eta_PbPb15 = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[2]), *(h2_electronTrigEffDen_pt_eta[2]));

  TEfficiency* t_electronIDEff_pt_eta_pp = new TEfficiency (*(h2_electronIDEffNum_pt_eta[0]), *(h2_electronIDEffDen_pt_eta[0]));
  TEfficiency* t_electronIDEff_pt_eta_PbPb18 = new TEfficiency (*(h2_electronIDEffNum_pt_eta[1]), *(h2_electronIDEffDen_pt_eta[1]));
  TEfficiency* t_electronIDEff_pt_eta_PbPb15 = new TEfficiency (*(h2_electronIDEffNum_pt_eta[2]), *(h2_electronIDEffDen_pt_eta[2]));



  /////////////////////////////////////////////////////////////////////////////////
  //**** Calculate residual centrality dependent efficiencies ****//
  /////////////////////////////////////////////////////////////////////////////////

  // Loop over events
  if (isPbPb) {
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << "Info: In TagAndProbe.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      tree->GetEntry (iEvt);

      //**** Event selection ****//
      if (collisionSystem == PbPb18 && t->BlayerDesyn)
        continue;
      if (isPbPb && t->isOOTPU)
        continue;

      bool hasPrimaryVert = false;
      //bool hasPileupVert = false;
      for (int iVert = 0; iVert < t->nvert; iVert++) {
        if (t->vert_type[iVert] == 1) {
          if (hasPrimaryVert) {
            hasPrimaryVert = false;
            break;
          }
          hasPrimaryVert = (t->vert_type[iVert] == 1);
          vz = t->vert_z[iVert];
        }
        //else if (t->vert_type[iVert] == 3) {
        //  hasPileupVert = true;
        //}
      }
      if (!hasPrimaryVert)// || hasPileupVert)
        continue;

      if (isPbPb) {
        fcal_et = t->fcalA_et + t->fcalC_et;

        const float zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
        const float nNeutrons = (zdcEnergy) / (2.51);
        const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3);
        if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
          continue; // Zdc-based in-time pile-up cut
      }

      //**** Loop over tag electrons ****//
      if (isMC || electronTrigger->trigBool) {
        for (int iTag = 0; iTag < t->electron_n; iTag++) {
          if (!isMC && !t->electron_matched[iTag])
            continue;
          if (!t->electron_lhmedium_hi[iTag])
            continue;
          if (!InEMCal (t->electron_eta[iTag]))
            continue;
          if (t->electron_pt[iTag] < 30)
            continue;

          if (isMC) {
            bool hasTruthMatch = false;
            for (int iTE = 0; iTE < t->truth_electron_n; iTE++) {
              if (!InEMCal (t->truth_electron_eta[iTE]))
                continue;
              if (t->truth_electron_barcode[iTE] > 10000)
                continue;
              if (DeltaR (t->truth_electron_eta[iTE], t->electron_eta[iTag], t->truth_electron_phi[iTE], t->electron_phi[iTag]) < 0.3) {
                hasTruthMatch = true;
                break;
              }
            }
            if (!hasTruthMatch)
              continue;
          }

          TLorentzVector etag;
          etag.SetPtEtaPhiM (t->electron_pt[iTag], t->electron_eta[iTag], t->electron_phi[iTag], electron_mass);


          //**** Fill electron trigger efficiency hists. ****//
          for (int iProbe = 0; iProbe < t->electron_n; iProbe++) {
            if (iProbe == iTag)
              continue;

            if (!t->electron_lhloose_hi[iProbe])
              continue;
            if (!InEMCal (t->electron_eta[iProbe]))
              continue;

            if (t->electron_pt[iProbe] < 20)
              continue;

            if (t->electron_charge[iProbe] == t->electron_charge[iTag])
              continue;

            TLorentzVector eprobe;
            eprobe.SetPtEtaPhiM (t->electron_pt[iProbe], t->electron_eta[iProbe], t->electron_phi[iProbe], electron_mass);

            TLorentzVector z = etag+eprobe;
            if (z.M () < 76 || 106 < z.M ())
              continue;

            if (DeltaR (t->electron_eta[iTag], t->electron_eta[iProbe], t->electron_phi[iTag], t->electron_phi[iProbe]) < 0.3)
              continue; // avoids ambiguity of electrons being too close together

            const float effprobe = (is2015data ? t_electronTrigEff_pt_eta_PbPb15->GetEfficiency (t_electronTrigEff_pt_eta_PbPb15->FindFixBin (t->electron_eta[iProbe], t->electron_pt[iProbe])) : t_electronTrigEff_pt_eta_PbPb18->GetEfficiency (t_electronTrigEff_pt_eta_PbPb18->FindFixBin (t->electron_eta[iProbe], t->electron_pt[iProbe])));
            if (effprobe == 0)
              continue;

            h_electronTrigEffDen_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
            if (isMC || t->electron_matched[iProbe])
              h_electronTrigEffNum_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
          } // end loop over probe


          ////**** Fill electron ID efficiency hists. ****//
          //for (int iProbe = 0; iProbe < t->Cluster_n; iProbe++) {

          //  if (!InEMCal (t->Cluster_etaBE->at (iProbe)))
          //    continue; // pointing cuts
          //  if (DeltaR (t->electron_eta[iTag], t->Cluster_etaBE->at (iProbe), t->electron_phi[iTag], t->Cluster_phi->at (iProbe)) < 0.3)
          //    continue; // avoids tag electron cluster

          //  float plotPt = t->Cluster_pt->at (iProbe);
          //  float plotEta = t->Cluster_etaBE->at (iProbe);

          //  if (isMC) {
          //    bool hasTruthMatch = false;
          //    int iTE = 0;
          //    for (; iTE < t->truth_electron_n; iTE++) {
          //      if (!InEMCal (t->truth_electron_eta[iTE]))
          //        continue;
          //      if (t->truth_electron_barcode[iTE] > 10000)
          //        continue;
          //      if (t->truth_electron_pt[iTE] < 20)
          //        continue;
          //      if (DeltaR (t->truth_electron_eta[iTE], t->Cluster_etaBE->at (iProbe), t->truth_electron_phi[iTE], t->Cluster_phi->at (iProbe)) < 0.3) {
          //        hasTruthMatch = true;
          //        break;
          //      }
          //    }
          //    if (!hasTruthMatch)
          //      continue;
          //    plotPt = t->truth_electron_pt[iTE];
          //    plotEta = t->truth_electron_eta[iTE];
          //  }

          //  TLorentzVector eprobe;
          //  eprobe.SetPtEtaPhiM (t->Cluster_pt->at (iProbe), t->Cluster_etaBE->at (iProbe), t->Cluster_phi->at (iProbe), electron_mass);

          //  TLorentzVector z = etag+eprobe;
          //  if (z.M () < 76 || 106 < z.M ())
          //    continue; // reduce background by cutting on Z mass

          //  bool hasMatch = false;
          //  for (int iElectron = 0; iElectron < t->electron_n; iElectron++) {
          //    if (iElectron == iTag)
          //      continue; // don't consider tag electron
          //    //if (!InEMCal (t->electron_eta[iElectron]))
          //    //  continue; // pointing cuts
          //    if (!t->electron_lhloose_hi[iElectron])
          //      continue; // ID cuts
          //    hasMatch = (DeltaR (t->electron_eta[iElectron], t->Cluster_etaBE->at (iProbe), t->electron_phi[iElectron], t->Cluster_phi->at (iProbe)) < 0.3); // cluster-electron matching
          //  }

          //  const float effprobe = (is2015data ? t_electronIDEff_pt_eta_PbPb15->GetEfficiency (t_electronIDEff_pt_eta_PbPb15->FindFixBin (plotEta, plotPt)) : t_electronIDEff_pt_eta_PbPb18->GetEfficiency (t_electronIDEff_pt_eta_PbPb18->FindFixBin (plotEta, plotPt)));
          //  if (effprobe == 0)
          //    continue;

          //  h_electronIDEffDen_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
          //  if (hasMatch)
          //    h_electronIDEffNum_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
          //} // end loop over probe

        } // end loop over tag
      } // end if electron trigger fired


      //**** Loop over truth muons ****//
      if (isMC) {
        for (int iTM = 0; iTM < t->truth_muon_n; iTM++) {
          if (t->truth_muon_barcode[iTM] > 10000)
            continue;
          if (!InEMCal (t->truth_muon_eta[iTM]))
            continue;

          if (t->truth_muon_pt[iTM] < 20)
            continue;

          bool hasRecoMatch = false;
          for (int iM = 0; iM < t->muon_n; iM++) {
            if (!t->muon_medium[iM])
              continue;
            if (!InEMCal (t->muon_eta[iM]))
              continue;
            if (DeltaR (t->truth_muon_eta[iTM], t->muon_eta[iM], t->truth_muon_phi[iTM], t->muon_phi[iM]) < 0.3) {
              hasRecoMatch = true;
              break;
            }
          }

          TEfficiency* h = (is2015data ? t_muonIDEff_eta_phi_PbPb15 : t_muonIDEff_eta_phi_PbPb18);
          const float effprobe = h->GetEfficiency (h->FindFixBin (t->truth_muon_eta[iTM], InTwoPi (t->truth_muon_phi[iTM])));
          if (effprobe == 0)
            continue;

          h_muonIDEffDen_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
          if (hasRecoMatch)
            h_muonIDEffNum_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
        } // end loop over truth muons
      } // end if isMC


      //**** Loop over truth electrons ****//
      if (isMC) {
        for (int iTE = 0; iTE < t->truth_electron_n; iTE++) {
          if (t->truth_electron_barcode[iTE] > 10000)
            continue;
          if (!InEMCal (t->truth_electron_eta[iTE]))
            continue;

          if (t->truth_electron_pt[iTE] < 20)
            continue;

          bool hasRecoMatch = false;
          for (int iE = 0; iE < t->electron_n; iE++) {
            if ((isPbPb && !t->electron_lhloose_hi[iE]) || (!isPbPb && !t->electron_lhloose[iE]))
              continue;
            if (!InEMCal (t->electron_eta[iE]))
              continue;
            if (DeltaR (t->truth_electron_eta[iTE], t->electron_eta[iE], t->truth_electron_phi[iTE], t->electron_phi[iE]) < 0.3) {
              hasRecoMatch = true;
              break;
            }
          }

          TEfficiency* h = (is2015data ? t_electronIDEff_pt_eta_PbPb15 : t_electronIDEff_pt_eta_PbPb18);
          const float effprobe = h->GetEfficiency (h->FindFixBin (t->truth_electron_eta[iTE], fmin (t->truth_electron_pt[iTE], h->GetPassedHistogram ()->GetYaxis ()->GetBinCenter (h->GetPassedHistogram ()->GetNbinsY ()))));
          if (effprobe == 0)
            continue;

          h_electronIDEffDen_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
          if (hasRecoMatch)
            h_electronIDEffNum_fcal[is2015data]->Fill (fcal_et*1e-3, 1./effprobe);
        } // end loop over truth electrons
      } // end if isMC
    } // end loop over events
  } // end if isPbPb


  TH1D* h_muonIDEff_fcal[2];
  TH1D* h_electronTrigEff_fcal[2];
  TH1D* h_electronIDEff_fcal[2];
  h_muonIDEff_fcal[0] = (TH1D*) h_muonIDEffNum_fcal[0]->Clone ("h_muonIDEff_fcal_PbPb18");
  h_muonIDEff_fcal[1] = (TH1D*) h_muonIDEffNum_fcal[1]->Clone ("h_muonIDEff_fcal_PbPb15");
  h_muonIDEff_fcal[0]->Divide (h_muonIDEffDen_fcal[0]);
  h_muonIDEff_fcal[1]->Divide (h_muonIDEffDen_fcal[1]);
  h_electronTrigEff_fcal[0] = (TH1D*) h_electronTrigEffNum_fcal[0]->Clone ("h_electronTrigEff_fcal_PbPb18");
  h_electronTrigEff_fcal[1] = (TH1D*) h_electronTrigEffNum_fcal[1]->Clone ("h_electronTrigEff_fcal_PbPb15");
  h_electronTrigEff_fcal[0]->Divide (h_electronTrigEffDen_fcal[0]);
  h_electronTrigEff_fcal[1]->Divide (h_electronTrigEffDen_fcal[1]);
  h_electronIDEff_fcal[0] = (TH1D*) h_electronIDEffNum_fcal[0]->Clone ("h_electronIDEff_fcal_PbPb18");
  h_electronIDEff_fcal[1] = (TH1D*) h_electronIDEffNum_fcal[1]->Clone ("h_electronIDEff_fcal_PbPb15");
  h_electronIDEff_fcal[0]->Divide (h_electronIDEffDen_fcal[0]);
  h_electronIDEff_fcal[1]->Divide (h_electronIDEffDen_fcal[1]);




  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //**** Calculate Z boson trigger and reconstruction efficiencies via reweighted distributions ****//
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  // Loop over events
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TagAndProbe.cxx: Event loop " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    tree->GetEntry (iEvt);

    //**** Event selection ****//
    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;

    bool hasPrimaryVert = false;
    //bool hasPileupVert = false;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      if (t->vert_type[iVert] == 1) {
        if (hasPrimaryVert) {
          hasPrimaryVert = false;
          break;
        }
        hasPrimaryVert = (t->vert_type[iVert] == 1);
        vz = t->vert_z[iVert];
      }
      //else if (t->vert_type[iVert] == 3) {
      //  hasPileupVert = true;
      //}
    }
    if (!hasPrimaryVert)// || hasPileupVert)
      continue;

    if (isPbPb) {
      fcal_et = t->fcalA_et + t->fcalC_et;

      const float zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3);
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

    if ((isMC || muonTrigger->trigBool) && t->passes_toroid) {
      for (int iM1 = 0; iM1 < t->muon_n; iM1++) {
        if (!t->muon_medium[iM1])
          continue;

        const float reco_eff_cent = (isPbPb ? h_muonIDEff_fcal[is2015data]->GetBinContent (h_muonIDEff_fcal[is2015data]->FindFixBin (fcal_et*1e-3)) : 1);

        TEfficiency* h = (isPbPb ? (is2015data ? t_muonTrigEff_eta_phi_PbPb15 : t_muonTrigEff_eta_phi_PbPb18) : t_muonTrigEff_eta_phi_pp);
        const float trig_eff_m1 = h->GetEfficiency (h->FindFixBin (t->muon_eta[iM1], InTwoPi (t->muon_phi[iM1])));
        h = (isPbPb ? (is2015data ? t_muonIDEff_eta_phi_PbPb15 : t_muonIDEff_eta_phi_PbPb18) : t_muonIDEff_eta_phi_pp);
        const float reco_eff_m1 = h->GetEfficiency (h->FindFixBin (t->muon_eta[iM1], InTwoPi (t->muon_phi[iM1]))) * reco_eff_cent;

        TLorentzVector m1_tlv;
        m1_tlv.SetPtEtaPhiM (t->muon_pt[iM1], t->muon_eta[iM1], t->muon_phi[iM1], muon_mass);

        for (int iM2 = 0; iM2 < iM1; iM2++) {

          if (!t->muon_medium[iM2])
            continue;

          h = (isPbPb ? (is2015data ? t_muonTrigEff_eta_phi_PbPb15 : t_muonTrigEff_eta_phi_PbPb18) : t_muonTrigEff_eta_phi_pp);
          const float trig_eff_m2 = h->GetEfficiency (h->FindFixBin (t->muon_eta[iM2], InTwoPi (t->muon_phi[iM2])));
          h = (isPbPb ? (is2015data ? t_muonIDEff_eta_phi_PbPb15 : t_muonIDEff_eta_phi_PbPb18) : t_muonIDEff_eta_phi_pp);
          const float reco_eff_m2 = h->GetEfficiency (h->FindFixBin (t->muon_eta[iM2], InTwoPi (t->muon_phi[iM2]))) * reco_eff_cent;

          TLorentzVector m2_tlv;
          m2_tlv.SetPtEtaPhiM (t->muon_pt[iM2], t->muon_eta[iM2], t->muon_phi[iM2], muon_mass);

          if (t->muon_charge[iM1] == t->muon_charge[iM2])
            continue;

          TLorentzVector z_tlv = m1_tlv+m2_tlv;
          if (z_tlv.M () < 76 || 106 < z_tlv.M ())
            continue;

          const float trig_eff_z = 1. - (1.-trig_eff_m1)*(1.-trig_eff_m2);
          const float reco_eff_z = reco_eff_m1 * reco_eff_m2;
         
          if (!isPbPb) {
            h2_zmumuTrigEffNum_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), trig_eff_z);
            h2_zmumuTrigEffDen_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
            h2_zmumuIDEffNum_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), reco_eff_z);
            h2_zmumuIDEffDen_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
          } 
          else if (isPbPb && !is2015data) {
            h2_zmumuTrigEffNum_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), trig_eff_z);
            h2_zmumuTrigEffDen_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
            h2_zmumuIDEffNum_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), reco_eff_z);
            h2_zmumuIDEffDen_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
          } 
          else if (isPbPb && is2015data) {
            h2_zmumuTrigEffNum_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), trig_eff_z);
            h2_zmumuTrigEffDen_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
            h2_zmumuIDEffNum_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), reco_eff_z);
            h2_zmumuIDEffDen_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
          } 
          
        } // end loop over iM2
      } // end loop over iM1
    } // end if muon trigger fired (and passes_toroid)

    if (isMC || electronTrigger->trigBool) {
      for (int iE1 = 0; iE1 < t->electron_n; iE1++) {
        if ((isPbPb && !t->electron_lhloose_hi[iE1]) || (!isPbPb && !t->electron_lhloose[iE1]))
          continue;
        if (!InEMCal (t->electron_eta[iE1]))
          continue;

        if (t->electron_pt[iE1] < 20)
          continue;

        //int iC1 = -1;
        //while (iC1 == -1 || (iC1 < t->Cluster_n && DeltaR (t->electron_eta[iE1], t->Cluster_etaBE->at (iC1), t->electron_phi[iE1], t->Cluster_phi->at (iC1)) > 0.3)) {
        //  iC1++;
        //}
        //if (iC1 == -1 || iC1 == t->Cluster_n) {
        //  cout << "Warning: Cannot find best cluster for this electron! Continuing..." << endl;
        //  continue;
        //}

        const float trig_eff_cent = (isPbPb ? h_electronTrigEff_fcal[is2015data]->GetBinContent (h_electronTrigEff_fcal[is2015data]->FindFixBin (fcal_et*1e-3)) : 1);
        const float reco_eff_cent = (isPbPb ? h_electronIDEff_fcal[is2015data]->GetBinContent (h_electronIDEff_fcal[is2015data]->FindFixBin (fcal_et*1e-3)) : 1);

        TEfficiency* h = (isPbPb ? (is2015data ? t_electronTrigEff_pt_eta_PbPb15 : t_electronTrigEff_pt_eta_PbPb18) : t_electronTrigEff_pt_eta_pp);
        const float trig_eff_e1 = h->GetEfficiency (h->FindFixBin (t->electron_eta[iE1], fmin (t->electron_pt[iE1], h->GetPassedHistogram ()->GetYaxis ()->GetBinCenter (h->GetPassedHistogram ()->GetNbinsY ())))) * trig_eff_cent;
        h = (isPbPb ? (is2015data ? t_electronIDEff_pt_eta_PbPb15 : t_electronIDEff_pt_eta_PbPb18) : t_electronIDEff_pt_eta_pp);
        const float reco_eff_e1 = h->GetEfficiency (h->FindFixBin (t->electron_eta[iE1], fmin (t->electron_pt[iE1], h->GetPassedHistogram ()->GetYaxis ()->GetBinCenter (h->GetPassedHistogram ()->GetNbinsY ())))) * reco_eff_cent;

        TLorentzVector e1;
        e1.SetPtEtaPhiM (t->electron_pt[iE1], t->electron_eta[iE1], t->electron_phi[iE1], electron_mass);

        for (int iE2 = 0; iE2 < iE1; iE2++) {
          if ((isPbPb && !t->electron_lhloose_hi[iE2]) || (!isPbPb && !t->electron_lhloose[iE2]))
            continue;
          if (!InEMCal (t->electron_eta[iE2]))
            continue;

          if (t->electron_pt[iE2] < 20)
            continue;

          //int iC2 = -1;
          //while (iC2 == -1 || (iC2 < t->Cluster_n && DeltaR (t->electron_eta[iE1], t->Cluster_etaBE->at (iC2), t->electron_phi[iE1], t->Cluster_phi->at (iC2)) > 0.3)) {
          //  iC2++;
          //}
          //if (iC2 == -1 || iC2 == t->Cluster_n) {
          //  cout << "Warning: Cannot find best cluster for this electron! Continuing..." << endl;
          //  continue;
          //}

          if (t->electron_charge[iE2] == t->electron_charge[iE1])
            continue;

          h = (isPbPb ? (is2015data ? t_electronTrigEff_pt_eta_PbPb15 : t_electronTrigEff_pt_eta_PbPb18) : t_electronTrigEff_pt_eta_pp);
          const float trig_eff_e2 = h->GetEfficiency (h->FindFixBin (t->electron_eta[iE2], fmin (t->electron_pt[iE2], h->GetPassedHistogram ()->GetYaxis ()->GetBinCenter (h->GetPassedHistogram ()->GetNbinsY ())))) * trig_eff_cent;
          h = (isPbPb ? (is2015data ? t_electronIDEff_pt_eta_PbPb15 : t_electronIDEff_pt_eta_PbPb18) : t_electronIDEff_pt_eta_pp);
          const float reco_eff_e2 = h->GetEfficiency (h->FindFixBin (t->electron_eta[iE2], fmin (t->electron_pt[iE2], h->GetPassedHistogram ()->GetYaxis ()->GetBinCenter (h->GetPassedHistogram ()->GetNbinsY ())))) * reco_eff_cent;

          TLorentzVector e2;
          e2.SetPtEtaPhiM (t->electron_pt[iE2], t->electron_eta[iE2], t->electron_phi[iE2], electron_mass);

          TLorentzVector z_tlv = e1+e2;
          if (z_tlv.M () < 76 || 106 < z_tlv.M ())
            continue;

          const float trig_eff_z = (1. - (1.-trig_eff_e1)*(1.-trig_eff_e2));
          const float reco_eff_z = reco_eff_e1 * reco_eff_e2;

          if (!isPbPb) {
            h2_zeeTrigEffNum_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), trig_eff_z);
            h2_zeeTrigEffDen_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
            h2_zeeIDEffNum_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), reco_eff_z);
            h2_zeeIDEffDen_pt_y[0]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
          } 
          else if (isPbPb && !is2015data) {
            h2_zeeTrigEffNum_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), trig_eff_z);
            h2_zeeTrigEffDen_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
            h2_zeeIDEffNum_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), reco_eff_z);
            h2_zeeIDEffDen_pt_y[1]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
          } 
          else if (isPbPb && is2015data) {
            h2_zeeTrigEffNum_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), trig_eff_z);
            h2_zeeTrigEffDen_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
            h2_zeeIDEffNum_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt (), reco_eff_z);
            h2_zeeIDEffDen_pt_y[2]->Fill (z_tlv.Rapidity (), z_tlv.Pt ());
          } 

        } // end loop over iE2
      } // end loop over iE1
    } // end if electron trigger fired

  }

  outFile->cd ();

  if (!isMC) {
    h_electronTrigEffNum_pt[0]->Write ();
    h_electronTrigEffDen_pt[0]->Write ();
    h_electronTrigEffNum_pt[1]->Write ();
    h_electronTrigEffDen_pt[1]->Write ();
    h_electronTrigEffNum_pt[2]->Write ();
    h_electronTrigEffDen_pt[2]->Write ();
    h_electronTrigEffNum_eta[0]->Write ();
    h_electronTrigEffDen_eta[0]->Write ();
    h_electronTrigEffNum_eta[1]->Write ();
    h_electronTrigEffDen_eta[1]->Write ();
    h_electronTrigEffNum_eta[2]->Write ();
    h_electronTrigEffDen_eta[2]->Write ();
    h2_electronTrigEffNum_pt_eta[0]->Write ();
    h2_electronTrigEffDen_pt_eta[0]->Write ();
    h2_electronTrigEffNum_pt_eta[1]->Write ();
    h2_electronTrigEffDen_pt_eta[1]->Write ();
    h2_electronTrigEffNum_pt_eta[2]->Write ();
    h2_electronTrigEffDen_pt_eta[2]->Write ();
    h_electronTrigEffNum_fcal[0]->Write ();
    h_electronTrigEffDen_fcal[0]->Write ();
    h_electronTrigEffNum_fcal[1]->Write ();
    h_electronTrigEffDen_fcal[1]->Write ();

    t_electronTrigEff_pt_eta_pp->Write ("t_electronTrigEff_pt_eta_pp");
    t_electronTrigEff_pt_eta_PbPb18->Write ("t_electronTrigEff_pt_eta_PbPb18");
    t_electronTrigEff_pt_eta_PbPb15->Write ("t_electronTrigEff_pt_eta_PbPb15");
    h_electronTrigEff_fcal[0]->Write ();
    h_electronTrigEff_fcal[1]->Write ();

    h2_zeeTrigEffNum_pt_y[0]->Write ();
    h2_zeeTrigEffDen_pt_y[0]->Write ();
    h2_zeeTrigEffNum_pt_y[1]->Write ();
    h2_zeeTrigEffDen_pt_y[1]->Write ();
    h2_zeeTrigEffNum_pt_y[2]->Write ();
    h2_zeeTrigEffDen_pt_y[2]->Write ();

    h_muonTrigEffNum_pt[0]->Write ();
    h_muonTrigEffDen_pt[0]->Write ();
    h_muonTrigEffNum_pt[1]->Write ();
    h_muonTrigEffDen_pt[1]->Write ();
    h_muonTrigEffNum_pt[2]->Write ();
    h_muonTrigEffDen_pt[2]->Write ();
    h2_muonTrigEffNum_eta_phi[0]->Write ();
    h2_muonTrigEffDen_eta_phi[0]->Write ();
    h2_muonTrigEffNum_eta_phi[1]->Write ();
    h2_muonTrigEffDen_eta_phi[1]->Write ();
    h2_muonTrigEffNum_eta_phi[2]->Write ();
    h2_muonTrigEffDen_eta_phi[2]->Write ();
    h_muonTrigEffNum_fcal[0]->Write ();
    h_muonTrigEffDen_fcal[0]->Write ();
    h_muonTrigEffNum_fcal[1]->Write ();
    h_muonTrigEffDen_fcal[1]->Write ();

    t_muonTrigEff_eta_phi_pp->Write ("t_muonTrigEff_eta_phi_pp");
    t_muonTrigEff_eta_phi_PbPb18->Write ("t_muonTrigEff_eta_phi_PbPb18");
    t_muonTrigEff_eta_phi_PbPb15->Write ("t_muonTrigEff_eta_phi_PbPb15");

    h2_zmumuTrigEffNum_pt_y[0]->Write ();
    h2_zmumuTrigEffDen_pt_y[0]->Write ();
    h2_zmumuTrigEffNum_pt_y[1]->Write ();
    h2_zmumuTrigEffDen_pt_y[1]->Write ();
    h2_zmumuTrigEffNum_pt_y[2]->Write ();
    h2_zmumuTrigEffDen_pt_y[2]->Write ();

  }

  if (isMC) {

    h_muonIDEffNum_pt[0]->Write ();
    h_muonIDEffDen_pt[0]->Write ();
    h_muonIDEffNum_pt[1]->Write ();
    h_muonIDEffDen_pt[1]->Write ();
    h_muonIDEffNum_pt[2]->Write ();
    h_muonIDEffDen_pt[2]->Write ();
    h2_muonIDEffNum_eta_phi[0]->Write ();
    h2_muonIDEffDen_eta_phi[0]->Write ();
    h2_muonIDEffNum_eta_phi[1]->Write ();
    h2_muonIDEffDen_eta_phi[1]->Write ();
    h2_muonIDEffNum_eta_phi[2]->Write ();
    h2_muonIDEffDen_eta_phi[2]->Write ();
    h_muonIDEffNum_fcal[0]->Write ();
    h_muonIDEffDen_fcal[0]->Write ();
    h_muonIDEffNum_fcal[1]->Write ();
    h_muonIDEffDen_fcal[1]->Write ();

    //for (int iSign : {0, 1}) {
    //  h_muonMed_IDEffNum_pt[0][iSign]->Write ();
    //  h_muonMed_IDEffDen_pt[0][iSign]->Write ();
    //  h_muonMed_IDEffNum_pt[1][iSign]->Write ();
    //  h_muonMed_IDEffDen_pt[1][iSign]->Write ();
    //  h_muonMed_IDEffNum_pt[2][iSign]->Write ();
    //  h_muonMed_IDEffDen_pt[2][iSign]->Write ();
    //  h_muonMed_IDEffNum_eta[0][iSign]->Write ();
    //  h_muonMed_IDEffDen_eta[0][iSign]->Write ();
    //  h_muonMed_IDEffNum_eta[1][iSign]->Write ();
    //  h_muonMed_IDEffDen_eta[1][iSign]->Write ();
    //  h_muonMed_IDEffNum_eta[2][iSign]->Write ();
    //  h_muonMed_IDEffDen_eta[2][iSign]->Write ();
    //  h_muonMed_IDEffNum_fcal[0][iSign]->Write ();
    //  h_muonMed_IDEffDen_fcal[0][iSign]->Write ();
    //  h_muonMed_IDEffNum_fcal[1][iSign]->Write ();
    //  h_muonMed_IDEffDen_fcal[1][iSign]->Write ();

    //  h_muonID_MSEffNum_pt[0][iSign]->Write ();
    //  h_muonID_MSEffDen_pt[0][iSign]->Write ();
    //  h_muonID_MSEffNum_pt[1][iSign]->Write ();
    //  h_muonID_MSEffDen_pt[1][iSign]->Write ();
    //  h_muonID_MSEffNum_pt[2][iSign]->Write ();
    //  h_muonID_MSEffDen_pt[2][iSign]->Write ();
    //  h_muonID_MSEffNum_eta[0][iSign]->Write ();
    //  h_muonID_MSEffDen_eta[0][iSign]->Write ();
    //  h_muonID_MSEffNum_eta[1][iSign]->Write ();
    //  h_muonID_MSEffDen_eta[1][iSign]->Write ();
    //  h_muonID_MSEffNum_eta[2][iSign]->Write ();
    //  h_muonID_MSEffDen_eta[2][iSign]->Write ();
    //  h_muonID_MSEffNum_fcal[0][iSign]->Write ();
    //  h_muonID_MSEffDen_fcal[0][iSign]->Write ();
    //  h_muonID_MSEffNum_fcal[1][iSign]->Write ();
    //  h_muonID_MSEffDen_fcal[1][iSign]->Write ();
    //}

    //for (int iSign : {0, 1}) {
    //  t_muonMed_IDEff_eta_pp[iSign]->Write (Form ("t_muonMed_IDEff_eta_pp_sign%i", iSign));
    //  t_muonMed_IDEff_eta_PbPb18[iSign]->Write (Form ("t_muonMed_IDEff_eta_PbPb18_sign%i", iSign));
    //  t_muonMed_IDEff_eta_PbPb15[iSign]->Write (Form ("t_muonMed_IDEff_eta_PbPb15_sign%i", iSign));
    //  t_muonID_MSEff_eta_pp[iSign]->Write (Form ("t_muonID_MSEff_eta_pp_sign%i", iSign));
    //  t_muonID_MSEff_eta_PbPb18[iSign]->Write (Form ("t_muonID_MSEff_eta_PbPb18_sign%i", iSign));
    //  t_muonID_MSEff_eta_PbPb15[iSign]->Write (Form ("t_muonID_MSEff_eta_PbPb15_sign%i", iSign));
    //}
    t_muonIDEff_eta_phi_pp->Write ("t_muonIDEff_eta_phi_pp");
    t_muonIDEff_eta_phi_PbPb18->Write ("t_muonIDEff_eta_phi_PbPb18");
    t_muonIDEff_eta_phi_PbPb15->Write ("t_muonIDEff_eta_phi_PbPb15");

    h2_zmumuIDEffNum_pt_y[0]->Write ();
    h2_zmumuIDEffDen_pt_y[0]->Write ();
    h2_zmumuIDEffNum_pt_y[1]->Write ();
    h2_zmumuIDEffDen_pt_y[1]->Write ();
    h2_zmumuIDEffNum_pt_y[2]->Write ();
    h2_zmumuIDEffDen_pt_y[2]->Write ();

    h_electronIDEffNum_pt[0]->Write ();
    h_electronIDEffDen_pt[0]->Write ();
    h_electronIDEffNum_pt[1]->Write ();
    h_electronIDEffDen_pt[1]->Write ();
    h_electronIDEffNum_pt[2]->Write ();
    h_electronIDEffDen_pt[2]->Write ();
    h_electronIDEffNum_eta[0]->Write ();
    h_electronIDEffDen_eta[0]->Write ();
    h_electronIDEffNum_eta[1]->Write ();
    h_electronIDEffDen_eta[1]->Write ();
    h_electronIDEffNum_eta[2]->Write ();
    h_electronIDEffDen_eta[2]->Write ();
    h2_electronIDEffNum_pt_eta[0]->Write ();
    h2_electronIDEffDen_pt_eta[0]->Write ();
    h2_electronIDEffNum_pt_eta[1]->Write ();
    h2_electronIDEffDen_pt_eta[1]->Write ();
    h2_electronIDEffNum_pt_eta[2]->Write ();
    h2_electronIDEffDen_pt_eta[2]->Write ();
    h_electronIDEffNum_fcal[0]->Write ();
    h_electronIDEffDen_fcal[0]->Write ();
    h_electronIDEffNum_fcal[1]->Write ();
    h_electronIDEffDen_fcal[1]->Write ();

    t_electronIDEff_pt_eta_pp->Write ("t_electronIDEff_pt_eta_pp");
    t_electronIDEff_pt_eta_PbPb18->Write ("t_electronIDEff_pt_eta_PbPb18");
    t_electronIDEff_pt_eta_PbPb15->Write ("t_electronIDEff_pt_eta_PbPb15");
    h_electronIDEff_fcal[0]->Write ();
    h_electronIDEff_fcal[1]->Write ();

    h2_zeeIDEffNum_pt_y[0]->Write ();
    h2_zeeIDEffDen_pt_y[0]->Write ();
    h2_zeeIDEffNum_pt_y[1]->Write ();
    h2_zeeIDEffDen_pt_y[1]->Write ();
    h2_zeeIDEffNum_pt_y[2]->Write ();
    h2_zeeIDEffDen_pt_y[2]->Write ();
  }


  outFile->Close ();

  return true;
}



bool GetMuonTrackCut (const TreeVariables* t, const int muon, const string cutLevel) {
  for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
    if (t->trk_pt[iTrk] == t->muon_id_track_pt[muon] && t->trk_eta[iTrk] == t->muon_id_track_eta[muon] && InTwoPi (t->trk_phi[iTrk]) == InTwoPi (t->muon_id_track_phi[muon])) {
      if (cutLevel == "HILoose")
        return t->trk_HIloose[iTrk];
      else if (cutLevel == "HITight")
        return t->trk_HItight[iTrk];
      else
        return false;
    }
  }
  return false;
}

} // end namespace

#endif
