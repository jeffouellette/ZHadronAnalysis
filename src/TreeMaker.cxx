#ifndef __TreeMaker_cxx__
#define __TreeMaker_cxx__

#include "TreeMaker.h"
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

bool TreeMaker (const char* directory,
                const int dataSet,
                const char* inFileName) {
 
  cout << "Info: In TreeMaker.cxx: Entered TreeMaker routine." << endl;
  cout << "Info: In TreeMaker.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoElectronPtUpVar: " << doElectronPtUpVar;
  cout << "\n\tdoElectronPtDownVar: " << doElectronPtDownVar;
  cout << "\n\tdoMuonPtUpVar: " << doMuonPtUpVar;
  cout << "\n\tdoMuonPtDownVar: " << doMuonPtDownVar;
  cout << "\n\tdoElectronLHMediumVar: " << doElectronLHMediumVar;
  cout << "\n\tdoMuonTightVar: " << doMuonTightVar;
  cout << "\n\tdoHITightVar: " << doHITightVar;
  cout << endl;

  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  if (isHijing)
    cout << "Info: In TreeMaker.cxx: File detected as Hijing overlay" << endl;
  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TreeMaker.cxx: File detected as data overlay, will check data conditions" << endl;

  if (isMC && isHijing)
    SetupDirectories ("MCAnalysis/Hijing", false);
  else if (isMC)
    SetupDirectories ("MCAnalysis");
  else
    SetupDirectories ("DataAnalysis");
    

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TreeMaker.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TreeMaker.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In TreeMaker.C: Cannot identify this MC file! Quitting." << endl;
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
    cout << "Error: In TreeMaker.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }




  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In TreeMaker.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }




  TFile* f_pp_muonTrigEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/muontrigger_sf_2017_mc16d_v01.root", "read");
  TH2F* h2_muonTrigEff_barrel_eta_phi_pp = (TH2F*) f_pp_muonTrigEff->Get ("Medium/PeriodK/HLT_mu14/eff_etaphi_fine_barrel_data_nominal")->Clone ("h2_muonTrigEff_barrel_eta_phi_pp");
  TH2F* h2_muonTrigEff_endcap_eta_phi_pp = (TH2F*) f_pp_muonTrigEff->Get ("Medium/PeriodK/HLT_mu14/eff_etaphi_fine_endcap_data_nominal")->Clone ("h2_muonTrigEff_endcap_eta_phi_pp");
  //f_pp_muonTrigEff->Close ();
  if (!h2_muonTrigEff_barrel_eta_phi_pp || !h2_muonTrigEff_barrel_eta_phi_pp)
    cout << "Error: In TreeMaker.cxx: Muon trigger efficiency in pp not found!" << endl;

  TFile* f_muonIDEff_pp = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/Reco_Medium_Z.root", "read");
  TH3D* h3_muonIDEff_eta_phi_pt_pp = (TH3D*) f_muonIDEff_pp->Get ("Eff_2017")->Clone ("h3_muonIDEff_eta_phi_pt_pp");
  //f_muonIDEff_pp->Close ();
  if (!h3_muonIDEff_eta_phi_pt_pp)
    cout << "Error: In TreeMaker.cxx: Muon reconstruction efficiency in pp not found!" << endl;

  TFile* f_myEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/Nominal/outFile.root", "read");

  TH2D* h2_muonTrigEff_eta_phi_PbPb18 = (TH2D*) f_myEff->Get ("h2_muonTrigEffNum_eta_phi_PbPb18")->Clone ("h2_muonTrigEff_eta_phi_PbPb18");
  h2_muonTrigEff_eta_phi_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_muonTrigEffDen_eta_phi_PbPb18"));
  if (!h2_muonTrigEff_eta_phi_PbPb18)
    cout << "Error: In TreeMaker.cxx: Muon trigger efficiency in PbPb (2018) not found!" << endl;

  TH2D* h2_muonIDEff_eta_phi_PbPb18 = (TH2D*) f_myEff->Get ("h2_muonIDEffNum_eta_phi_PbPb18")->Clone ("h2_muonIDEff_eta_phi_PbPb18");
  h2_muonIDEff_eta_phi_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_muonIDEffDen_eta_phi_PbPb18"));
  TH1D* h_muonIDEff_fcal_PbPb18 = (TH1D*) f_myEff->Get ("h_muonIDEffNum_fcal_PbPb18")->Clone ("h_muonIDEff_fcal_PbPb18");
  h_muonIDEff_fcal_PbPb18->Divide ((TH1D*) f_myEff->Get ("h_muonIDEffDen_fcal_PbPb18"));
  if (!h2_muonIDEff_eta_phi_PbPb18 || !h_muonIDEff_fcal_PbPb18)
    cout << "Error: In TreeMaker.cxx: Muon reconstruction efficiency in PbPb (2018) not found!" << endl;
  TF1* f_muonIDEff_fcal_PbPb18 = new TF1 ("f_muonIDeff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_muonIDEff_fcal_PbPb18->SetParameter (0, 1);
  f_muonIDEff_fcal_PbPb18->SetParameter (1, 0);
  h_muonIDEff_fcal_PbPb18->Fit (f_muonIDEff_fcal_PbPb18, "RN0Q");

  TH2D* h2_electronTrigEff_pt_eta_pp = (TH2D*) f_myEff->Get ("h2_electronTrigEffNum_pt_eta_pp")->Clone ("h2_electronTrigEff_pt_eta_pp");
  h2_electronTrigEff_pt_eta_pp->Divide ((TH2D*) f_myEff->Get ("h2_electronTrigEffDen_pt_eta_pp"));
  if (!h2_electronTrigEff_pt_eta_pp)
    cout << "Error: In TreeMaker.cxx: Electron trigger efficiency in pp not found!" << endl;

  TH2D* h2_electronIDEff_pt_eta_pp = (TH2D*) f_myEff->Get ("h2_electronIDEffNum_pt_eta_pp")->Clone ("h2_electronIDEff_pt_eta_pp");
  h2_electronIDEff_pt_eta_pp->Divide ((TH2D*) f_myEff->Get ("h2_electronIDEffDen_pt_eta_pp"));
  if (!h2_electronIDEff_pt_eta_pp)
    cout << "Error: In TreeMaker.cxx: Electron reconstruction efficiency in pp not found!" << endl;

  TH2D* h2_electronTrigEff_pt_eta_PbPb18 = (TH2D*) f_myEff->Get ("h2_electronTrigEffNum_pt_eta_PbPb18")->Clone ("h2_electronTrigEff_pt_eta_PbPb18");
  h2_electronTrigEff_pt_eta_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_electronTrigEffDen_pt_eta_PbPb18"));
  TH1D* h_electronTrigEff_fcal_PbPb18 = (TH1D*) f_myEff->Get ("h_electronTrigEffNum_fcal_PbPb18")->Clone ("h_electronTrigEff_fcal_PbPb18");
  h_electronTrigEff_fcal_PbPb18->Divide ((TH1D*) f_myEff->Get ("h_electronTrigEffDen_fcal_PbPb18"));
  if (!h2_electronTrigEff_pt_eta_PbPb18 || !h_electronTrigEff_fcal_PbPb18)
    cout << "Error: In TreeMaker.cxx: Electron trigger efficiency in PbPb (2018) not found!" << endl;
  TF1* f_electronTrigEff_fcal_PbPb18 = new TF1 ("f_electronTrigeff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_electronTrigEff_fcal_PbPb18->SetParameter (0, 1);
  f_electronTrigEff_fcal_PbPb18->SetParameter (1, 0);
  h_electronTrigEff_fcal_PbPb18->Fit (f_electronTrigEff_fcal_PbPb18, "RN0Q");

  TH2D* h2_electronIDEff_pt_eta_PbPb18 = (TH2D*) f_myEff->Get ("h2_electronIDEffNum_pt_eta_PbPb18")->Clone ("h2_electronIDEff_pt_eta_PbPb18");
  h2_electronIDEff_pt_eta_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_electronIDEffDen_pt_eta_PbPb18"));
  TH1D* h_electronIDEff_fcal_PbPb18 = (TH1D*) f_myEff->Get ("h_electronIDEffNum_fcal_PbPb18")->Clone ("h_electronIDEff_fcal_PbPb18");
  h_electronIDEff_fcal_PbPb18->Divide ((TH1D*) f_myEff->Get ("h_electronIDEffDen_fcal_PbPb18"));
  if (!h2_electronIDEff_pt_eta_PbPb18 || !h_electronIDEff_fcal_PbPb18)
    cout << "Error: In TreeMaker.cxx: Electron reconstruction efficiency in PbPb (2018) not found!" << endl;
  //TF1* f_electronIDEff_fcal_PbPb18 = new TF1 ("f_electronIDeff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  //f_electronIDEff_fcal_PbPb18->SetParameter (0, 1);
  //f_electronIDEff_fcal_PbPb18->SetParameter (1, 0);
  //h_electronIDEff_fcal_PbPb18->Fit (f_electronIDEff_fcal_PbPb18, "RN0Q");


  // Input tree
  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetMuons ();
  t->SetGetTracks ();
  if (!isMC) {
    if (!isPbPb) t->SetGetJets ();
    else if (isPbPb) t->SetGetZdc ();
  }
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

  Trigger* electronTrigger = nullptr, *muonTrigger = nullptr;
  std::string electronTrigName, muonTrigName;
  if (collisionSystem == PbPb15) {
    electronTrigName = "HLT_e15_loose_ion_L1EM12";
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
    cout << "Error: In TreeMaker.cxx: Invalid collision system, quitting." << endl;
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
      break;
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
      break;
    }
  }


  // Load files for output
  const int numFileBins = (isHijing ? numFileIPBins : numFileCentBins);
  TFile* outFiles[numFileBins][numRunGroups];
  OutTree* outTrees[numFileBins][numRunGroups];
  const char* outTreeName = (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree");
  if (isPbPb) {
    const short iRG = (isMC ? -1: GetRunGroup (dataSet));
    for (int iCent = 1; iCent < numFileBins; iCent++) {
      if (iRG != -1) {
        outFiles[iCent][iRG] = new TFile (Form ("%s/%s/%s_iCent%i.root", rootPath.Data (), GetRunGroupTString (dataSet).Data (), identifier.Data (), iCent), "recreate");
        outFiles[iCent][iRG]->Delete (Form ("%s;*", outTreeName));
        outTrees[iCent][iRG] = new OutTree (outTreeName, outFiles[iCent][iRG]);
        outTrees[iCent][iRG]->SetBranchEventInfo ();
        outTrees[iCent][iRG]->SetBranchLeptons ();
        outTrees[iCent][iRG]->SetBranchZs ();
        if (!isPbPb) outTrees[iCent][iRG]->SetBranchJets ();
        outTrees[iCent][iRG]->SetBranchTracks ();
        outTrees[iCent][iRG]->SetBranches ();
      }
      else if (isHijing) {
        outFiles[iCent][0] = new TFile (Form ("%s/%s/%s_iCent%i.root", rootPath.Data (), identifier.Data (), identifier.Data (), iCent), "recreate");
        outFiles[iCent][0]->Delete (Form ("%s;*", outTreeName));
        outTrees[iCent][0] = new OutTree (outTreeName, outFiles[iCent][0]);
        outTrees[iCent][0]->SetBranchEventInfo ();
        outTrees[iCent][0]->SetBranchLeptons ();
        outTrees[iCent][0]->SetBranchZs ();
        if (!isPbPb) outTrees[iCent][0]->SetBranchJets ();
        outTrees[iCent][0]->SetBranchTracks ();
        outTrees[iCent][0]->SetBranches ();
        outTrees[iCent][0]->tree->Branch ("impactParameter", &ip, "impactParameter/F");
        outTrees[iCent][0]->tree->Branch ("eventPlane", &eventPlane, "eventPlane/F");
      }
      else {
        for (short iRGMC = 0; iRGMC < numRunGroups; iRGMC++) {
          outFiles[iCent][iRGMC] = new TFile (Form ("%s/%s/%s_iCent%i.root", rootPath.Data (), identifier.Data (), runGroups[iRGMC].first.c_str (), iCent), "recreate");
          outFiles[iCent][iRGMC]->Delete (Form ("%s;*", outTreeName));
          outTrees[iCent][iRGMC] = new OutTree (outTreeName, outFiles[iCent][iRGMC]);
          outTrees[iCent][iRGMC]->SetBranchEventInfo ();
          outTrees[iCent][iRGMC]->SetBranchLeptons ();
          outTrees[iCent][iRGMC]->SetBranchZs ();
          if (!isPbPb) outTrees[iCent][iRGMC]->SetBranchJets ();
          outTrees[iCent][iRGMC]->SetBranchTracks ();
          outTrees[iCent][iRGMC]->SetBranches ();
        }
      }
    }
  }
  else {
    outFiles[0][0] = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
    outFiles[0][0]->Delete (Form ("%s;*", outTreeName));
    outTrees[0][0] = new OutTree (outTreeName, outFiles[0][0]);
    outTrees[0][0]->SetBranchEventInfo ();
    outTrees[0][0]->SetBranchLeptons ();
    outTrees[0][0]->SetBranchZs ();
    if (!isPbPb) outTrees[0][0]->SetBranchJets ();
    outTrees[0][0]->SetBranchTracks ();
    outTrees[0][0]->SetBranches ();
  }


  const int nEvts = tree->GetEntries ();
  TLorentzVector l1, l2;
  float trigEff1, trigEff2, recoEff1, recoEff2, trigEff, recoEff;

  // First loop over events looking for Z->ee candidates
  isEE = true;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TreeMaker.cxx: Electrons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;

    event_number = t->event_number;
    run_number = t->run_number;

    if (!isHijing) {
      bool hasPrimary = false;
      //bool hasPileup = false;
      vz = -999;
      for (int iVert = 0; iVert < t->nvert; iVert++) {
        const bool isPrimary = (t->vert_type[iVert] == 1);
        hasPrimary = hasPrimary || isPrimary;
        //hasPileup = hasPileup || (t->vert_type[iVert] == 3);
        if (isPrimary)
          vz = t->vert_z[iVert];
      }
      //if (!hasPrimary || (!isPbPb && hasPileup) || fabs (vz) > 150)
      if (!hasPrimary || fabs (vz) > 150)
        continue;
    }
    else {
      if (t->nvert > 0) vz = t->vert_z[0];
      else vz = 0;
    }

    event_weight = -1;
    if (!isMC) {
      if (electronTrigger->trigBool)
        event_weight = electronTrigger->trigPrescale;
    }
    else {
      event_weight = crossSectionPicoBarns * mcFilterEfficiency * (isPbPb ? 0.00171786 : 258.4); // sigma * f * L_int
      //event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency * 0.00171786; // sigma * f * L_int
      if (isHijing) {
        assert (t->nTruthEvt > 0);
        ip = t->impactParameter[0];
        eventPlane = t->eventPlaneAngle[0];
      }
    }
    if (event_weight == -1)
      continue; // trigger requirement

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

    if (isPbPb) {
      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

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
    short iRG = (isPbPb && !isHijing ? GetRunGroup (run_number) : 0);

    const float _event_weight = event_weight; // for storage, in case of multiple Z's in one event

    for (int iE1 = 0; iE1 < t->electron_n; iE1++) {
      l1_pt = t->electron_pt[iE1];
      if (doElectronPtDownVar)    l1_pt = t->electron_pt[iE1] - t->electron_pt_sys[iE1];
      else if (doElectronPtUpVar) l1_pt = t->electron_pt[iE1] + t->electron_pt_sys[iE1];
      l1_eta = t->electron_eta[iE1];
      l1_phi = t->electron_phi[iE1];
      l1_charge = t->electron_charge[iE1];
      l1_trk_pt = t->electron_id_track_pt[iE1];
      l1_trk_eta = t->electron_id_track_eta[iE1];
      l1_trk_phi = t->electron_id_track_phi[iE1];
      l1_pt_sys = t->electron_pt_sys[iE1];
      l1_eta_sys = t->electron_eta_sys[iE1];
      l1_phi_sys = t->electron_phi_sys[iE1];

      if (isPbPb) {
        l1_pt *= GetZmassSF_MC (fcal_et, l1_eta);
        if (!isMC)
          l1_pt *= GetZmassSF_PbPb (fcal_et, l1_eta);
      }

      if (l1_pt < electron_pt_cut)
        continue; // basic electron pT cut
      if (!InEMCal (l1_eta))
        continue; // reject electrons reconstructed outside the EMCal

      if (!doElectronLHMediumVar) {
        if (!isPbPb && !t->electron_lhloose[iE1])
          continue; // reject non-loose electrons
        else if (isPbPb && !t->electron_lhloose_hi[iE1])
          continue; // alternative for HI collisions
      }
      else {
        if (!isPbPb && !t->electron_lhmedium[iE1])
          continue; // reject non-medium electrons
        else if (isPbPb && !t->electron_lhmedium_hi[iE1])
          continue; // alternative for HI collisions
      }

      if (fabs (t->electron_id_track_d0sig[iE1]) > 5)
        continue; // electron d0 significance vertex compatibility cut
      l1_d0sig = t->electron_id_track_d0sig[iE1];

      if (fabs ( (t->electron_id_track_z0[iE1] + t->electron_id_track_vz[iE1] - vz) * sin (t->electron_id_track_theta[iE1])) > 0.5)
        continue; // electron z0 vertex compatibility cut
      l1_z0 = t->electron_id_track_z0[iE1] + t->electron_id_track_vz[iE1]; // z0 in detector coordinates

      for (int iE2 = 0; iE2 < iE1; iE2++) {
        l2_pt = t->electron_pt[iE2];
        if (doElectronPtDownVar)    l2_pt = t->electron_pt[iE2] - t->electron_pt_sys[iE2];
        else if (doElectronPtUpVar) l2_pt = t->electron_pt[iE2] + t->electron_pt_sys[iE2];
        l2_eta = t->electron_eta[iE2];
        l2_phi = t->electron_phi[iE2];
        l2_charge = t->electron_charge[iE2];
        l2_trk_pt = t->electron_id_track_pt[iE2];
        l2_trk_eta = t->electron_id_track_eta[iE2];
        l2_trk_phi = t->electron_id_track_phi[iE2];
        l2_pt_sys = t->electron_pt_sys[iE2];
        l2_eta_sys = t->electron_eta_sys[iE2];
        l2_phi_sys = t->electron_phi_sys[iE2];

        if (isPbPb) {
          l2_pt *= GetZmassSF_MC (fcal_et, l2_eta);
          if (!isMC)
            l2_pt *= GetZmassSF_PbPb (fcal_et, l2_eta);
        }

        if (l2_pt < electron_pt_cut)
          continue; // basic electron pT cut
        if (!InEMCal (l2_eta))
          continue; // reject electrons reconstructed outside the EMCal

        if (!doElectronLHMediumVar) {
          if (!isPbPb && !t->electron_lhloose[iE2])
            continue; // reject non-loose electrons
          else if (isPbPb && !t->electron_lhloose_hi[iE2])
            continue; // alternative for HI collisions
        }
        else {
          if (!isPbPb && !t->electron_lhmedium[iE2])
            continue; // reject non-medium electrons
          else if (isPbPb && !t->electron_lhmedium_hi[iE2])
            continue; // alternative for HI collisions
        }

        if (fabs (t->electron_id_track_d0sig[iE2]) > 5)
          continue; // electron d0 significance vertex compatibility cut
        l2_d0sig = t->electron_id_track_d0sig[iE2];

        if (fabs ( (t->electron_id_track_z0[iE2] + t->electron_id_track_vz[iE2] - vz) * sin (t->electron_id_track_theta[iE2])) > 0.5)
          continue; // electron z0 vertex compatibility cut
        l2_z0 = t->electron_id_track_z0[iE2] + t->electron_id_track_vz[iE2]; // z0 in detector coordinates

        l1.SetPtEtaPhiM (t->electron_pt[iE1], l1_eta, l1_phi, electron_mass);
        l2.SetPtEtaPhiM (t->electron_pt[iE2], l2_eta, l2_phi, electron_mass);

        const float nom_z_pt = (l1+l2).Pt ();
        const float nom_z_y = (l1+l2).Rapidity ();

        l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, electron_mass);
        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, electron_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged electrons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        if (!isMC && !t->electron_matched[iE1] && !t->electron_matched[iE2])
          continue;

        if (!isPbPb) {
          TH2D* h = h2_electronTrigEff_pt_eta_pp;
          trigEff1 = (!isMC ? h->GetBinContent (h->FindFixBin (l1_eta, fmin (l1_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ())))) : 1);
          trigEff2 = (!isMC ? h->GetBinContent (h->FindFixBin (l2_eta, fmin (l2_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ())))) : 1);
          h = h2_electronIDEff_pt_eta_pp;
          recoEff1 = (h->GetBinContent (h->FindFixBin (l1_eta, fmin (l1_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ())))));
          recoEff2 = (h->GetBinContent (h->FindFixBin (l2_eta, fmin (l2_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ())))));
        }
        else {
          float fcalEff = (isMC ? 1 : f_electronTrigEff_fcal_PbPb18->Eval (fcal_et * 1e-3));
          TH2D* h = h2_electronTrigEff_pt_eta_PbPb18;
          trigEff1 = (isMC ? 1 : h->GetBinContent (h->FindFixBin (l1_eta, fmin (l1_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ()))))) * fcalEff;
          trigEff2 = (isMC ? 1 : h->GetBinContent (h->FindFixBin (l2_eta, fmin (l2_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ()))))) * fcalEff;

          fcalEff = (isMC ? 1 : h_electronIDEff_fcal_PbPb18->GetBinContent (h_electronIDEff_fcal_PbPb18->FindFixBin (fcal_et * 1e-3)));
          h = h2_electronIDEff_pt_eta_PbPb18;
          recoEff1 = h->GetBinContent (h->FindFixBin (l1_eta, fmin (l1_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ())))) * fcalEff;
          recoEff2 = h->GetBinContent (h->FindFixBin (l2_eta, fmin (l2_pt, h->GetYaxis ()->GetBinCenter (h->GetYaxis ()->GetNbins ())))) * fcalEff;
        }
        trigEff = (1.-(1.-trigEff1)*(1.-trigEff2));
        recoEff = (recoEff1)*(recoEff2);

        if (trigEff == 0) {
          cout << "trigEff == 0" << endl; 
          continue;
        }
        if (recoEff == 0) {
          cout << "recoEff == 0" << endl;
          continue;
        }
        event_weight = _event_weight / (trigEff * recoEff);

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

          if (IsElectronTrack (t, iTrk, iE1, iE2))
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

        if (sumptccw > sumptcw) {
          phi_transmax = z_phi + pi/2.;
          phi_transmin = z_phi - pi/2.;
        }
        else {
          phi_transmax = z_phi - pi/2.;
          phi_transmin = z_phi + pi/2.;
        }

        if (!isPbPb) {
          njet = 0;
          for (int iJet = 0; iJet < t->akt4emtopo_jet_n; iJet++) {
            float _jpt = t->akt4emtopo_jet_pt[iJet];
            float _jeta = t->akt4emtopo_jet_eta[iJet];
            float _jphi = t->akt4emtopo_jet_phi[iJet];
            float _je = t->akt4emtopo_jet_e[iJet];

            if (_jpt < 15)
              continue;
            if (DeltaR (_jeta, l1_eta, _jphi, l1_phi) < 0.2)
              continue;
            if (DeltaR (_jeta, l2_eta, _jphi, l2_phi) < 0.2)
              continue;

            jet_pt[njet] = _jpt;
            jet_eta[njet] = _jeta;
            jet_phi[njet] = _jphi;
            jet_e[njet] = _je;
            njet++;
          }
        }

        ntrk_perp *= 2;
        //if (isPbPb && isHijing) {
        //  iCent = GetFileNtrkBin (ntrk_perp);
        //  if (isPbPb && (iCent < 1 || iCent > numFileBins-1))
        //    continue;
        //}

        outTrees[iCent][iRG]->Fill ();
        
      } // end iE2 loop
    } // end iE1 loop

  } // end electron selection
  cout << endl << "Info: In TreeMaker.cxx: Finished processing electrons." << endl;


  // Now loop over events looking for Z->mumu candidates
  isEE = false;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TreeMaker.cxx: Muons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    event_number = t->event_number;
    run_number = t->run_number;

    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;
    if ((!isMC || isOverlayMC) && !t->passes_toroid)
      continue; // additional check for muon quality

    if (!isHijing) {
      bool hasPrimary = false;
      //bool hasPileup = false;
      vz = -999;
      for (int iVert = 0; iVert < t->nvert; iVert++) {
        const bool isPrimary = (t->vert_type[iVert] == 1);
        hasPrimary = hasPrimary || isPrimary;
        //hasPileup = hasPileup || (t->vert_type[iVert] == 3);
        if (isPrimary)
          vz = t->vert_z[iVert];
      }
      //if (!hasPrimary || (!isPbPb && hasPileup) || fabs (vz) > 150)
      if (!hasPrimary || fabs (vz) > 150)
        continue;
    }
    else {
      if (t->nvert > 0) vz = t->vert_z[0];
      else vz = 0;
    }

    event_weight = -1;
    if (!isMC) {
      if (muonTrigger->trigBool)
        event_weight = muonTrigger->trigPrescale;
    }
    else {
      event_weight = crossSectionPicoBarns * mcFilterEfficiency * (isPbPb ? 0.00143844 : 258.4); // sigma * f * L_int
      //event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency * (isPbPb ? 0.00143844 : 258.4); // sigma * f * L_int
      if (isHijing) {
        assert (t->nTruthEvt > 0);
        ip = t->impactParameter[0];
        eventPlane = t->eventPlaneAngle[0];
      }
    }
    if (event_weight == -1)
      continue; // trigger requirement

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

    if (isPbPb) {
      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

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
    short iRG = (isPbPb && !isHijing ? GetRunGroup (run_number) : 0);

    const float _event_weight = event_weight; // for storage, in case of multiple Z's in one event

    for (int iM1 = 0; iM1 < t->muon_n; iM1++) {
      l1_pt = (isPbPb ? t->muon_ms_pt[iM1] : t->muon_pt[iM1]);
      l1_eta = t->muon_eta[iM1];
      l1_phi = t->muon_phi[iM1];
      l1_charge = t->muon_charge[iM1];
      l1_trk_pt = t->muon_id_track_pt[iM1];
      l1_trk_eta = t->muon_id_track_eta[iM1];
      l1_trk_phi = t->muon_id_track_phi[iM1];

      //int iMaxSys1 = 0;
      //for (int iSys = 0; iSys < (int)(t->muon_pt_sys->at (iM1).size ()); iSys++)
      //  if (fabs (t->muon_pt_sys->at (iM1).at (iSys) - l1_pt) > fabs (t->muon_pt_sys->at (iM1).at (iMaxSys1) - l1_pt))
      //    iMaxSys1 = iSys;

      //l1_pt_sys = t->muon_pt_sys->at (iM1).at (iMaxSys1);
      //l1_eta_sys = t->muon_eta_sys->at (iM1).at (iMaxSys1);
      //l1_phi_sys = t->muon_phi_sys->at (iM1).at (iMaxSys1);

      //if (doMuonPtDownVar) l1_pt = l1_pt - fabs (l1_pt - l1_pt_sys);
      //else if (doMuonPtUpVar)   l1_pt = l1_pt + fabs (l1_pt - l1_pt_sys);

      if (l1_pt < muon_pt_cut)
        continue; // basic muon pT cut
      if (fabs (l1_eta) > 2.5)
        continue; // reject muons reconstructed outside the detector

      if (!doMuonTightVar && !t->muon_medium[iM1])
        continue; // reject non-medium muons
      else if (doMuonTightVar && !t->muon_tight[iM1])
        continue; // reject non-tight muons

      if (fabs (t->muon_id_track_d0sig[iM1]) > 3)
        continue; // muon d0 cut
      l1_d0sig = t->muon_id_track_d0sig[iM1];

      if (fabs ( (t->muon_id_track_z0[iM1] + t->muon_id_track_vz[iM1] - vz) * sin (t->muon_id_track_theta[iM1])) > 0.5)
        continue; // muon z0 to vertex cut
      l1_z0 = t->muon_id_track_z0[iM1] + t->muon_id_track_vz[iM1]; // z0 in detector coordinates

      for (int iM2 = 0; iM2 < iM1; iM2++) {
        l2_pt = (isPbPb ? t->muon_ms_pt[iM2] : t->muon_pt[iM2]);
        l2_eta = t->muon_eta[iM2];
        l2_phi = t->muon_phi[iM2];
        l2_charge = t->muon_charge[iM2];
        l2_trk_pt = t->muon_id_track_pt[iM2];
        l2_trk_eta = t->muon_id_track_eta[iM2];
        l2_trk_phi = t->muon_id_track_phi[iM2];

        //int iMaxSys2 = 0;
        //for (int iSys = 0; iSys < (int)(t->muon_pt_sys->at (iM2).size ()); iSys++)
        //  if (fabs (t->muon_pt_sys->at (iM2).at (iSys) - l2_pt) > fabs (t->muon_pt_sys->at (iM2).at (iMaxSys2) - l2_pt))
        //    iMaxSys2 = iSys;

        //l2_pt_sys = t->muon_pt_sys->at (iM2).at (iMaxSys2);
        //l2_eta_sys = t->muon_eta_sys->at (iM2).at (iMaxSys2);
        //l2_phi_sys = t->muon_phi_sys->at (iM2).at (iMaxSys2);

        //if (doMuonPtDownVar) l2_pt = l2_pt - fabs (l2_pt - l2_pt_sys);
        //else if (doMuonPtUpVar)   l2_pt = l2_pt + fabs (l2_pt - l2_pt_sys);

        if (l2_pt < muon_pt_cut)
          continue; // basic muon pT cut
        if (fabs (l2_eta) > 2.5)
          continue; // reject muons reconstructed outside the detector

        if (!doMuonTightVar && !t->muon_medium[iM2])
          continue; // reject non-medium muons
        else if (doMuonTightVar && !t->muon_tight[iM2])
          continue; // reject non-tight muons

        if (fabs (t->muon_id_track_d0sig[iM2]) > 3)
          continue; // muon d0 significance vertex compatibility cut
        l2_d0sig = t->muon_id_track_d0sig[iM2];

        if (fabs ( (t->muon_id_track_z0[iM2] + t->muon_id_track_vz[iM2] - vz) * sin (t->muon_id_track_theta[iM2])) > 0.5)
          continue; // muon z0 vertex compatibility cut
        l2_z0 = t->muon_id_track_z0[iM2] + t->muon_id_track_vz[iM2]; // z0 in detector coordinates

        l1.SetPtEtaPhiM (isPbPb ? t->muon_ms_pt[iM1] : t->muon_pt[iM1], l1_eta, l1_phi, muon_mass);
        l2.SetPtEtaPhiM (isPbPb ? t->muon_ms_pt[iM2] : t->muon_pt[iM2], l2_eta, l2_phi, muon_mass);

        const float nom_z_pt = (l1+l2).Pt ();
        const float nom_z_y = (l1+l2).Rapidity ();

        l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, muon_mass);
        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, muon_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged muons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        if (!isMC && !t->muon_matched[iM1] && !t->muon_matched[iM2])
          continue;

        if (!isPbPb) {
          const float _l1_phi = (InTwoPi (l1_phi) > pi ? InTwoPi (l1_phi) - 2*pi : InTwoPi (l1_phi));
          const float _l2_phi = (InTwoPi (l2_phi) > pi ? InTwoPi (l2_phi) - 2*pi : InTwoPi (l2_phi));

          TH2F* th2f = (fabs (l1_eta) < 1.05 ? h2_muonTrigEff_barrel_eta_phi_pp : h2_muonTrigEff_endcap_eta_phi_pp); 
          trigEff1 = (isMC ? 1 : th2f->GetBinContent (th2f->FindFixBin (l1_eta, _l1_phi)));
          th2f = (fabs (l2_eta) < 1.05 ? h2_muonTrigEff_barrel_eta_phi_pp : h2_muonTrigEff_endcap_eta_phi_pp); 
          trigEff2 = (isMC ? 1 : th2f->GetBinContent (th2f->FindFixBin (l2_eta, _l2_phi)));

          TH3D* th3d = h3_muonIDEff_eta_phi_pt_pp;
          recoEff1 = th3d->GetBinContent (th3d->FindFixBin (l1_eta, _l1_phi, l1_pt));
          recoEff2 = th3d->GetBinContent (th3d->FindFixBin (l2_eta, _l2_phi, l2_pt));
        }
        else {
          trigEff1 = (isMC ? 1 : h2_muonTrigEff_eta_phi_PbPb18->GetBinContent (h2_muonTrigEff_eta_phi_PbPb18->FindFixBin (l1_eta, InTwoPi (l1_phi))));
          trigEff2 = (isMC ? 1 : h2_muonTrigEff_eta_phi_PbPb18->GetBinContent (h2_muonTrigEff_eta_phi_PbPb18->FindFixBin (l2_eta, InTwoPi (l2_phi))));

          const float fcalEff = (isMC ? 1 : f_muonIDEff_fcal_PbPb18->Eval (fcal_et * 1e-3));
          recoEff1 = fcalEff * h2_muonIDEff_eta_phi_PbPb18->GetBinContent (h2_muonIDEff_eta_phi_PbPb18->FindFixBin (l1_eta, InTwoPi (l1_phi)));
          recoEff2 = fcalEff * h2_muonIDEff_eta_phi_PbPb18->GetBinContent (h2_muonIDEff_eta_phi_PbPb18->FindFixBin (l2_eta, InTwoPi (l2_phi)));
        }

        trigEff = (1.-(1.-trigEff1)*(1.-trigEff2));
        recoEff = (recoEff1)*(recoEff2);

        if (trigEff == 0) {
          cout << "trigEff == 0" << endl; 
          continue;
        }
        if (recoEff == 0) {
          cout << "recoEff == 0" << endl;
          continue;
        }
        event_weight = _event_weight / (trigEff * recoEff);

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

          if (IsMuonTrack (t, iTrk, iM1, iM2))
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

        if (sumptccw > sumptcw) {
          phi_transmax = z_phi + pi/2.;
          phi_transmin = z_phi - pi/2.;
        }
        else {
          phi_transmax = z_phi - pi/2.;
          phi_transmin = z_phi + pi/2.;
        }

        if (!isPbPb) {
          njet = 0;
          for (int iJet = 0; iJet < t->akt4emtopo_jet_n; iJet++) {
            float _jpt = t->akt4emtopo_jet_pt[iJet];
            float _jeta = t->akt4emtopo_jet_eta[iJet];
            float _jphi = t->akt4emtopo_jet_phi[iJet];
            float _je = t->akt4emtopo_jet_e[iJet];

            if (_jpt < 15)
              continue;
            if (DeltaR (_jeta, l1_eta, _jphi, l1_phi) < 0.2)
              continue;
            if (DeltaR (_jeta, l2_eta, _jphi, l2_phi) < 0.2)
              continue;

            jet_pt[njet] = _jpt;
            jet_eta[njet] = _jeta;
            jet_phi[njet] = _jphi;
            jet_e[njet] = _je;
            njet++;
          }
        }

        ntrk_perp *= 2;
        //if (isPbPb && isHijing) {
        //  iCent = GetFileNtrkBin (ntrk_perp);
        //  if (isPbPb && (iCent < 1 || iCent > numFileBins-1))
        //    continue;
        //}

        outTrees[iCent][iRG]->Fill ();
        
      } // end iM2 loop
    } // end iM1 loop

  } // end muon selection
  cout << endl << "Info: In TreeMaker.cxx: Finished processing muons." << endl;

  SaferDelete (&electronTrigger);
  SaferDelete (&muonTrigger);

  SaferDelete (&t);

  if (isPbPb) {
    const short iRG = (isMC ? -1: GetRunGroup (dataSet));
    for (int iCent = 1; iCent < numFileBins; iCent++) {
      if (iRG != -1) {
        outFiles[iCent][iRG]->Write (0, TObject::kOverwrite);
        outFiles[iCent][iRG]->Close ();
        SaferDelete (&(outFiles[iCent][iRG]));
      }
      else if (isHijing) {
        outFiles[iCent][0]->Write (0, TObject::kOverwrite);
        outFiles[iCent][0]->Close ();
        SaferDelete (&(outFiles[iCent][0]));
      }
      else {
        for (int iRGMC = 0; iRGMC < numRunGroups; iRGMC++) {
          outFiles[iCent][iRGMC]->Write (0, TObject::kOverwrite);
          outFiles[iCent][iRGMC]->Close ();
          SaferDelete (&(outFiles[iCent][iRGMC]));
        }
      }
    }
  }
  else {
    outFiles[0][0]->Write ();//0, TObject::kOverwrite);
    outFiles[0][0]->Close ();
    SaferDelete (&(outFiles[0][0]));
  }

  SaferDelete (&h_zdcCuts);

  return true;
}

} // end namespace

#endif
