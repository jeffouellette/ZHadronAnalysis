#ifndef __BkgEstimator_cxx__
#define __BkgEstimator_cxx__

#include "BkgEstimator.h"
#include "Params.h"
#include "TreeVariables.h"
#include "OutTree.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"
#include "Cuts.h"

#include <Utilities.h>
#include <AtlasUtils.h>

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

bool BkgEstimator (const char* directory,
                   const int dataSet,
                   const char* inFileName) {
 
  cout << "Info: In BkgEstimator.cxx: Entered BkgEstimator routine." << endl;
  cout << "Info: In BkgEstimator.cxx: No configurable systematics available";
  cout << endl;

  SetupDirectories ("BkgEstimator");

  const bool isHijing = (isMC && strstr (inFileName, "PbPb") != NULL && strstr(inFileName, "Hijing") != NULL);
  if (isHijing)
    cout << "Info: In BkgEstimator.cxx: File detected as Hijing overlay" << endl;
  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In BkgEstimator.cxx: File detected as data overlay, will check data conditions" << endl;
    

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In BkgEstimator.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In BkgEstimator.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In BkgEstimator.C: Cannot identify this MC file! Quitting." << endl;
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
    cout << "Error: In BkgEstimator.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }




  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In BkgEstimator.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }




  TFile* f_pp_muonTrigEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/muontrigger_sf_2017_mc16d_v01.root", "read");
  TH2F* h2_muonTrigEff_barrel_eta_phi_pp = (TH2F*) f_pp_muonTrigEff->Get ("Medium/PeriodK/HLT_mu14/eff_etaphi_fine_barrel_data_nominal")->Clone ("h2_muonTrigEff_barrel_eta_phi_pp");
  TH2F* h2_muonTrigEff_endcap_eta_phi_pp = (TH2F*) f_pp_muonTrigEff->Get ("Medium/PeriodK/HLT_mu14/eff_etaphi_fine_endcap_data_nominal")->Clone ("h2_muonTrigEff_endcap_eta_phi_pp");
  //f_pp_muonTrigEff->Close ();

  TFile* f_muonIDEff_pp = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/Reco_Medium_Z.root", "read");
  TH3D* h3_muonIDEff_eta_phi_pt_pp = (TH3D*) f_muonIDEff_pp->Get ("Eff_2017")->Clone ("h3_muonIDEff_eta_phi_pt_pp");
  //f_muonIDEff_pp->Close ();

  TFile* f_myEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/Nominal/outFile.root", "read");

  TH2D* h2_muonTrigEff_eta_phi_PbPb18 = (TH2D*) f_myEff->Get ("h2_muonTrigEffNum_eta_phi_PbPb18")->Clone ("h2_muonTrigEff_eta_phi_PbPb18");
  h2_muonTrigEff_eta_phi_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_muonTrigEffDen_eta_phi_PbPb18"));

  TH2D* h2_muonIDEff_eta_phi_PbPb18 = (TH2D*) f_myEff->Get ("h2_muonIDEffNum_eta_phi_PbPb18")->Clone ("h2_muonIDEff_eta_phi_PbPb18");
  h2_muonIDEff_eta_phi_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_muonIDEffDen_eta_phi_PbPb18"));
  TH1D* h_muonIDEff_fcal_PbPb18 = (TH1D*) f_myEff->Get ("h_muonIDEffNum_fcal_PbPb18")->Clone ("h_muonIDEff_fcal_PbPb18");
  h_muonIDEff_fcal_PbPb18->Divide ((TH1D*) f_myEff->Get ("h_muonIDEffDen_fcal_PbPb18"));
  TF1* f_muonIDEff_fcal_PbPb18 = new TF1 ("f_muonIDeff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_muonIDEff_fcal_PbPb18->SetParameter (0, 1);
  f_muonIDEff_fcal_PbPb18->SetParameter (1, 0);
  h_muonIDEff_fcal_PbPb18->Fit (f_muonIDEff_fcal_PbPb18, "RN0Q");

  TH2D* h2_electronTrigEff_pt_eta_pp = (TH2D*) f_myEff->Get ("h2_electronTrigEffNum_pt_eta_pp")->Clone ("h2_electronTrigEff_pt_eta_pp");
  h2_electronTrigEff_pt_eta_pp->Divide ((TH2D*) f_myEff->Get ("h2_electronTrigEffDen_pt_eta_pp"));

  TH2D* h2_electronIDEff_pt_eta_pp = (TH2D*) f_myEff->Get ("h2_electronIDEffNum_pt_eta_pp")->Clone ("h2_electronIDEff_pt_eta_pp");
  h2_electronIDEff_pt_eta_pp->Divide ((TH2D*) f_myEff->Get ("h2_electronIDEffDen_pt_eta_pp"));

  TH2D* h2_electronTrigEff_pt_eta_PbPb18 = (TH2D*) f_myEff->Get ("h2_electronTrigEffNum_pt_eta_PbPb18")->Clone ("h2_electronTrigEff_pt_eta_PbPb18");
  h2_electronTrigEff_pt_eta_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_electronTrigEffDen_pt_eta_PbPb18"));
  TH1D* h_electronTrigEff_fcal_PbPb18 = (TH1D*) f_myEff->Get ("h_electronTrigEffNum_fcal_PbPb18")->Clone ("h_electronTrigEff_fcal_PbPb18");
  h_electronTrigEff_fcal_PbPb18->Divide ((TH1D*) f_myEff->Get ("h_electronTrigEffDen_fcal_PbPb18"));
  TF1* f_electronTrigEff_fcal_PbPb18 = new TF1 ("f_electronTrigeff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_electronTrigEff_fcal_PbPb18->SetParameter (0, 1);
  f_electronTrigEff_fcal_PbPb18->SetParameter (1, 0);
  h_electronTrigEff_fcal_PbPb18->Fit (f_electronTrigEff_fcal_PbPb18, "RN0Q");

  TH2D* h2_electronIDEff_pt_eta_PbPb18 = (TH2D*) f_myEff->Get ("h2_electronIDEffNum_pt_eta_PbPb18")->Clone ("h2_electronIDEff_pt_eta_PbPb18");
  h2_electronIDEff_pt_eta_PbPb18->Divide ((TH2D*) f_myEff->Get ("h2_electronIDEffDen_pt_eta_PbPb18"));
  TH1D* h_electronIDEff_fcal_PbPb18 = (TH1D*) f_myEff->Get ("h_electronIDEffNum_fcal_PbPb18")->Clone ("h_electronIDEff_fcal_PbPb18");
  h_electronIDEff_fcal_PbPb18->Divide ((TH1D*) f_myEff->Get ("h_electronIDEffDen_fcal_PbPb18"));
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
  if (!isMC) {
    if (!isPbPb) t->SetGetJets ();
    else if (isPbPb) t->SetGetZdc ();
  }
  t->SetBranchAddresses ();

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
    cout << "Error: In BkgEstimator.cxx: Invalid collision system, quitting." << endl;
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
  TFile* outFile;
  OutTree* outTree;
  const char* outTreeName = "PbPbZTree";
  outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");
  outFile->Delete (Form ("%s;*", outTreeName));
  outTree = new OutTree (outTreeName, outFile);
  outTree->SetBranchEventInfo ();
  outTree->SetBranchLeptons ();
  outTree->SetBranchZs ();
  if (!isPbPb) outTree->SetBranchJets ();
  outTree->SetBranchTracks ();
  outTree->SetBranches ();


  const int nEvts = tree->GetEntries ();
  TLorentzVector l1, l2;
  float trigEff1, trigEff2, recoEff1, recoEff2, trigEff, recoEff;

  // First loop over events looking for Z->ee candidates
  isEE = true;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In BkgEstimator.cxx: Electrons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;

    event_number = t->event_number;
    run_number = t->run_number;

    bool hasPrimaryVert = false;
    //bool hasPileupVert = false;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      if (t->vert_type[iVert] == 1) {
        if (hasPrimaryVert) {
          cout << "Info: In BkgEstimator.cxx: Found multiple primary vertices in this event" << endl;
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
    if (!hasPrimaryVert || fabs (vz) > 150)
      continue;

    event_weight = -1;
    if (!isMC) {
      if (electronTrigger->trigBool)
        event_weight = electronTrigger->trigPrescale;
    }
    else {
      //event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency / (mcNumberEvents);
      //event_weight = crossSectionPicoBarns * mcFilterEfficiency / (mcNumberEvents);
      event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency * 0.00171786; // sigma * f * L_int
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
    if (isPbPb) {
      iCent = GetFileCentBin (fcal_et);
      if (iCent < 1 || iCent > numFileCentBins-1)
        continue;
    }

    const float _event_weight = event_weight; // for storage, in case of multiple Z's in one event

    for (int iE1 = 0; iE1 < t->electron_n; iE1++) {
      l1_pt = t->electron_pt[iE1];
      l1_eta = t->electron_eta[iE1];
      l1_phi = t->electron_phi[iE1];
      l1_charge = t->electron_charge[iE1];
      l1_trk_pt = t->electron_id_track_pt[iE1];
      l1_trk_eta = t->electron_id_track_eta[iE1];
      l1_trk_phi = t->electron_id_track_phi[iE1];

      if (isPbPb) {
        l1_pt *= GetZmassSF_MC (fcal_et, l1_eta);
        if (!isMC)
          l1_pt *= GetZmassSF_PbPb (fcal_et, l1_eta);
      }

      if (l1_pt < electron_pt_cut)
        continue; // basic electron pT cut
      if (!InEMCal (l1_eta))
        continue; // reject electrons reconstructed outside the EMCal

      if (!isPbPb && !t->electron_lhloose[iE1])
        continue; // reject non-loose electrons
      else if (isPbPb && !t->electron_lhloose_hi[iE1])
        continue; // alternative for HI collisions

      if (fabs (t->electron_id_track_d0sig[iE1]) > 5)
        continue; // electron d0 significance vertex compatibility cut
      l1_d0sig = t->electron_id_track_d0sig[iE1];

      if (fabs ( (t->electron_id_track_z0[iE1] + t->electron_id_track_vz[iE1] - vz) * sin (t->electron_id_track_theta[iE1])) > 0.5)
        continue; // electron z0 vertex compatibility cut
      l1_z0 = t->electron_id_track_z0[iE1] + t->electron_id_track_vz[iE1]; // z0 in detector coordinates

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, electron_mass);

      for (int iE2 = 0; iE2 < iE1; iE2++) {
        l2_pt = t->electron_pt[iE2];
        l2_eta = t->electron_eta[iE2];
        l2_phi = t->electron_phi[iE2];
        l2_charge = t->electron_charge[iE2];
        l2_trk_pt = t->electron_id_track_pt[iE2];
        l2_trk_eta = t->electron_id_track_eta[iE2];
        l2_trk_phi = t->electron_id_track_phi[iE2];

        if (isPbPb) {
          l2_pt *= GetZmassSF_MC (fcal_et, l2_eta);
          if (!isMC)
            l2_pt *= GetZmassSF_PbPb (fcal_et, l2_eta);
        }

        if (l2_pt < electron_pt_cut)
          continue; // basic electron pT cut
        if (!InEMCal (l2_eta))
          continue; // reject electrons reconstructed outside the EMCal

        if (!isPbPb && !t->electron_lhloose[iE2])
          continue; // reject non-loose electrons
        else if (isPbPb && !t->electron_lhloose_hi[iE2])
          continue; // alternative for HI collisions

        if (fabs (t->electron_id_track_d0sig[iE2]) > 5)
          continue; // electron d0 significance vertex compatibility cut
        l2_d0sig = t->electron_id_track_d0sig[iE2];

        if (fabs ( (t->electron_id_track_z0[iE2] + t->electron_id_track_vz[iE2] - vz) * sin (t->electron_id_track_theta[iE2])) > 0.5)
          continue; // electron z0 vertex compatibility cut
        l2_z0 = t->electron_id_track_z0[iE2] + t->electron_id_track_vz[iE2]; // z0 in detector coordinates

        if (!isMC && !t->electron_matched[iE1] && !t->electron_matched[iE2])
          continue;

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, electron_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        isSameSign = (l1_charge == l2_charge);

        if (z_m < 60 || 150 < z_m)
          continue; // require Z to be in wider mass window for bkg. studies

    
        // for electrons, QCD bkg. corresponds to loose-non-medium and non-isolated electron pairs
        if (isPbPb)
          isQCDBkg = ((!t->electron_lhmedium_hi[iE1] && t->electron_etcone20[iE1]/l1_pt > 0.10*fcal_et*1e-3+0.050) || (!t->electron_lhmedium_hi[iE2] && t->electron_etcone20[iE2]/l2_pt > 0.10*fcal_et*1e-3+0.050));
          //isQCDBkg = (!t->electron_lhmedium_hi[iE1] && !t->electron_lhmedium_hi[iE2] && fmin (t->electron_etcone20[iE1]/l1_pt, t->electron_etcone20[iE2]/l2_pt) > 0.10*fcal_et*1e-3 +0.050);
        else
          isQCDBkg = ((!t->electron_lhmedium[iE1] && t->electron_etcone20[iE1]/l1_pt > 0.050) || (!t->electron_lhmedium[iE2] && t->electron_etcone20[iE2]/l2_pt > 0.050));
          //isQCDBkg = (!t->electron_lhmedium[iE1] && !t->electron_lhmedium[iE2] && fmin (t->electron_etcone20[iE1]/l1_pt, t->electron_etcone20[iE2]/l2_pt) > 0.050);
    

        // find correct Z boson weighting factor
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
            if (z_pt < 25 && fabs (z_y) < 1.75)       event_weight = event_weight / (415887+638346+638321+976912);
            else if (z_pt >= 25 && fabs (z_y) < 1.75) event_weight = event_weight / (415887+638346+638321+976912+415862+638366+638366+976816);
            else if (z_pt < 25 && fabs (z_y) >= 1.75) event_weight = event_weight / (415887+638346+638321+976912+415862+638366+637401+976816);
            else {
              assert (z_pt >= 25 && fabs (z_y) >= 1.75);
              event_weight = event_weight / (415887+638346+638321+976912+415862+638366+638366+976816+415862+638366+637401+976816);
            }
          }
          else event_weight = event_weight / 1000000;
        }

        outTree->Fill ();
        
      } // end iE2 loop
    } // end iE1 loop

  } // end electron selection
  cout << endl << "Info: In BkgEstimator.cxx: Finished processing electrons." << endl;


  // Now loop over events looking for Z->mumu candidates
  isEE = false;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In BkgEstimator.cxx: Muons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    event_number = t->event_number;
    run_number = t->run_number;

    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;
    if ((!isMC || isOverlayMC) && !t->passes_toroid)
      continue; // additional check for muon quality

    bool hasPrimaryVert = false;
    //bool hasPileupVert = false;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      if (t->vert_type[iVert] == 1) {
        if (hasPrimaryVert) {
          cout << "Info: In BkgEstimator.cxx: Found multiple primary vertices in this event" << endl;
          hasPrimaryVert = false;
          break;
        }
        hasPrimaryVert = (t->vert_type[iVert] == 1);
        vz = t->vert_z[iVert];
      }
      //else if (t->vert_type[iVert] == 3)
      //  hasPileupVert = true;
    }
    if (!hasPrimaryVert || fabs (vz) > 150)
      continue; // check for primary vertex

    event_weight = -1;
    if (!isMC) {
      if (muonTrigger->trigBool)
        event_weight = muonTrigger->trigPrescale;
    }
    else {
      //event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency / (mcNumberEvents);
      //event_weight = crossSectionPicoBarns * mcFilterEfficiency / (mcNumberEvents);
      event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency * 0.00143844; // sigma * f * L_int
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
    if (isPbPb) {
      iCent = GetFileCentBin (fcal_et);
      if (iCent < 1 || iCent > numFileCentBins-1)
        continue;
    }

    const float _event_weight = event_weight; // for storage, in case of multiple Z's in one event

    for (int iM1 = 0; iM1 < t->muon_n; iM1++) {
      if (!isPbPb)
        l1_pt = t->muon_pt[iM1];
      else
        l1_pt = t->muon_ms_pt[iM1];
      l1_eta = t->muon_eta[iM1];
      l1_phi = t->muon_phi[iM1];
      l1_charge = t->muon_charge[iM1];
      l1_trk_pt = t->muon_id_track_pt[iM1];
      l1_trk_eta = t->muon_id_track_eta[iM1];
      l1_trk_phi = t->muon_id_track_phi[iM1];

      if (l1_pt < muon_pt_cut)
        continue; // basic muon pT cut
      if (fabs (l1_eta) > 2.5)
        continue; // reject muons reconstructed outside the detector

      if (!t->muon_medium[iM1])
        continue; // reject non-medium muons

      if (fabs (t->muon_id_track_d0sig[iM1]) > 3)
        continue; // muon d0 cut
      l1_d0sig = t->muon_id_track_d0sig[iM1];

      if (fabs ( (t->muon_id_track_z0[iM1] + t->muon_id_track_vz[iM1] - vz) * sin (t->muon_id_track_theta[iM1])) > 0.5)
        continue; // muon z0 to vertex cut
      l1_z0 = t->muon_id_track_z0[iM1] + t->muon_id_track_vz[iM1]; // z0 in detector coordinates

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, muon_mass);

      for (int iM2 = 0; iM2 < iM1; iM2++) {
        if (!isPbPb)
          l2_pt = t->muon_pt[iM2];
        else 
          l2_pt = t->muon_ms_pt[iM2];
        l2_eta = t->muon_eta[iM2];
        l2_phi = t->muon_phi[iM2];
        l2_charge = t->muon_charge[iM2];
        l2_trk_pt = t->muon_id_track_pt[iM2];
        l2_trk_eta = t->muon_id_track_eta[iM2];
        l2_trk_phi = t->muon_id_track_phi[iM2];

        if (l2_pt < muon_pt_cut)
          continue; // basic muon pT cut
        if (fabs (l2_eta) > 2.5)
          continue; // reject muons reconstructed outside the detector

        if (!t->muon_medium[iM2])
          continue; // reject non-medium muons

        if (fabs (t->muon_id_track_d0sig[iM2]) > 3)
          continue; // muon d0 significance vertex compatibility cut
        l2_d0sig = t->muon_id_track_d0sig[iM2];

        if (fabs ( (t->muon_id_track_z0[iM2] + t->muon_id_track_vz[iM2] - vz) * sin (t->muon_id_track_theta[iM2])) > 0.5)
          continue; // muon z0 vertex compatibility cut
        l2_z0 = t->muon_id_track_z0[iM2] + t->muon_id_track_vz[iM2]; // z0 in detector coordinates

        if (!isMC && !t->muon_matched[iM1] && !t->muon_matched[iM2])
          continue;

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, muon_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        isSameSign = (l1_charge == l2_charge);

        if (z_m < 60 || 150 < z_m)
          continue; // require Z to be in mass window for bkg. studies


        isQCDBkg = (isSameSign); // for muons, QCD bkg. corresponds to same-sign pairs since charge mis-identification rate is very low


        // find correct Z boson weighting factor
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
            if (z_pt < 25 && fabs (z_y) < 1.75)       event_weight = event_weight / (357861+541661+561006+840470);
            else if (z_pt >= 25 && fabs (z_y) < 1.75) event_weight = event_weight / (357861+541661+561006+840470+357894+560953+560953+841487);
            else if (z_pt < 25 && fabs (z_y) >= 1.75) event_weight = event_weight / (357861+541661+561006+840470+357900+560051+561091+839951);
            else {
              assert (z_pt >= 25 && fabs (z_y) >= 1.75);
              event_weight = event_weight / (357861+541661+561006+840470+357894+560953+560953+841487+357900+560051+561091+839951);
            }
          }
          else event_weight = event_weight / 3360000;
        }

        outTree->Fill ();
        
      } // end iM2 loop
    } // end iM1 loop

  } // end muon selection
  cout << endl << "Info: In BkgEstimator.cxx: Finished processing muons." << endl;

  SaferDelete (&electronTrigger);
  SaferDelete (&muonTrigger);

  SaferDelete (&t);

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (&outFile);

  SaferDelete (&h_zdcCuts);

  return true;
}

} // end namespace

#endif
